# ==============================================================================
# ANÁLISIS MULTIVARIADO CON INTERVALOS MÉTRICOS DE 1 CM (MÁXIMA RESOLUCIÓN)
# ==============================================================================

library(tidyverse)
library(vegan)
library(permute)

# ------------------------------------------------------------------------------
# 1. CARGA DE DATOS
# ------------------------------------------------------------------------------
datos_originales <- read_csv("BASE DE DATOS PLAYA - 18_05_2026.csv") %>% 
  rename_all(str_trim) %>% 
  mutate(
    ZONA = str_trim(ZONA),
    `M/T` = str_trim(`M/T`),
    FECHA = str_trim(FECHA),
    NUCLEO = str_trim(NUCLEO),
    LONGITUD = as.numeric(LONGITUD)
  )

# ------------------------------------------------------------------------------
# 2. TRANSFORMACIÓN A MATRIZ DE INTERVALOS CONTINUOS DE 1 CM
# ------------------------------------------------------------------------------

# Determinamos el valor máximo para cerrar el último intervalo de forma exacta
max_longitud <- max(datos_originales$LONGITUD, na.rm = TRUE)
limite_superior <- ceiling(max_longitud) # Redondea 

datos_por_intervalos <- datos_originales %>% 
  mutate(
    # Creamos los intervalos de 1 en 1 cm. El argumento 'right = FALSE' hace que sean [0,1), [1,2), etc.
    Intervalo_cm = cut(LONGITUD, 
                       breaks = seq(0, limite_superior, by = 1), 
                       right = FALSE,
                       labels = paste0(seq(0, limite_superior - 1, by = 1), "-", seq(1, limite_superior, by = 1)))
  )

# C. Construimos la matriz de comunidad donde CADA COLUMNA ES UN INTERVALO DE 1 CM
matriz_comunidad_intervalos <- datos_por_intervalos %>% 
  group_by(FECHA, ZONA, `M/T`, NUCLEO, Intervalo_cm) %>% 
  summarise(Abundancia = n(), .groups = 'drop') %>% 
  pivot_wider(names_from = Intervalo_cm, values_from = Abundancia, values_fill = 0)

# Aseguramos el orden correcto de las columnas de menor a mayor intervalo
columnas_ordenadas <- sort(colnames(matriz_comunidad_intervalos %>% select(-FECHA, -ZONA, -`M/T`, -NUCLEO)))

# Separamos metadatos de la nueva matriz numérica por centímetros
metadatos <- matriz_comunidad_intervalos %>% select(FECHA, ZONA, `M/T`, NUCLEO)
matriz_numerica_cm <- matriz_comunidad_intervalos %>% select(all_of(columnas_ordenadas))

factor_combinado <- paste(metadatos$ZONA, metadatos$`M/T`, sep = "_")

# Ver en la consola la estructura de las columnas
cat("\n=== MATRIZ MULTIVARIADA POR INTERVALOS DE 1 CM ===\n")
print(head(matriz_numerica_cm, 5))

# ------------------------------------------------------------------------------
# 3. EJECUCIÓN DE LOS ANÁLISIS
# ------------------------------------------------------------------------------
# NMDS con los intervalos métricos
set.seed(123)
nmds_intervalos <- metaMDS(matriz_numerica_cm, distance = "bray", k = 2, trymax = 100)

# PERMANOVA2 con los intervalos métricos
cat("\n=== RESULTADOS: PERMANOVA2 (Estructura por Intervalos de 1 cm) ===\n")
resultado_custom_intervalos <- PERMANOVA2(x = matriz_numerica_cm, factor = factor_combinado, distancia = "bray", nperm = 999)
print(resultado_custom_intervalos)

# ------------------------------------------------------------------------------
# 4. GRAFICOS 
# ------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------
#FIGURA 1: Ordenación bidimensional mediante Escalamiento Multidimensional No Métrico (NMDS)
#--------------------------------------------------------------------------------------------

nmds_coords <- as.data.frame(scores(nmds_intervalos, display = "sites")) %>% bind_cols(metadatos)

ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = `M/T`, shape = ZONA)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = interaction(ZONA, `M/T`)), linetype = 2, linewidth = 0.5) +
  labs(
    title = "NMDS: Estructura de Población por Intervalos Métricos de 1 cm",
    subtitle = paste("Distancia Bray-Curtis | Stress =", round(nmds_intervalos$stress, 4)),
    x = "NMDS Dimensión 1",
    y = "NMDS Dimensión 2",
    color = "Horario (M/T)",
    shape = "Zona de Muestreo"
  ) +
  scale_color_manual(values = c("M" = "#1f77b4", "T" = "#ff7f0e"), labels = c("Mañana", "Tarde")) +
  theme_bw(base_size = 12)

#---------------------------------------------------------------------
#FIGURA 2: Histograma de frecuencias por escenarios espaciotemporales
#---------------------------------------------------------------------
  
# 1. Preparar los datos en formato largo para ggplot
datos_histograma <- datos_por_intervalos %>% 
  filter(!is.na(Intervalo_cm)) %>% 
  group_by(ZONA, `M/T`, Intervalo_cm) %>% 
  summarise(Abundancia_Total = n(), .groups = 'drop')

# 2. Crear el gráfico con barras interactivas por escenario
grafico_histograma <- ggplot(datos_histograma, aes(x = Intervalo_cm, y = Abundancia_Total, fill = `M/T`)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.8) +
  # Dividir el gráfico en páneles por Zona y Horario
  facet_grid(ZONA ~ `M/T`, scales = "free_y", 
             labeller = labeller(`M/T` = c("M" = "Mañana", "T" = "Tarde"))) +
  labs(
    title = "Distribución de Frecuencias de Tallas de Oligoquetos",
    subtitle = "Estructura poblacional por intervalos métricos de 1 cm en Dzilam de Bravo",
    x = "Intervalos de Longitud Corporal (cm)",
    y = "Abundancia Total",
    fill = "Horario"
  ) +
  # Colores: azul para mañana y  naranja para tarde
  scale_color_manual(values = c("M" = "#1f77b4", "T" = "#ff7f0e")) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1), # Gira las etiquetas de los cm para que no se amontonen
    strip.text = element_text(face = "bold", size = 11), # Pone los títulos de los cuadros en negrita
    legend.position = "none" # Quitamos la leyenda porque los paneles ya dicen Mañana y Tarde
  )

# 3. Ver grafico
print(grafico_histograma)

#---------------------------------------------------------------------
#FIGURA 3: Diagrama de cajas (BOXPLOT) por intervalos clave
#---------------------------------------------------------------------

# 1. Identificar automáticamente cuáles son los 4 intervalos con mayor abundancia total
intervalos_top4 <- datos_por_intervalos %>% 
  filter(!is.na(Intervalo_cm)) %>% 
  group_by(Intervalo_cm) %>% 
  summarise(Total = n()) %>% 
  arrange(desc(Total)) %>% 
  slice(1:4) %>% 
  pull(Intervalo_cm)

# 2. Filtrar la matriz original para quedarnos solo con esos 4 intervalos top
datos_boxplot <- datos_por_intervalos %>% 
  filter(Intervalo_cm %in% intervalos_top4) %>% 
  # Contamos cuántos organismos de cada talla cayeron en cada núcleo específico
  group_by(ZONA, `M/T`, NUCLEO, Intervalo_cm) %>% 
  summarise(Abundancia_Nucleo = n(), .groups = 'drop')

# 3. Crear el gráfico de cajas
grafico_boxplot <- ggplot(datos_boxplot, aes(x = `M/T`, y = Abundancia_Nucleo, fill = `M/T`)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2, alpha = 0.7) +
  # geom_jitter añade los puntos de cada núcleo para ver cómo se distribuyen los datos reales
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 1.2) +
  # Dividir por Zona y por Intervalo de Talla
  facet_grid(ZONA ~ Intervalo_cm, scales = "free_y") +
  labs(
    title = "Variabilidad de Abundancia en Intervalos Dominantes",
    subtitle = "Comparativa de núcleos individuales para las cuatro clases métricas más abundantes",
    x = "Horario de Muestreo (M: Mañana / T: Tarde)",
    y = "Abundancia por Núcleo (unidades)",
    fill = "Horario"
  ) +
  scale_fill_manual(values = c("M" = "#1f77b4", "T" = "#ff7f0e"), labels = c("Mañana", "Tarde")) +
  scale_x_discrete(labels = c("M" = "Mañana", "T" = "Tarde")) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  )

# 4. Ver grafico
print(grafico_boxplot)  
  
  
