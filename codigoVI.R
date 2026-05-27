
# ==============================================================================
# ANÁLISIS MULTIVARIADO: PERMANOVA DE DOS VÍAS (ZONA Y HORA) - TODOS LOS DÍAS
# ==============================================================================

# 1. Cargar librerías necesarias
library(tidyverse)
library(vegan)       # Paquete esencial para ecología (contiene adonis2)

# 2. Cargar y limpiar la base de datos (Inlcuye el 22 de mayo)
datos_playa <- read_csv("BASE DE DATOS PLAYA - 18_05_2026.csv")

datos_limpios <- datos_playa %>% 
  mutate(
    ZONA = str_trim(ZONA),
    NUCLEO = str_trim(NUCLEO),
    `M/T` = str_trim(`M/T`),
    FECHA = str_trim(FECHA)
  )

# ------------------------------------------------------------------------------
# 3. PREPARACIÓN DE LAS MATRICES
# Separamos la variable de respuesta de tus factores de diseño.
# ------------------------------------------------------------------------------

# Matriz de la variable biológica (Longitud de los oligoquetos)
matriz_variables <- datos_limpios %>% 
  select(LONGITUD)

# Matriz de diseño/factores ambientales
matriz_factores <- datos_limpios %>% 
  select(ZONA, `M/T`) %>% 
  mutate(
    ZONA = as.factor(ZONA),
    `M/T` = as.factor(`M/T`)
  )

# ------------------------------------------------------------------------------
# 4. EJECUCIÓN DEL PERMANOVA (adonis2)
# ------------------------------------------------------------------------------

# Al ser una sola variable continua (LONGITUD), se utiliza "euclidean" 
# como método de distancia para el análisis de permutaciones.
resultado_permanova <- adonis2(
  matriz_variables ~ ZONA * `M/T`,  # Factores principales y su interacción
  data = matriz_factores, 
  method = "euclidean",             # Distancia idónea para medidas lineales continuas
  permutations = 999                # Número de permutaciones libres
)

# Mostrar la tabla de resultados en tu reporte de Quarto
print(resultado_permanova)


# ==============================================================================
# VISUALIZACIÓN DE DATOS: BOXPLOT E INTERACIÓN (TODOS LOS DÍAS)
# ==============================================================================

# 1. Cargar librerías
library(tidyverse)

# 2. Cargar y limpiar datos
datos_playa <- read_csv("BASE DE DATOS PLAYA - 18_05_2026.csv") %>% 
  mutate(
    ZONA = str_trim(ZONA),
    `M/T` = str_trim(`M/T`)
  )

# ------------------------------------------------------------------------------
# GRÁFICO 1: Boxplot de la Longitud por Zona y Horario
# Ideal para ver la dispersión, la mediana y los datos atípicos
# ------------------------------------------------------------------------------
grafico_boxplot <- ggplot(datos_playa, aes(x = ZONA, y = LONGITUD, fill = `M/T`)) +
  geom_boxplot(outlier.color = "red", alpha = 0.7, width = 0.6) +
  # Añadir los puntos de los datos reales con un poco de dispersión (jitter)
  geom_jitter(aes(group = `M/T`), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6), 
              alpha = 0.3, size = 1, color = "black") +
  labs(
    title = "Distribución de la Longitud de Oligoquetos",
    subtitle = "Comparativa por Zona Intermareal y Horario (Muestreo Completo)",
    x = "Zona de Muestreo (Z1: Supralitoral | Z2: Mesolitoral)",
    y = "Longitud de los organismos (mm)",
    fill = "Horario"
  ) +
  scale_fill_manual(values = c("M" = "#4ea8de", "T" = "#f77f00"), labels = c("M" = "Mañana", "T" = "Tarde")) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(grafico_boxplot)

# ------------------------------------------------------------------------------
# GRÁFICO 2: Gráfico de Interacción (Líneas de Tendencia)
# ¡Este gráfico explica directamente el resultado del PERMANOVA!
# ------------------------------------------------------------------------------
# Calcular los promedios primero
promedios <- datos_playa %>%
  group_by(ZONA, `M/T`) %>%
  summarise(Longitud_Promedio = mean(LONGITUD, na.rm = TRUE), .groups = 'drop')

grafico_interaccion <- ggplot(promedios, aes(x = ZONA, y = Longitud_Promedio, color = `M/T`, group = `M/T`)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  labs(
    title = "Gráfico de Interacción: Zona x Horario",
    subtitle = "Tendencia de la longitud promedio",
    x = "Zona de Muestreo",
    y = "Longitud Promedio (mm)",
    color = "Horario"
  ) +
  scale_color_manual(values = c("M" = "#4ea8de", "T" = "#f77f00"), labels = c("M" = "Mañana", "T" = "Tarde")) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

print(grafico_interaccion)


