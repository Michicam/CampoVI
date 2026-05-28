# ==============================================================================
# ANÁLISIS MULTIVARIADO DE DOS VÍAS: VERSIÓN ULTRA-ROBUSTA
# ==============================================================================

# 1. Cargar librerías necesarias
library(tidyverse)
library(vegan)
library(permute)

# ------------------------------------------------------------------------------
# 2. DEFINICIÓN DE FUNCIONES OPTIMIZADAS (PERMANOVA 2)
# ------------------------------------------------------------------------------
SS <- function(d) {
  n <- dim(as.matrix(d))[1]
  ss <- sum(d ^ 2) / n
  return(ss)
}

v <- function(d) {
  n <- dim(as.matrix(d))[1]
  ss <- sum(d ^ 2) / n
  v <- ss / (n - 1)
  return(v)
}

pseudo.F.optimizada <- function(x, factor, distancia) {
  x_mat <- as.data.frame(x)
  d <- vegdist(x_mat, method = distancia)
  TSS <- SS(d)
  
  factor <- as.factor(factor)
  lab <- levels(factor)
  lev <- table(factor)
  
  CR <- numeric(length = nlevels(factor))
  for (i in 1:nlevels(factor)) {
    sub_grupo <- x_mat[factor == lab[i], , drop = FALSE]
    CR[i] <- SS(vegdist(sub_grupo, method = distancia))
  }
  RSS <- sum(CR)
  
  Var <- numeric(length = nlevels(factor))
  d.res <- as.data.frame(matrix(nrow = nlevels(factor), ncol = 3))
  colnames(d.res) <- c("n", "Var", "Var_ponderada")
  
  for (i in 1:nlevels(factor)) {
    sub_grupo <- x_mat[factor == lab[i], , drop = FALSE]
    Var[i] <- v(vegdist(sub_grupo, method = distancia))
    d.res[i, ] <- c(lev[i], Var[i], (1 - (lev[i] / sum(lev))) * Var[i])
  }
  
  den <- sum(d.res$Var_ponderada)
  ASS <- TSS - RSS
  Fobs <- ASS / den
  return(Fobs)
}

PERMANOVA2 <- function(x, factor, distancia, nperm = 999) {
  control <- how(nperm = nperm, within = Within(type = "free"))
  Fobs <- pseudo.F.optimizada(x, factor, distancia = distancia)
  Nobs <- nobs(x)
  
  F.permu <- numeric(length = control$nperm + 1)
  F.permu[1] <- Fobs
  
  for (i in 1:control$nperm) {
    want <- permute(i, Nobs, control)
    F.permu[i + 1] <- pseudo.F.optimizada(x[want, , drop = FALSE], factor, distancia = distancia)
  }
  
  pval <- sum(abs(F.permu) >= abs(F.permu[1])) / (control$nperm + 1)
  return(data.frame("Pseudo.F" = F.permu[1], "p.perm" = pval))
}

# ------------------------------------------------------------------------------
# 3. CARGA Y LIMPIEZA DE LA BASE DE DATOS (Con seguro anti-espacios en columnas)
# ------------------------------------------------------------------------------

datos_playa <- read_csv("BASE DE DATOS PLAYA - 18_05_2026.csv")

# ¡PASO CLAVE! Limpiamos espacios tanto en las filas como en los nombres de columnas
colnames(datos_playa) <- str_trim(colnames(datos_playa))

datos_limpios <- datos_playa %>% 
  mutate(
    ZONA = str_trim(ZONA),
    NUCLEO = str_trim(NUCLEO),
    `M/T` = str_trim(`M/T`),
    FECHA = str_trim(FECHA),
    LONGITUD = as.numeric(LONGITUD) # Nos aseguramos de que R la reconozca como numérica
  )

# ------------------------------------------------------------------------------
# 4. PREPARACIÓN ESTRICTA DE MATRICES (Evita errores de entorno en adonis2)
# ------------------------------------------------------------------------------

# Creamos los vectores y matrices de forma independiente para que adonis2 no se confunda
matriz_respuesta <- datos_limpios %>% 
  select(LONGITUD) %>% 
  as.data.frame()

factor_zona <- as.factor(datos_limpios$ZONA)
factor_hora <- as.factor(datos_limpios$`M/T`)

# Factor combinado para tu función PERMANOVA2
factor_combinado <- datos_limpios %>% 
  mutate(ZONA_HORA = paste(ZONA, `M/T`, sep = "_")) %>% 
  pull(ZONA_HORA)

# ------------------------------------------------------------------------------
# 5. EJECUCIÓN DE LOS ANÁLISIS
# ------------------------------------------------------------------------------

cat("\n=== RUNNING: PERMANOVA CLÁSICO (adonis2) ===\n")
# Al pasarle la matriz y los vectores directamente mediante el signo $, 
# R ya no tiene que buscar en ambientes flotantes.
resultado_permanova_classic <- adonis2(
  matriz_respuesta ~ factor_zona * factor_hora, 
  method = "euclidean", 
  permutations = 999
)
print(resultado_permanova_classic)


cat("\n=== RUNNING: PERMANOVA2 MODIFICADO (Anderson et al., 2017) ===\n")
resultado_permanova_custom <- PERMANOVA2(
  x = matriz_respuesta, 
  factor = factor_combinado, 
  distancia = "euclidean", 
  nperm = 999
)
print(resultado_permanova_custom)



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


# ==============================================================================
# VISUALIZACIONES EN R: RESPONDIDENDO A LA HIPÓTESIS DE DZILAM DE BRAVO
# ==============================================================================

# 1. Cargar librerías
library(tidyverse)


# 2. Cargar y limpiar datos (Asegúrate de que el nombre del CSV sea el correcto)
datos_playa <- read_csv("BASE DE DATOS PLAYA - 18_05_2026.csv") %>% 
  mutate(
    ZONA = str_trim(ZONA),
    `M/T` = str_trim(`M/T`)
  )

# ------------------------------------------------------------------------------
# GRÁFICO 1: COMPARATIVA DE ABUNDANCIA (¿Hay más organismos en la mañana?)
# Este gráfico calcula el número total de oligoquetos encontrados por combinación.
# ------------------------------------------------------------------------------
abundancia_zona_hora <- datos_playa %>% 
  group_by(ZONA, `M/T`) %>% 
  summarise(Abundancia_Total = n(), .groups = 'drop') # Cuenta cuántos organismos hay

grafico_abundancia <- ggplot(abundancia_zona_hora, aes(x = ZONA, y = Abundancia_Total, fill = `M/T`)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black") +
  geom_text(aes(label = Abundancia_Total), position = position_dodge(0.8), vjust = -0.5, fontface = "bold") +
  labs(
    title = "Abundancia Total de Oligoquetos Bentónicos en Dzilam de Bravo",
    subtitle = "Comparativa temporal para evaluar la hipótesis de migración",
    x = "Zona de Muestreo (Z1: Supralitoral | Z2: Barrido/Mesolitoral)",
    y = "Número Total de Organismos Recolectados (Abundancia)",
    fill = "Horario de Muestreo"
  ) +
  scale_fill_manual(values = c("M" = "#1f77b4", "T" = "#ff7f0e"), labels = c("M" = "Mañana (M)", "T" = "Tarde (T)")) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(grafico_abundancia)

# ------------------------------------------------------------------------------
# GRÁFICO 2: ESTRUCTURA DE TALLAS (Densidad de frecuencias superpuestas)
# Este gráfico responde a cómo se distribuyen los tamaños de los organismos.
# ------------------------------------------------------------------------------
grafico_tallas <- ggplot(datos_playa, aes(x = LONGITUD, fill = `M/T`)) +
  geom_density(alpha = 0.5, color = "black") +
  facet_wrap(~ZONA, labeller = as_labeller(c("Z1" = "Zona Supralitoral (Z1)", "Z2" = "Zona de Barrido (Z2)"))) +
  labs(
    title = "Estructura de Tallas de Oligoquetos según el Horario",
    subtitle = "Análisis de distribución de densidad de longitudes",
    x = "Longitud de los organismos (mm)",
    y = "Densidad de Frecuencia",
    fill = "Horario de Muestreo"
  ) +
  scale_fill_manual(values = c("M" = "#1f77b4", "T" = "#ff7f0e"), labels = c("M" = "Mañana (M)", "T" = "Tarde (T)")) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey20"),
    strip.text = element_text(color = "white", face = "bold"),
    legend.position = "bottom"
  )

print(grafico_tallas)
