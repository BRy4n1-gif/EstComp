# escalado.R
# Aplicación de técnicas de escalado a un conjunto de datos multivariado

# Cargar librerías
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("gridExtra")) install.packages("gridExtra")
library(dplyr)
library(ggplot2)
library(gridExtra)

# Generar datos de ejemplo: variables con distintas escalas
set.seed(123)
n <- 500
datos <- data.frame(
  ingresos = rlnorm(n, meanlog = 10, sdlog = 0.8),   # Sesgado positivo
  edad = rnorm(n, mean = 40, sd = 12),                # Aprox. normal
  pH_agua = runif(n, min = 6.5, max = 9.0),           # Uniforme
  captura_kg = rpois(n, lambda = 50)                  # Discreta, sesgada
)

# Eliminar valores no válidos
datos <- datos %>% filter(if_all(everything(), ~ !is.na(.x) & is.finite(.x)))

# Funciones de escalado
escalar_normalizar <- function(x) (x - min(x)) / (max(x) - min(x))
escalar_estandarizar <- function(x) (x - mean(x)) / sd(x)
escalar_robusto <- function(x) (x - median(x)) / IQR(x)
escalar_log <- function(x) {
  if (min(x, na.rm = TRUE) <= 0) {
    x <- x - min(x, na.rm = TRUE) + 1e-6
  }
  log(x)
}

# Aplicar transformaciones
datos_esc <- datos %>%
  mutate(
    ingresos_log = escalar_log(ingresos),
    ingresos_norm = escalar_normalizar(ingresos),
    ingresos_std = escalar_estandarizar(ingresos),
    ingresos_rob = escalar_robusto(ingresos),
    edad_std = escalar_estandarizar(edad),
    pH_norm = escalar_normalizar(pH_agua),
    captura_std = escalar_estandarizar(captura_kg)
  )

# Guardar datos escalados
write.csv(datos_esc, "datos_escalados.csv", row.names = FALSE)

# Función para graficar distribuciones
graficar_distribucion <- function(var, titulo) {
  ggplot(data.frame(x = var), aes(x = x)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    labs(title = titulo, x = "Valor", y = "Densidad") +
    theme_minimal()
}

# Generar gráficos
p1 <- graficar_distribucion(datos$ingresos, "Ingresos (original)")
p2 <- graficar_distribucion(datos_esc$ingresos_log, "Ingresos (log)")
p3 <- graficar_distribucion(datos_esc$ingresos_norm, "Ingresos (normalizado)")
p4 <- graficar_distribucion(datos_esc$ingresos_std, "Ingresos (estandarizado)")

# Combinar gráficos
g_ingresos <- grid.arrange(p1, p2, p3, p4, ncol = 2, top = "Transformaciones de Ingresos")

# Guardar gráfico
ggsave("escalado_ingresos.png", plot = g_ingresos, width = 10, height = 8)

# Resumen estadístico
resumen <- datos_esc %>%
  select(ingresos, ingresos_log, ingresos_norm, ingresos_std, ingresos_rob) %>%
  summarise_all(list(media = ~mean(.), sd = ~sd(.), mediana = ~median(.), iqr = ~IQR(.)))

write.csv(resumen, "resumen_escalado.csv", row.names = FALSE)
