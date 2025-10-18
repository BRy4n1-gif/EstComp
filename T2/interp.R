# Cargar librerías
library(nortest)
library(car)

# ================================
# DATOS
# ================================
ph_agua <- c(8.2, 8.5, 7.9, 8.3, 7.8, 8.6, 8.1, 8.4, 7.7, 8.7,
             8.3, 8.0, 8.8, 7.6, 8.5, 8.2, 8.1, 8.6, 7.9, 8.4,
             8.5, 8.0, 8.3, 8.2, 8.1)

oxigeno_disuelto <- c(7.2, 8.1, 6.5, 7.8, 6.2, 8.5, 7.0, 7.9, 6.8, 8.3,
                      7.5, 7.1, 8.7, 6.4, 8.0, 7.3, 7.2, 8.4, 6.9, 7.8,
                      8.2, 7.0, 7.6, 7.4, 7.3)

ingresos_mensuales <- c(850, 920, 750, 1100, 890, 980, 820, 1050, 870, 940, 
                        960, 810, 1080, 780, 1020, 900, 950, 830, 1120, 800,
                        1060, 880, 930, 820, 1140, 770, 1000, 890, 960, 840,
                        1090, 790, 1030, 910, 940, 850, 1110, 760, 1010, 920,
                        980, 860, 1070, 810, 1040, 930, 890, 1000, 820, 1130,
                        780, 960, 900, 1050, 840, 1080, 790, 980, 920, 950,
                        870, 1100, 750, 1010, 880, 960, 830, 1090, 800, 1040,
                        930, 910, 980, 840, 1120, 790, 1020, 890, 960, 850,
                        1080, 820, 1010, 920, 950, 860, 1110, 770, 980, 900,
                        1060, 840, 1090, 810, 960, 930, 880, 1050, 820, 1130)

captura_semanal_kg <- c(45, 52, 38, 58, 42, 55, 40, 56, 44, 53, 51, 39, 59, 36, 54,
                        46, 50, 41, 60, 37, 57, 43, 52, 40, 61, 35, 54, 45, 50, 42,
                        58, 38, 55, 47, 51, 43, 59, 34, 53, 48, 55, 44, 57, 40, 56,
                        50, 45, 54, 41, 60, 37, 51, 46, 56, 42, 58, 38, 52, 47, 50,
                        44, 59, 35, 53, 45, 50, 41, 57, 39, 56, 50, 46, 54, 42, 60,
                        38, 55, 45, 50, 43, 58, 40, 53, 47, 51, 44, 59, 36, 52, 46,
                        57, 42, 58, 40, 50, 49, 45, 56, 41, 61)

tipo_embarcacion <- factor(c(rep("Traditional", 40), rep("Motor", 35), rep("Vela", 25)))

datos_pescadores <- data.frame(
  ingresos = ingresos_mensuales,
  captura = captura_semanal_kg,
  embarcacion = tipo_embarcacion
)

# ================================
# RESULTADOS PRINCIPALES
# ================================

# 1. Comparación de ingresos entre Traditional vs Motor
datos_comparacion <- datos_pescadores[datos_pescadores$embarcacion %in% c("Traditional", "Motor"), ]
datos_comparacion$embarcacion <- droplevels(datos_comparacion$embarcacion)

test_ingresos <- wilcox.test(ingresos ~ embarcacion, data = datos_comparacion)

# 2. Correlación ingresos vs captura
spearman_cor <- cor.test(datos_pescadores$ingresos, datos_pescadores$captura, method = "spearman")

# 3. Oxígeno vs referencia
wilcox_oxigeno <- wilcox.test(oxigeno_disuelto, mu = 7.5)

# ================================
# MOSTRAR SOLO RESULTADOS
# ================================
cat("COMPARACIÓN INGRESOS (Traditional vs Motor):\n")
cat("Prueba:", class(test_ingresos)[1], "\n")
cat("p-valor:", round(test_ingresos$p.value, 4), "\n\n")

cat("CORRELACIÓN INGRESOS-CAPTURA:\n")
cat("Método: Spearman\n")
cat("rho:", round(spearman_cor$estimate, 4), "\n")
cat("p-valor:", round(spearman_cor$p.value, 4), "\n\n")

cat("OXÍGENO DISUELTO vs 7.5 mg/L:\n")
cat("Prueba: Wilcoxon una muestra\n")
cat("p-valor:", round(wilcox_oxigeno$p.value, 4), "\n")
