# Cargar librerías necesarias
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)

# Función para simular una ejecución del modelo
simular_una_vez <- function(total_dulces = 30, objetivo_chupetines = 10, max_iter = 1000) {
  # Generar distribución aleatoria de 30 dulces entre x, y, z
  # Usamos rmultinom para una partición aleatoria uniforme
  dulces_iniciales <- rmultinom(1, size = total_dulces, prob = c(1/3, 1/3, 1/3))
  x <- dulces_iniciales[1]
  y <- dulces_iniciales[2]
  z <- dulces_iniciales[3]
  c <- 0  # chupetines iniciales
  chupetines_producidos <- 0
  interacciones <- 0
  
  while (chupetines_producidos < objetivo_chupetines && interacciones < max_iter) {
    # Regla 1: Si tenemos al menos uno de cada dulce, producimos un chupetín
    if (x > 0 && y > 0 && z > 0) {
      x <- x - 1
      y <- y - 1
      z <- z - 1
      c <- c + 1
      chupetines_producidos <- chupetines_producidos + 1
      interacciones <- interacciones + 1
    }
    # Regla 2: Si no podemos producir, pero tenemos chupetines y algún dulce, reciclamos
    else if (c > 0 && (x > 0 || y > 0 || z > 0)) {
      # Elegir un dulce disponible al azar para consumir
      disponibles <- c()
      if (x > 0) disponibles <- c(disponibles, "x")
      if (y > 0) disponibles <- c(disponibles, "y")
      if (z > 0) disponibles <- c(disponibles, "z")
      
      dulce_usado <- sample(disponibles, 1)
      if (dulce_usado == "x") x <- x - 1
      if (dulce_usado == "y") y <- y - 1
      if (dulce_usado == "z") z <- z - 1
      
      c <- c - 1  # consumir un chupetín
      
      # Generar 4 dulces aleatorios
      nuevos_dulces <- sample(c("x", "y", "z"), size = 4, replace = TRUE)
      x <- x + sum(nuevos_dulces == "x")
      y <- y + sum(nuevos_dulces == "y")
      z <- z + sum(nuevos_dulces == "z")
      
      interacciones <- interacciones + 1
    }
    else {
      # No se puede hacer nada: atascado
      return(NA)
    }
  }
  
  if (chupetines_producidos >= objetivo_chupetines) {
    return(interacciones)
  } else {
    return(NA)
  }
}

# Ejecutar múltiples simulaciones
set.seed(123)  # Para reproducibilidad
n_sim <- 10000
resultados <- replicate(n_sim, simular_una_vez())

# Filtrar solo las simulaciones exitosas
exitosas <- resultados[!is.na(resultados)]
fracaso <- sum(is.na(resultados))

# Mostrar resultados
cat("Simulaciones exitosas:", length(exitosas), "de", n_sim, "\n")
cat("Tasa de éxito:", round(length(exitosas) / n_sim * 100, 2), "%\n\n")

if (length(exitosas) > 0) {
  cat("Estadísticas de interacciones (solo casos exitosos):\n")
  cat("Mínimo:", min(exitosas), "\n")
  cat("Mediana:", median(exitosas), "\n")
  cat("Promedio:", round(mean(exitosas), 2), "\n")
  cat("Máximo:", max(exitosas), "\n")
  
  # Opcional: histograma
  hist(exitosas, breaks = 20, col = "lightblue", 
       main = "Distribución de interacciones necesarias",
       xlab = "Número de interacciones", 
       ylab = "Frecuencia")
}
