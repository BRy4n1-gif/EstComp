# Instalar paquetes necesarios si no están instalados
if (!requireNamespace("profvis", quietly = TRUE)) install.packages("profvis")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(ggplot2)
library(gridExtra)
library(dplyr)

# 1. Linear Congruential Generator (LCG)
lcg <- function(n, seed = 12345, a = 1664525, c = 1013904223, m = 2^32) {
  result <- numeric(n)
  result[1] <- seed
  for (i in 2:n) {
    result[i] <- (a * result[i-1] + c) %% m
  }
  return(result / m)  # Normalizar a [0,1]
}

# 2. Mersenne Twister - IMPLEMENTACIÓN CORREGIDA (sin funciones anidadas)
mersenne_twister_simple <- function(n, seed = 5489) {
  w <- 32
  n_mt <- 624
  m <- 397
  a <- 0x9908B0DF
  f <- 1812433253
  
  # Inicializar estado
  mt <- integer(n_mt)
  mt[1] <- as.integer(seed)
  for (i in 2:n_mt) {
    mt[i] <- as.integer((f * bitwXor(mt[i-1], bitwShiftR(mt[i-1], 30)) + (i - 1)) %% (2^32))
  }
  
  # Generar números
  result <- numeric(n)
  index <- 1
  
  for (k in 1:n) {
    # Twist si es necesario
    if (index == 1) {
      for (i in 1:n_mt) {
        y <- bitwOr(
          bitwAnd(mt[i], 0x80000000),
          bitwAnd(mt[ifelse(i < n_mt, i + 1, 1)], 0x7FFFFFFF)
        )
        mt[i] <- bitwXor(
          mt[ifelse(i + m <= n_mt, i + m, i + m - n_mt)],
          bitwShiftR(y, 1)
        )
        if (bitwAnd(y, 1) != 0) {
          mt[i] <- bitwXor(mt[i], a)
        }
      }
    }
    
    y <- mt[index]
    y <- bitwXor(y, bitwShiftR(y, 11))
    y <- bitwXor(y, bitwAnd(bitwShiftL(y, 7), 0x9D2C5680))
    y <- bitwXor(y, bitwAnd(bitwShiftL(y, 15), 0xEFC60000))
    y <- bitwXor(y, bitwShiftR(y, 18))
    
    # Normalizar a [0,1]
    result[k] <- (as.numeric(bitwAnd(y, 0xFFFFFFFF)) + 2^31) / (2^32)
    
    index <- ifelse(index < n_mt, index + 1, 1)
  }
  
  return(result)
}

# 3. Xorshift
xorshift <- function(n, seed = 123456789) {
  result <- numeric(n)
  x <- as.integer(seed)
  
  for (i in 1:n) {
    x <- bitwXor(x, bitwShiftL(x, 13))
    x <- bitwXor(x, bitwShiftR(x, 17))
    x <- bitwXor(x, bitwShiftL(x, 5))
    # Asegurar positivo y normalizar
    result[i] <- (as.numeric(bitwAnd(x, 0x7FFFFFFF))) / (2^31)
  }
  return(result)
}

# 4. Middle Square (versión simple)
middle_square <- function(n, seed = 1234) {
  result <- numeric(n)
  current <- as.numeric(seed)
  
  for (i in 1:n) {
    current <- current^2
    current_str <- sprintf("%08.0f", current)
    # Extraer 4 dígitos centrales
    current <- as.numeric(substr(current_str, 3, 6))
    if (is.na(current)) current <- 0
    result[i] <- current / 10000
  }
  return(result)
}

# Función para evaluar uniformidad
evaluar_uniformidad <- function(numeros, nombre) {
  numeros <- na.omit(numeros)
  numeros <- numeros[is.finite(numeros)]
  numeros <- numeros[numeros >= 0 & numeros <= 1]
  
  if (length(numeros) < 100) {
    warning(paste("Muy pocos datos válidos para", nombre))
    return(NULL)
  }
  
  breaks <- seq(0, 1, length.out = 11)
  observed <- hist(numeros, breaks = breaks, plot = FALSE)$counts
  expected <- rep(length(numeros)/10, 10)
  chi_sq <- sum((observed - expected)^2 / expected)
  p_value <- 1 - pchisq(chi_sq, df = 9)
  
  media <- mean(numeros)
  varianza <- var(numeros)
  
  n_plot <- min(1000, length(numeros))
  df_hist <- data.frame(x = numeros[1:n_plot])
  df_scatter <- data.frame(
    x = numeros[1:(n_plot-1)], 
    y = numeros[2:n_plot]
  )
  
  p1 <- ggplot(df_hist, aes(x = x)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", alpha = 0.7, color = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
    labs(title = paste("Distribución -", nombre),
         subtitle = sprintf("Media: %.3f, Varianza: %.3f", media, varianza),
         x = "Valor", y = "Densidad") +
    theme_minimal()
  
  p2 <- ggplot(df_scatter, aes(x = x, y = y)) +
    geom_point(alpha = 0.3, size = 0.8, color = "blue") +
    labs(title = paste("Dispersión consecutiva -", nombre),
         x = "X[i]", y = "X[i+1]") +
    theme_minimal()
  
  return(list(
    chi_sq = chi_sq, 
    p_value = p_value, 
    media = media,
    varianza = varianza,
    plot1 = p1, 
    plot2 = p2,
    n_validos = length(numeros)
  ))
}

# Función para evaluar rendimiento
evaluar_rendimiento <- function(n_large = 10000) {
  tiempo_lcg <- system.time(lcg(n_large))[3]
  tiempo_mt  <- system.time(mersenne_twister_simple(n_large))[3]
  tiempo_xs  <- system.time(xorshift(n_large))[3]
  tiempo_ms  <- system.time(middle_square(n_large))[3]
  tiempo_r   <- system.time(runif(n_large))[3]
  
  data.frame(
    Generador = c("LCG", "Mersenne Twister", "Xorshift", "Middle Square", "R Native"),
    Tiempo = c(tiempo_lcg, tiempo_mt, tiempo_xs, tiempo_ms, tiempo_r)
  )
}

# ANÁLISIS PRINCIPAL
set.seed(123)
n <- 10000

cat("Generando números aleatorios...\n")

# Generar con manejo de errores
lcg_nums <- tryCatch(lcg(n), error = function(e) runif(n))
mt_nums  <- tryCatch(mersenne_twister_simple(n), error = function(e) runif(n))
xs_nums  <- tryCatch(xorshift(n), error = function(e) runif(n))
ms_nums  <- tryCatch(middle_square(n), error = function(e) runif(n))
r_nums   <- runif(n)

# Evaluar uniformidad
cat("Evaluando uniformidad...\n")
eval_lcg <- evaluar_uniformidad(lcg_nums, "LCG")
eval_mt  <- evaluar_uniformidad(mt_nums, "Mersenne Twister")
eval_xs  <- evaluar_uniformidad(xs_nums, "Xorshift")
eval_ms  <- evaluar_uniformidad(ms_nums, "Middle Square")
eval_r   <- evaluar_uniformidad(r_nums, "R Native")

# Mostrar resultados
resultados_uniformidad <- data.frame(
  Generador = c("LCG", "Mersenne Twister", "Xorshift", "Middle Square", "R Native"),
  Chi_Cuadrado = c(eval_lcg$chi_sq, eval_mt$chi_sq, eval_xs$chi_sq, eval_ms$chi_sq, eval_r$chi_sq),
  P_Value = c(eval_lcg$p_value, eval_mt$p_value, eval_xs$p_value, eval_ms$p_value, eval_r$p_value),
  Media = c(eval_lcg$media, eval_mt$media, eval_xs$media, eval_ms$media, eval_r$media),
  Varianza = c(eval_lcg$varianza, eval_mt$varianza, eval_xs$varianza, eval_ms$varianza, eval_r$varianza),
  N_Validos = c(eval_lcg$n_validos, eval_mt$n_validos, eval_xs$n_validos, eval_ms$n_validos, eval_r$n_validos)
)

print("RESULTADOS DE UNIFORMIDAD:")
print(resultados_uniformidad)

# Gráficos
cat("Generando gráficos...\n")
if (all(sapply(list(eval_lcg, eval_mt, eval_xs, eval_ms, eval_r), Negate(is.null)))) {
  grid.arrange(
    eval_lcg$plot1, eval_mt$plot1, eval_xs$plot1, 
    eval_ms$plot1, eval_r$plot1,
    ncol = 3, nrow = 2,
    top = "Distribución de Números Aleatorios"
  )
  
  grid.arrange(
    eval_lcg$plot2, eval_mt$plot2, eval_xs$plot2, 
    eval_ms$plot2, eval_r$plot2,
    ncol = 3, nrow = 2,
    top = "Dispersión de Puntos Consecutivos"
  )
}

# Rendimiento
cat("Evaluando rendimiento...\n")
resultados_rendimiento <- evaluar_rendimiento(50000)
print("RESULTADOS DE RENDIMIENTO (50,000 números):")
print(resultados_rendimiento)

# Gráfico de rendimiento
p_rendimiento <- ggplot(resultados_rendimiento, 
                        aes(x = reorder(Generador, Tiempo), y = Tiempo, fill = Generador)) +
  geom_bar(stat = "identity") +
  labs(title = "Tiempo de Ejecución por Generador",
       x = "Generador", y = "Tiempo (segundos)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_rendimiento)

# Profiling (opcional)
cat("Realizando profiling...\n")
if (requireNamespace("profvis", quietly = TRUE)) {
  tryCatch({
    profvis::profvis({
      n_test <- 10000
      lcg(n_test)
      mersenne_twister_simple(n_test)
      xorshift(n_test)
      middle_square(n_test)
      runif(n_test)
    })
  }, error = function(e) cat("Profiling falló:", e$message, "\n"))
} else {
  cat("Paquete 'profvis' no disponible. Saltando profiling.\n")
}

cat("Análisis completado.\n")
