set.seed(123) 

n <- 100
datos <- data.frame(
  Edad = round(rnorm(n, mean = 35, sd = 8)),
  Ingresos = round(rnorm(n, mean = 2500, sd = 800)),
  HorasdeTrabajo = round(rnorm(n, mean = 40, sd = 5)),
  Experiencia = round(rnorm(n, mean = 10, sd = 4)),
  Satisfaccion = round(runif(n, min = 1, max = 10)),
  Productividad = round(rnorm(n, mean = 85, sd = 12))
)

print(summary(datos))

cat("\nPRIMERAS 10 FILAS DEL DATASET:\n")
print(head(datos, 10))

calcular_cv <- function(x) {
  x <- na.omit(x)
  cv <- (sd(x) / mean(x)) * 100
  return(round(cv, 2))
}

cv_resultados <- sapply(datos, calcular_cv)

resultados_cv <- data.frame(
  Variable = names(cv_resultados),
  Coeficiente_Variacion = cv_resultados
)

print(resultados_cv)


for (var in names(datos)) {
  cat("\n==", var, "==\n")
  cat("Media:", round(mean(datos[[var]], na.rm = TRUE), 2), "\n")
  cat("Desviación estándar:", round(sd(datos[[var]], na.rm = TRUE), 2), "\n")
  cat("Coeficiente de variación:", calcular_cv(datos[[var]]), "%\n")
  cat("Mínimo:", min(datos[[var]], na.rm = TRUE), "\n")
  cat("Máximo:", max(datos[[var]], na.rm = TRUE), "\n")
}
