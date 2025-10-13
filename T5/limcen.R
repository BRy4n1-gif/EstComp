# =============================================================================
# TESTS PARA TENDENCIA CENTRAL - UNA MUESTRA
# =============================================================================

# One-sample t-test
satisfaccion <- c(6.5, 7.8, 8.1, 6.9, 7.3, 8.0, 6.7, 7.5, 8.2, 7.1)
t.test(satisfaccion, mu = 7.2)

# One-sample z-test (requiere instalar BSDA)
# install.packages("BSDA")  # Descomenta si no lo tienes instalado
library(BSDA)
pesos <- c(502, 498, 505, 495, 501, 499, 503, 497, 500, 504)
z.test(pesos, mu = 500, sigma.x = 3) # CORREGIDO: sigma.x más realista

# Wilcoxon signed-rank test
tiempo_reaccion <- c(245, 240, 255, 235, 248, 252, 238, 242, 250, 246)
wilcox.test(tiempo_reaccion, mu = 250)

# Sign test (usando binom.test)
mejoras <- c(1, 1, 0, 1, 1, 1, 0, 1, 1, 0)
binom.test(sum(mejoras), length(mejoras), p = 0.5)

# =============================================================================
# TESTS COMPARANDO DOS GRUPOS INDEPENDIENTES
# =============================================================================

# Independent samples t-test
control <- c(18, 20, 16, 22, 19, 17, 21, 18, 20, 19)
experimental <- c(14, 12, 15, 13, 11, 16, 12, 14, 13, 15)
t.test(control, experimental, var.equal = TRUE)

# Welch's t-test
universidad_A <- c(3.2, 4.1, 3.8, 4.5, 3.9, 4.2, 3.7, 4.0)
universidad_B <- c(2.8, 3.1, 2.9, 3.3, 2.7, 3.0, 2.9, 3.2, 2.8)
t.test(universidad_A, universidad_B, var.equal = FALSE)

# Mann-Whitney U test
sucursal_norte <- c(4, 5, 3, 4, 5, 4, 3, 5, 4)
sucursal_sur <- c(3, 2, 4, 3, 3, 2, 4, 3, 2, 3)
wilcox.test(sucursal_norte, sucursal_sur)

# =============================================================================
# TESTS PARA GRUPOS RELACIONADOS/PAREADOS
# =============================================================================

# Paired t-test
antes <- c(140, 135, 142, 138, 145, 139, 141, 136, 143, 137)
despues <- c(132, 128, 135, 130, 138, 131, 134, 129, 136, 130)
t.test(antes, despues, paired = TRUE)

# Wilcoxon signed-rank (pareado)
dolor_antes <- c(8, 7, 9, 6, 8, 7, 9, 8, 7, 6)
dolor_despues <- c(5, 4, 6, 3, 5, 4, 7, 5, 4, 3)
wilcox.test(dolor_antes, dolor_despues, paired = TRUE)

# =============================================================================
# TESTS COMPARANDO MÚLTIPLES GRUPOS
# =============================================================================

# One-way ANOVA
metodo_A <- c(85, 88, 82, 87, 84, 86, 83, 89, 85, 87)
metodo_B <- c(78, 80, 76, 82, 79, 77, 81, 78, 80, 79)
metodo_C <- c(92, 94, 90, 93, 91, 95, 89, 92, 94, 91)

datos <- data.frame(
  rendimiento = c(metodo_A, metodo_B, metodo_C),
  metodo = factor(rep(c("A", "B", "C"), each = 10))
)
resultado_anova <- aov(rendimiento ~ metodo, data = datos)
summary(resultado_anova)
TukeyHSD(resultado_anova)

# Kruskal-Wallis test
ventas <- c(5, 6, 4, 5, 6, 5, 4, 6, 5)
marketing <- c(3, 4, 3, 4, 3, 4, 3, 4)
finanzas <- c(6, 7, 6, 7, 6, 7, 6, 7, 6)
rrhh <- c(4, 5, 4, 5, 4, 5, 4, 5, 4, 5)

datos_estres <- data.frame(
  estres = c(ventas, marketing, finanzas, rrhh),
  departamento = factor(rep(c("Ventas", "Marketing", "Finanzas", "RRHH"), 
                            c(length(ventas), length(marketing), 
                              length(finanzas), length(rrhh)))) # CORREGIDO: usa length()
)
kruskal.test(estres ~ departamento, data = datos_estres)

# Repeated measures ANOVA
participante <- factor(rep(1:8, 3))
musica <- factor(rep(c("Clasica", "Ambiental", "Silencio"), each = 8))
concentracion <- c(
  c(25, 28, 22, 30, 26, 29, 24, 27),
  c(30, 32, 28, 35, 31, 33, 29, 31),
  c(20, 22, 18, 24, 21, 23, 19, 22)
)

datos_rm <- data.frame(participante, musica, concentracion) # CORREGIDO: orden correcto
resultado_rm <- aov(concentracion ~ musica + Error(participante/musica), data = datos_rm)
summary(resultado_rm)

# Friedman test
preferencias <- matrix(c(
  4, 5, 2, 3, 1,
  5, 4, 1, 3, 2,
  3, 5, 2, 4, 1,
  4, 5, 1, 2, 3,
  5, 4, 2, 3, 1,
  4, 5, 1, 3, 2
), nrow = 6, byrow = TRUE)

friedman.test(preferencias)

# =============================================================================
# TESTS PARA VARIANZA
# =============================================================================

# F-test para igualdad de varianzas
maquina_1 <- c(12.5, 12.8, 12.3, 12.7, 12.4, 12.6, 12.5, 12.9)
maquina_2 <- c(12.1, 13.2, 11.8, 13.5, 12.0, 13.1, 11.9, 13.3)
var.test(maquina_1, maquina_2)

# Levene's test (requiere car)
# install.packages("car")  # Descomenta si no lo tienes instalado
library(car)
region_norte <- c(150, 145, 155, 148, 152, 149, 151)
region_centro <- c(120, 135, 115, 140, 125, 130, 118)
region_sur <- c(180, 165, 190, 175, 185, 170, 195)

datos_ventas <- data.frame(
  ventas = c(region_norte, region_centro, region_sur),
  region = factor(rep(c("Norte", "Centro", "Sur"), each = 7))
)
leveneTest(ventas ~ region, data = datos_ventas)

# Bartlett's test
grupo_20_30 <- c(75, 78, 72, 76, 74, 77, 73, 79)
grupo_31_40 <- c(68, 70, 65, 72, 69, 71, 67, 73)
grupo_41_50 <- c(62, 65, 60, 67, 63, 66, 61, 68)
grupo_51_60 <- c(58, 60, 55, 62, 59, 61, 57, 63)

datos_creatividad <- data.frame(
  creatividad = c(grupo_20_30, grupo_31_40, grupo_41_50, grupo_51_60),
  edad = factor(rep(c("20-30", "31-40", "41-50", "51-60"), each = 8))
)
bartlett.test(creatividad ~ edad, data = datos_creatividad)
