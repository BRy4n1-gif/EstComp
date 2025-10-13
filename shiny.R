library(shiny)
library(ggplot2)
library(dplyr)
library(nortest)

ui <- fluidPage(
  titlePanel("Análisis con Prueba t de Student"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Configuración de Datos"),
      
      selectInput("data_source", "Fuente de datos:",
                  choices = c("Cargar archivo CSV" = "csv",
                              "Datos de ejemplo" = "example",
                              "Ingresar manualmente" = "manual")),
      
      conditionalPanel(
        condition = "input.data_source == 'csv'",
        fileInput("file", "Seleccionar archivo CSV:",
                  accept = c(".csv"))
      ),
      
      conditionalPanel(
        condition = "input.data_source == 'manual'",
        textAreaInput("manual_data", "Ingrese datos (separados por comas):",
                      value = "23,25,28,30,22,26,29,24,27,25"),
        textAreaInput("manual_data2", "Grupo 2 (opcional, para t de 2 muestras):",
                      value = "")
      ),
      
      hr(),
      
      selectInput("test_type", "Tipo de prueba:",
                  choices = c("Una muestra" = "one",
                              "Dos muestras independientes" = "two",
                              "Muestras pareadas" = "paired")),
      
      checkboxInput("use_nonparametric", "Usar prueba no paramétrica", value = FALSE),
      
      conditionalPanel(
        condition = "input.test_type == 'one'",
        numericInput("mu", "Valor de referencia (μ₀):", value = 25)
      ),
      
      conditionalPanel(
        condition = "input.test_type == 'two' && input.data_source == 'csv'",
        uiOutput("var_selector")
      ),
      
      selectInput("alternative", "Hipótesis alternativa:",
                  choices = c("Bilateral (≠)" = "two.sided",
                              "Mayor que (>)" = "greater",
                              "Menor que (<)" = "less")),
      
      sliderInput("conf_level", "Nivel de confianza:",
                  min = 0.80, max = 0.99, value = 0.95, step = 0.01),
      
      actionButton("analyze", "Analizar", class = "btn-primary btn-block")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Resultados",
                 h3("Diagnóstico de Supuestos"),
                 verbatimTextOutput("assumptions"),
                 plotOutput("qq_plot", height = "300px"),
                 
                 hr(),
                 conditionalPanel(
                   condition = "!input.use_nonparametric",
                   h3("Resultados de la Prueba t (Paramétrica)"),
                   verbatimTextOutput("test_results")
                 ),
                 
                 conditionalPanel(
                   condition = "input.use_nonparametric",
                   h3("Resultados de la Prueba No Paramétrica"),
                   verbatimTextOutput("nonparam_results")
                 ),
                 
                 hr(),
                 h3("Comparación de Métodos"),
                 verbatimTextOutput("method_comparison"),
                 
                 hr(),
                 h3("Recomendaciones"),
                 uiOutput("recommendations")
        ),
        
        tabPanel("Visualización",
                 plotOutput("dist_plot", height = "400px"),
                 plotOutput("box_plot", height = "300px")
        ),
        
        tabPanel("Datos",
                 h4("Resumen de los datos"),
                 verbatimTextOutput("data_summary"),
                 hr(),
                 tableOutput("data_table")
        ),
        
        tabPanel("Ayuda",
                 h4("Guía de uso"),
                 HTML("
                 <h5>Tipos de Prueba t:</h5>
                 <ul>
                   <li><b>Una muestra:</b> Compara la media de un grupo con un valor de referencia</li>
                   <li><b>Dos muestras independientes:</b> Compara las medias de dos grupos diferentes</li>
                   <li><b>Muestras pareadas:</b> Compara dos mediciones del mismo grupo (antes/después)</li>
                 </ul>
                 
                 <h5>Pruebas No Paramétricas (alternativas):</h5>
                 <ul>
                   <li><b>Wilcoxon (una muestra):</b> Alternativa cuando los datos no son normales</li>
                   <li><b>Mann-Whitney (dos muestras):</b> Compara medianas de grupos independientes</li>
                   <li><b>Wilcoxon pareado:</b> Para datos pareados sin normalidad</li>
                 </ul>
                 
                 <h5>Supuestos:</h5>
                 <ul>
                   <li>Normalidad: Los datos deben seguir una distribución normal</li>
                   <li>Independencia: Las observaciones deben ser independientes</li>
                   <li>Homocedasticidad: Varianzas iguales (para dos muestras independientes)</li>
                 </ul>
                 
                 <h5>Interpretación del p-valor:</h5>
                 <ul>
                   <li>p < 0.05: Evidencia significativa contra H₀</li>
                   <li>p ≥ 0.05: No hay evidencia suficiente para rechazar H₀</li>
                 </ul>
                 ")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  data <- reactiveVal(NULL)
  
  observeEvent(input$data_source, {
    if(input$data_source == "example") {
      set.seed(123)
      df <- data.frame(
        grupo1 = rnorm(30, mean = 25, sd = 3),
        grupo2 = rnorm(30, mean = 27, sd = 3)
      )
      data(df)
    }
  })
  
  observeEvent(input$file, {
    req(input$file)
    df <- read.csv(input$file$datapath)
    data(df)
  })
  
  output$var_selector <- renderUI({
    req(data())
    df <- data()
    num_vars <- names(df)[sapply(df, is.numeric)]
    
    tagList(
      selectInput("var1", "Variable grupo 1:", choices = num_vars),
      selectInput("var2", "Variable grupo 2:", choices = num_vars, selected = num_vars[2])
    )
  })
  
  get_data_for_analysis <- reactive({
    if(input$data_source == "manual") {
      g1 <- as.numeric(unlist(strsplit(input$manual_data, ",")))
      g1 <- g1[!is.na(g1)]
      
      if(input$test_type != "one" && input$manual_data2 != "") {
        g2 <- as.numeric(unlist(strsplit(input$manual_data2, ",")))
        g2 <- g2[!is.na(g2)]
        list(group1 = g1, group2 = g2)
      } else {
        list(group1 = g1, group2 = NULL)
      }
    } else {
      req(data())
      df <- data()
      if(input$test_type == "one") {
        list(group1 = df[[1]], group2 = NULL)
      } else {
        req(input$var1, input$var2)
        list(group1 = df[[input$var1]], group2 = df[[input$var2]])
      }
    }
  })
  
  observeEvent(input$analyze, {
    
    output$assumptions <- renderPrint({
      dat <- get_data_for_analysis()
      g1 <- dat$group1
      g2 <- dat$group2
      
      cat("=== EVALUACIÓN DE SUPUESTOS ===\n\n")
      
      # Normalidad grupo 1
      cat("1. Normalidad - Grupo 1:\n")
      if(length(g1) >= 3) {
        sw_test <- shapiro.test(g1)
        cat(sprintf("   Shapiro-Wilk: W = %.4f, p-valor = %.4f\n", 
                    sw_test$statistic, sw_test$p.value))
        if(sw_test$p.value > 0.05) {
          cat("   ✓ Los datos parecen seguir una distribución normal\n\n")
        } else {
          cat("   ✗ Los datos NO siguen una distribución normal\n\n")
        }
      }
      
      # Normalidad grupo 2
      if(!is.null(g2) && length(g2) >= 3) {
        cat("2. Normalidad - Grupo 2:\n")
        sw_test2 <- shapiro.test(g2)
        cat(sprintf("   Shapiro-Wilk: W = %.4f, p-valor = %.4f\n", 
                    sw_test2$statistic, sw_test2$p.value))
        if(sw_test2$p.value > 0.05) {
          cat("   ✓ Los datos parecen seguir una distribución normal\n\n")
        } else {
          cat("   ✗ Los datos NO siguen una distribución normal\n\n")
        }
      }
      
      # Homocedasticidad
      if(!is.null(g2) && input$test_type == "two") {
        cat("3. Homocedasticidad (igualdad de varianzas):\n")
        var_test <- var.test(g1, g2)
        cat(sprintf("   F-test: F = %.4f, p-valor = %.4f\n", 
                    var_test$statistic, var_test$p.value))
        if(var_test$p.value > 0.05) {
          cat("   ✓ Las varianzas son homogéneas\n")
          cat("   Recomendación: Usar var.equal = TRUE\n\n")
        } else {
          cat("   ✗ Las varianzas NO son homogéneas\n")
          cat("   Recomendación: Usar var.equal = FALSE (Welch)\n\n")
        }
      }
    })
    
    output$qq_plot <- renderPlot({
      dat <- get_data_for_analysis()
      g1 <- dat$group1
      g2 <- dat$group2
      
      if(!is.null(g2)) {
        par(mfrow = c(1, 2))
        qqnorm(g1, main = "Q-Q Plot Grupo 1")
        qqline(g1, col = "red")
        qqnorm(g2, main = "Q-Q Plot Grupo 2")
        qqline(g2, col = "red")
      } else {
        qqnorm(g1, main = "Q-Q Plot")
        qqline(g1, col = "red")
      }
    })
    
    output$test_results <- renderPrint({
      dat <- get_data_for_analysis()
      g1 <- dat$group1
      g2 <- dat$group2
      
      cat("=== RESULTADOS DE LA PRUEBA t (PARAMÉTRICA) ===\n\n")
      
      if(input$test_type == "one") {
        result <- t.test(g1, mu = input$mu, 
                         alternative = input$alternative,
                         conf.level = input$conf_level)
        cat(sprintf("Prueba t de una muestra\n"))
        cat(sprintf("H₀: μ = %.2f\n", input$mu))
        
      } else if(input$test_type == "two") {
        # Decidir si usar var.equal basado en test de varianzas
        var_equal <- var.test(g1, g2)$p.value > 0.05
        
        result <- t.test(g1, g2, 
                         alternative = input$alternative,
                         conf.level = input$conf_level,
                         var.equal = var_equal)
        cat(sprintf("Prueba t de dos muestras independientes\n"))
        cat(sprintf("Varianzas iguales: %s\n", ifelse(var_equal, "Sí", "No (Welch)")))
        
      } else {
        result <- t.test(g1, g2, 
                         paired = TRUE,
                         alternative = input$alternative,
                         conf.level = input$conf_level)
        cat(sprintf("Prueba t de muestras pareadas\n"))
      }
      
      cat(sprintf("\nEstadístico t = %.4f\n", result$statistic))
      cat(sprintf("Grados de libertad = %.2f\n", result$parameter))
      cat(sprintf("p-valor = %.5f\n", result$p.value))
      cat(sprintf("\nIntervalo de confianza al %d%%: [%.4f, %.4f]\n", 
                  input$conf_level * 100, result$conf.int[1], result$conf.int[2]))
      cat(sprintf("Media estimada: %.4f\n", result$estimate[1]))
      
      if(!is.null(g2) && input$test_type != "paired") {
        cat(sprintf("Diferencia de medias: %.4f\n", result$estimate[1] - result$estimate[2]))
      }
      
      cat("\n")
      if(result$p.value < 0.05) {
        cat("CONCLUSIÓN: Rechazamos H₀ (p < 0.05)\n")
        cat("Hay evidencia estadísticamente significativa.\n")
      } else {
        cat("CONCLUSIÓN: No rechazamos H₀ (p ≥ 0.05)\n")
        cat("No hay evidencia estadísticamente significativa.\n")
      }
    })
    
    output$nonparam_results <- renderPrint({
      dat <- get_data_for_analysis()
      g1 <- dat$group1
      g2 <- dat$group2
      
      cat("=== RESULTADOS DE LA PRUEBA NO PARAMÉTRICA ===\n\n")
      
      if(input$test_type == "one") {
        result <- wilcox.test(g1, mu = input$mu, 
                              alternative = input$alternative,
                              conf.level = input$conf_level,
                              conf.int = TRUE)
        cat(sprintf("Prueba de Wilcoxon de una muestra\n"))
        cat(sprintf("H₀: mediana = %.2f\n", input$mu))
        cat(sprintf("\nEstadístico V = %.0f\n", result$statistic))
        
      } else if(input$test_type == "two") {
        result <- wilcox.test(g1, g2, 
                              alternative = input$alternative,
                              conf.level = input$conf_level,
                              conf.int = TRUE)
        cat(sprintf("Prueba de Mann-Whitney (Wilcoxon de rango)\n"))
        cat(sprintf("H₀: las medianas son iguales\n"))
        cat(sprintf("\nEstadístico W = %.0f\n", result$statistic))
        
      } else {
        result <- wilcox.test(g1, g2, 
                              paired = TRUE,
                              alternative = input$alternative,
                              conf.level = input$conf_level,
                              conf.int = TRUE)
        cat(sprintf("Prueba de Wilcoxon pareada\n"))
        cat(sprintf("H₀: la mediana de las diferencias es cero\n"))
        cat(sprintf("\nEstadístico V = %.0f\n", result$statistic))
      }
      
      cat(sprintf("p-valor = %.5f\n", result$p.value))
      
      if(!is.null(result$conf.int)) {
        cat(sprintf("\nIntervalo de confianza al %d%%: [%.4f, %.4f]\n", 
                    input$conf_level * 100, result$conf.int[1], result$conf.int[2]))
      }
      
      if(!is.null(result$estimate)) {
        cat(sprintf("Estimación: %.4f\n", result$estimate))
      }
      
      cat("\n")
      if(result$p.value < 0.05) {
        cat("CONCLUSIÓN: Rechazamos H₀ (p < 0.05)\n")
        cat("Hay evidencia estadísticamente significativa.\n")
      } else {
        cat("CONCLUSIÓN: No rechazamos H₀ (p ≥ 0.05)\n")
        cat("No hay evidencia estadísticamente significativa.\n")
      }
      
      cat("\n")
      cat("NOTA: Las pruebas no paramétricas trabajan con rangos y\n")
      cat("medianas en lugar de medias, siendo más robustas ante\n")
      cat("valores atípicos y distribuciones no normales.\n")
    })
    
    output$method_comparison <- renderPrint({
      dat <- get_data_for_analysis()
      g1 <- dat$group1
      g2 <- dat$group2
      
      cat("=== COMPARACIÓN DE MÉTODOS ===\n\n")
      
      # Prueba paramétrica
      if(input$test_type == "one") {
        t_result <- t.test(g1, mu = input$mu, alternative = input$alternative)
        w_result <- wilcox.test(g1, mu = input$mu, alternative = input$alternative)
        cat(sprintf("Prueba t (paramétrica):      p-valor = %.5f\n", t_result$p.value))
        cat(sprintf("Wilcoxon (no paramétrica):   p-valor = %.5f\n\n", w_result$p.value))
        
      } else if(input$test_type == "two") {
        t_result <- t.test(g1, g2, alternative = input$alternative)
        w_result <- wilcox.test(g1, g2, alternative = input$alternative)
        cat(sprintf("Prueba t (paramétrica):      p-valor = %.5f\n", t_result$p.value))
        cat(sprintf("Mann-Whitney (no paramétrica): p-valor = %.5f\n\n", w_result$p.value))
        
      } else {
        t_result <- t.test(g1, g2, paired = TRUE, alternative = input$alternative)
        w_result <- wilcox.test(g1, g2, paired = TRUE, alternative = input$alternative)
        cat(sprintf("Prueba t pareada (paramétrica):    p-valor = %.5f\n", t_result$p.value))
        cat(sprintf("Wilcoxon pareada (no paramétrica): p-valor = %.5f\n\n", w_result$p.value))
      }
      
      # Evaluación de normalidad
      sw_test <- shapiro.test(g1)
      
      cat("RECOMENDACIÓN ÓPTIMA:\n")
      cat("─────────────────────\n")
      
      if(sw_test$p.value > 0.05 && length(g1) >= 30) {
        cat("✓ Los datos son normales y el tamaño de muestra es adecuado\n")
        cat("→ USE LA PRUEBA PARAMÉTRICA (t de Student)\n")
        cat("  Razón: Mayor poder estadístico y precisión en estimaciones\n")
      } else if(sw_test$p.value < 0.05) {
        cat("✗ Los datos NO siguen una distribución normal\n")
        cat("→ USE LA PRUEBA NO PARAMÉTRICA (Wilcoxon/Mann-Whitney)\n")
        cat("  Razón: Es más robusta y no asume normalidad\n")
      } else if(length(g1) < 30) {
        cat("⚠ Tamaño de muestra pequeño (n < 30)\n")
        if(sw_test$p.value > 0.05) {
          cat("→ Ambas pruebas son válidas, pero considere:\n")
          cat("  • Prueba t: si está seguro de la normalidad\n")
          cat("  • Prueba no paramétrica: opción más conservadora\n")
        } else {
          cat("→ USE LA PRUEBA NO PARAMÉTRICA\n")
          cat("  Razón: Muestra pequeña Y falta de normalidad\n")
        }
      }
      
      # Comparación de conclusiones
      cat("\n")
      if((t_result$p.value < 0.05) == (w_result$p.value < 0.05)) {
        cat("✓ Ambos métodos llegan a la MISMA conclusión\n")
        cat("  Esto aumenta la confianza en los resultados\n")
      } else {
        cat("⚠ Los métodos llegan a conclusiones DIFERENTES\n")
        cat("  Considere factores como:\n")
        cat("  • Presencia de valores atípicos\n")
        cat("  • Forma de la distribución\n")
        cat("  • Tamaño del efecto vs significancia estadística\n")
      }
    })
    
    output$recommendations <- renderUI({
      dat <- get_data_for_analysis()
      g1 <- dat$group1
      g2 <- dat$group2
      
      recommendations <- list()
      
      # Verificar normalidad
      if(length(g1) >= 3) {
        sw_test <- shapiro.test(g1)
        if(sw_test$p.value < 0.05) {
          recommendations <- c(recommendations, 
                               "<li><b>⚠ Violación de normalidad:</b> Considere usar pruebas no paramétricas (Wilcoxon/Mann-Whitney)</li>")
        }
      }
      
      # Verificar tamaño de muestra
      if(length(g1) < 30) {
        recommendations <- c(recommendations,
                             "<li><b>Tamaño de muestra pequeño:</b> Los resultados son más sensibles a la violación de supuestos. Considere recolectar más datos.</li>")
      }
      
      # Verificar homocedasticidad
      if(!is.null(g2) && input$test_type == "two") {
        var_test <- var.test(g1, g2)
        if(var_test$p.value < 0.05) {
          recommendations <- c(recommendations,
                               "<li><b>Varianzas desiguales:</b> Se recomienda usar la corrección de Welch (var.equal = FALSE), que es más robusta.</li>")
        }
      }
      
      # Interpretación del tamaño del efecto
      result <- if(input$test_type == "one") {
        t.test(g1, mu = input$mu)
      } else if(input$test_type == "paired") {
        t.test(g1, g2, paired = TRUE)
      } else {
        t.test(g1, g2)
      }
      
      # Cohen's d aproximado
      if(!is.null(g2)) {
        pooled_sd <- sqrt((var(g1) + var(g2)) / 2)
        cohen_d <- abs(mean(g1) - mean(g2)) / pooled_sd
        effect_size <- if(cohen_d < 0.2) "trivial" else if(cohen_d < 0.5) "pequeño" else if(cohen_d < 0.8) "mediano" else "grande"
        recommendations <- c(recommendations,
                             sprintf("<li><b>Tamaño del efecto (Cohen's d = %.2f):</b> %s. Considere la relevancia práctica además de la significancia estadística.</li>", 
                                     cohen_d, effect_size))
      }
      
      if(length(recommendations) == 0) {
        recommendations <- "<li><b>✓ Todos los supuestos se cumplen razonablemente.</b> Los resultados de la prueba t son válidos.</li>"
      }
      
      HTML(paste0(
        "<div class='alert alert-info'>",
        "<h4>Análisis y Recomendaciones:</h4>",
        "<ul>",
        paste(recommendations, collapse = ""),
        "</ul>",
        "</div>"
      ))
    })
    
    output$dist_plot <- renderPlot({
      dat <- get_data_for_analysis()
      g1 <- dat$group1
      g2 <- dat$group2
      
      if(!is.null(g2) && input$test_type != "one") {
        df_plot <- data.frame(
          value = c(g1, g2),
          group = factor(rep(c("Grupo 1", "Grupo 2"), c(length(g1), length(g2))))
        )
        
        ggplot(df_plot, aes(x = value, fill = group)) +
          geom_density(alpha = 0.5) +
          theme_minimal() +
          labs(title = "Distribución de los Datos",
               x = "Valor", y = "Densidad", fill = "Grupo") +
          scale_fill_manual(values = c("#3498db", "#e74c3c"))
      } else {
        df_plot <- data.frame(value = g1)
        ggplot(df_plot, aes(x = value)) +
          geom_density(fill = "#3498db", alpha = 0.5) +
          geom_vline(xintercept = mean(g1), color = "#e74c3c", linetype = "dashed", size = 1) +
          geom_vline(xintercept = input$mu, color = "black", linetype = "dotted", size = 1) +
          theme_minimal() +
          labs(title = "Distribución de los Datos",
               subtitle = "Línea roja: media observada | Línea negra: valor de referencia",
               x = "Valor", y = "Densidad")
      }
    })
    
    output$box_plot <- renderPlot({
      dat <- get_data_for_analysis()
      g1 <- dat$group1
      g2 <- dat$group2
      
      if(!is.null(g2) && input$test_type != "one") {
        df_plot <- data.frame(
          value = c(g1, g2),
          group = factor(rep(c("Grupo 1", "Grupo 2"), c(length(g1), length(g2))))
        )
        
        ggplot(df_plot, aes(x = group, y = value, fill = group)) +
          geom_boxplot(alpha = 0.7) +
          geom_jitter(width = 0.2, alpha = 0.3) +
          theme_minimal() +
          labs(title = "Diagrama de Cajas",
               x = "Grupo", y = "Valor") +
          scale_fill_manual(values = c("#3498db", "#e74c3c")) +
          theme(legend.position = "none")
      } else {
        df_plot <- data.frame(value = g1)
        ggplot(df_plot, aes(x = "", y = value)) +
          geom_boxplot(fill = "#3498db", alpha = 0.7) +
          geom_jitter(width = 0.2, alpha = 0.3) +
          theme_minimal() +
          labs(title = "Diagrama de Cajas",
               x = "", y = "Valor")
      }
    })
    
    output$data_summary <- renderPrint({
      dat <- get_data_for_analysis()
      g1 <- dat$group1
      g2 <- dat$group2
      
      cat("Grupo 1:\n")
      cat(sprintf("  n = %d\n", length(g1)))
      cat(sprintf("  Media = %.4f\n", mean(g1)))
      cat(sprintf("  Desviación estándar = %.4f\n", sd(g1)))
      cat(sprintf("  Mediana = %.4f\n", median(g1)))
      cat(sprintf("  Rango = [%.4f, %.4f]\n\n", min(g1), max(g1)))
      
      if(!is.null(g2)) {
        cat("Grupo 2:\n")
        cat(sprintf("  n = %d\n", length(g2)))
        cat(sprintf("  Media = %.4f\n", mean(g2)))
        cat(sprintf("  Desviación estándar = %.4f\n", sd(g2)))
        cat(sprintf("  Mediana = %.4f\n", median(g2)))
        cat(sprintf("  Rango = [%.4f, %.4f]\n", min(g2), max(g2)))
      }
    })
    
    output$data_table <- renderTable({
      dat <- get_data_for_analysis()
      g1 <- dat$group1
      g2 <- dat$group2
      
      if(!is.null(g2)) {
        max_len <- max(length(g1), length(g2))
        df <- data.frame(
          Grupo_1 = c(g1, rep(NA, max_len - length(g1))),
          Grupo_2 = c(g2, rep(NA, max_len - length(g2)))
        )
        head(df, 20)
      } else {
        data.frame(Valor = head(g1, 20))
      }
    })
  })
}

shinyApp(ui = ui, server = server)
