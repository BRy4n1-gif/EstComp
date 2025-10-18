# ============================================================================
# Aplicación Shiny: Análisis de Incertidumbre en Métricas de Cambio Climático
# Métodos: Monte Carlo, Inferencia Bayesiana, MCMC
# ============================================================================

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(ggplot2)
library(dplyr)
library(tidyr)
library(coda)
library(plotly)
library(DT)

# ============================================================================
# FUNCIONES AUXILIARES
# ============================================================================

# Generar datos de ejemplo (temperaturas globales)
generar_datos_ejemplo <- function() {
  set.seed(123)
  years <- 1900:2020
  temp <- 13.5 + 0.008 * (years - 1900) + 0.00002 * (years - 1900)^2 + 
    rnorm(length(years), 0, 0.15)
  data.frame(año = years, valor = temp)
}

# Simulación Monte Carlo MEJORADA
monte_carlo_sim <- function(data, n_sim = 1000, error_sd = 0.1, año_pred = 2030) {
  years <- data$año
  values <- data$valor
  n <- length(values)
  
  model <- lm(values ~ years)
  coef_est <- coef(model)
  
  # Simulaciones con bootstrap y propagación de errores
  sims <- replicate(n_sim, {
    noisy_values <- values + rnorm(n, 0, error_sd)
    sim_model <- lm(noisy_values ~ years)
    predict(sim_model, newdata = data.frame(years = año_pred))
  })
  
  # Calcular también la tendencia (pendiente)
  slopes <- replicate(n_sim, {
    noisy_values <- values + rnorm(n, 0, error_sd)
    sim_model <- lm(noisy_values ~ years)
    coef(sim_model)[2]
  })
  
  return(list(
    simulations = sims,
    slopes = slopes,
    mean_pred = mean(sims),
    median_pred = median(sims),
    sd_pred = sd(sims),
    ci_lower = quantile(sims, 0.05),
    ci_upper = quantile(sims, 0.95),
    slope_mean = mean(slopes),
    slope_sd = sd(slopes),
    slope_ci_lower = quantile(slopes, 0.05),
    slope_ci_upper = quantile(slopes, 0.95)
  ))
}

# Inferencia Bayesiana analítica
bayesian_inference <- function(data, prior_mean = 0.01, prior_sd = 0.005) {
  years <- data$año
  values <- data$valor
  years_scaled <- scale(years, center = TRUE, scale = FALSE)
  
  model <- lm(values ~ years_scaled)
  slope_mle <- coef(model)[2]
  slope_se <- summary(model)$coefficients[2, 2]
  
  prior_prec <- 1 / prior_sd^2
  like_prec <- 1 / slope_se^2
  
  post_prec <- prior_prec + like_prec
  post_sd <- sqrt(1 / post_prec)
  post_mean <- (prior_prec * prior_mean + like_prec * slope_mle) / post_prec
  
  cred_lower <- qnorm(0.05, post_mean, post_sd)
  cred_upper <- qnorm(0.95, post_mean, post_sd)
  
  # Calcular factor de Bayes aproximado
  bayes_factor <- prior_sd / post_sd
  
  return(list(
    prior_mean = prior_mean,
    prior_sd = prior_sd,
    likelihood_mean = slope_mle,
    likelihood_sd = slope_se,
    posterior_mean = post_mean,
    posterior_sd = post_sd,
    cred_lower = cred_lower,
    cred_upper = cred_upper,
    bayes_factor = bayes_factor,
    prior_precision = prior_prec,
    like_precision = like_prec,
    post_precision = post_prec
  ))
}

# Algoritmo Metropolis-Hastings MEJORADO
mcmc_regression <- function(data, n_iter = 5000, n_chains = 3, burn_in = 1000) {
  years <- data$año
  values <- data$valor
  n <- length(values)
  years_scaled <- scale(years, center = TRUE, scale = FALSE)[,1]
  
  log_posterior <- function(params) {
    intercept <- params[1]
    slope <- params[2]
    sigma <- params[3]
    
    if (sigma <= 0) return(-Inf)
    
    # Priors más informativos
    log_prior <- dnorm(intercept, mean(values), 5, log = TRUE) +
      dnorm(slope, 0, 0.1, log = TRUE) +
      dgamma(1/sigma^2, 0.01, 0.01, log = TRUE)
    
    mu <- intercept + slope * years_scaled
    log_lik <- sum(dnorm(values, mu, sigma, log = TRUE))
    
    return(log_prior + log_lik)
  }
  
  all_chains <- lapply(1:n_chains, function(chain_id) {
    # Inicialización con perturbación para cada cadena
    init_vals <- c(mean(values) + rnorm(1, 0, 0.1), 
                   0.01 + rnorm(1, 0, 0.002), 
                   sd(values) + rnorm(1, 0, 0.05))
    current <- init_vals
    samples <- matrix(NA, n_iter, 3)
    colnames(samples) <- c("intercept", "slope", "sigma")
    accepted <- 0
    
    # Adaptación de propuesta
    proposal_sd <- c(0.05, 0.002, 0.02)
    
    for (i in 1:n_iter) {
      proposal <- current + rnorm(3, 0, proposal_sd)
      log_alpha <- log_posterior(proposal) - log_posterior(current)
      
      if (log(runif(1)) < log_alpha) {
        current <- proposal
        accepted <- accepted + 1
        
        # Adaptación simple de la propuesta
        if (i %% 100 == 0 && i < burn_in) {
          accept_rate <- accepted / i
          if (accept_rate < 0.2) proposal_sd <- proposal_sd * 0.9
          if (accept_rate > 0.5) proposal_sd <- proposal_sd * 1.1
        }
      }
      
      samples[i, ] <- current
    }
    
    list(samples = samples, 
         acceptance_rate = accepted / n_iter,
         proposal_sd = proposal_sd)
  })
  
  chains_mcmc <- lapply(all_chains, function(x) {
    mcmc(x$samples[(burn_in + 1):n_iter, ])
  })
  
  combined <- do.call(rbind, lapply(chains_mcmc, as.matrix))
  
  r_hat <- tryCatch({
    gelman.diag(mcmc.list(chains_mcmc))$psrf[, 1]
  }, error = function(e) rep(NA, 3))
  
  ess <- tryCatch({
    effectiveSize(mcmc.list(chains_mcmc))
  }, error = function(e) rep(NA, 3))
  
  # Calcular autocorrelación
  autocorr <- tryCatch({
    lapply(chains_mcmc, function(chain) {
      acf(chain[, "slope"], plot = FALSE, lag.max = 50)$acf
    })
  }, error = function(e) NULL)
  
  return(list(
    chains = chains_mcmc,
    combined = combined,
    acceptance_rates = sapply(all_chains, function(x) x$acceptance_rate),
    r_hat = r_hat,
    ess = ess,
    summary = summary(mcmc.list(chains_mcmc)),
    autocorr = autocorr,
    n_iter = n_iter,
    burn_in = burn_in
  ))
}

# ============================================================================
# INTERFAZ DE USUARIO
# ============================================================================

ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(
    title = "Análisis de Incertidumbre Climática",
    titleWidth = 350
  ),
  
  dashboardSidebar(
    width = 280,
    sidebarMenu(
      id = "sidebar",
      menuItem("📊 Datos", tabName = "datos", icon = icon("database")),
      menuItem("🎲 Monte Carlo", tabName = "montecarlo", icon = icon("dice")),
      menuItem("📈 Bayesiano", tabName = "bayesiano", icon = icon("chart-line")),
      menuItem("⛓️ MCMC", tabName = "mcmc", icon = icon("link")),
      menuItem("🔍 Comparación", tabName = "comparacion", icon = icon("balance-scale")),
      menuItem("ℹ️ Ayuda", tabName = "ayuda", icon = icon("question-circle"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .skin-blue .main-header .logo { background-color: #2c3e50; }
        .skin-blue .main-header .navbar { background-color: #34495e; }
        .content-wrapper { background-color: #ecf0f1; }
        .box { border-radius: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        .info-box { border-radius: 5px; }
        .small-box { border-radius: 5px; }
        .btn-primary { background-color: #3498db; border-color: #2980b9; }
        .btn-success { background-color: #27ae60; border-color: #229954; }
        .btn-warning { background-color: #f39c12; border-color: #e67e22; }
        .interpretation-box { 
          background-color: #fff; 
          padding: 15px; 
          border-radius: 5px; 
          margin-bottom: 15px;
          border-left: 4px solid #3498db;
        }
      "))
    ),
    
    tabItems(
      # TAB: DATOS
      tabItem(
        tabName = "datos",
        fluidRow(
          box(
            title = "Cargar Datos", status = "primary", solidHeader = TRUE,
            width = 4, collapsible = TRUE,
            fileInput("file_input", "Subir archivo CSV",
                      accept = c(".csv"),
                      buttonLabel = "Buscar...",
                      placeholder = "Ningún archivo seleccionado"),
            helpText("El archivo debe contener dos columnas: 'año' y 'valor'"),
            hr(),
            actionButton("load_example", "📁 Cargar Datos de Ejemplo", 
                         class = "btn-success btn-block", 
                         icon = icon("file-upload")),
            hr(),
            h5(tags$b("Vista previa:")),
            DTOutput("data_preview")
          ),
          
          box(
            title = "Visualización de Datos", status = "info", solidHeader = TRUE,
            width = 8,
            plotlyOutput("data_plot", height = "350px")
          )
        ),
        
        fluidRow(
          box(
            title = "Estadísticas Descriptivas", status = "warning", solidHeader = TRUE,
            width = 12,
            verbatimTextOutput("data_summary")
          )
        )
      ),
      
      # TAB: MONTE CARLO
      tabItem(
        tabName = "montecarlo",
        fluidRow(
          box(
            title = "Configuración Monte Carlo", status = "primary", solidHeader = TRUE,
            width = 3, collapsible = TRUE,
            sliderInput("mc_nsim", "Número de simulaciones:", 
                        min = 100, max = 10000, value = 1000, step = 100),
            numericInput("mc_error_sd", "Desv. estándar del error:", 
                         value = 0.1, min = 0.01, max = 1, step = 0.01),
            numericInput("mc_year_pred", "Año para predicción:",
                         value = 2030, min = 2021, max = 2100, step = 1),
            hr(),
            actionButton("run_mc", "▶️ Ejecutar Simulación", 
                         class = "btn-success btn-block",
                         icon = icon("play")),
            br(),
            downloadButton("download_mc", "💾 Descargar Resultados", 
                           class = "btn-warning btn-block")
          ),
          
          box(
            title = "Distribución de Resultados", status = "success", solidHeader = TRUE,
            width = 9,
            fluidRow(
              column(6, plotlyOutput("mc_histogram", height = "300px")),
              column(6, plotlyOutput("mc_density", height = "300px"))
            )
          )
        ),
        
        fluidRow(
          valueBoxOutput("mc_mean_box", width = 3),
          valueBoxOutput("mc_median_box", width = 3),
          valueBoxOutput("mc_sd_box", width = 3),
          valueBoxOutput("mc_ci_box", width = 3)
        ),
        
        fluidRow(
          box(
            title = "Tabla de Estadísticas", status = "info", solidHeader = TRUE,
            width = 6,
            DTOutput("mc_summary_table")
          ),
          box(
            title = "Convergencia de Simulaciones", status = "success", solidHeader = TRUE,
            width = 6,
            plotOutput("mc_convergence", height = "250px")
          )
        ),
        
        fluidRow(
          box(
            title = "📖 Interpretación Monte Carlo", status = "warning", solidHeader = TRUE,
            width = 12, collapsible = TRUE,
            uiOutput("mc_interpretation")
          )
        )
      ),
      
      # TAB: BAYESIANO
      tabItem(
        tabName = "bayesiano",
        fluidRow(
          box(
            title = "Configuración Bayesiana", status = "primary", solidHeader = TRUE,
            width = 3, collapsible = TRUE,
            h5(tags$b("Prior (Distribución a priori):")),
            numericInput("bayes_prior_mean", "Media:", 
                         value = 0.01, step = 0.001),
            numericInput("bayes_prior_sd", "Desv. estándar:", 
                         value = 0.005, min = 0.001, max = 0.1, step = 0.001),
            hr(),
            actionButton("run_bayes", "▶️ Ejecutar Inferencia", 
                         class = "btn-success btn-block",
                         icon = icon("play")),
            br(),
            downloadButton("download_bayes", "💾 Descargar Resultados", 
                           class = "btn-warning btn-block")
          ),
          
          box(
            title = "Distribuciones: Prior, Verosimilitud y Posterior", 
            status = "success", solidHeader = TRUE,
            width = 9,
            plotlyOutput("bayes_distributions", height = "400px")
          )
        ),
        
        fluidRow(
          valueBoxOutput("bayes_prior_box", width = 4),
          valueBoxOutput("bayes_posterior_box", width = 4),
          valueBoxOutput("bayes_credible_box", width = 4)
        ),
        
        fluidRow(
          box(
            title = "Resumen de Inferencia", status = "info", solidHeader = TRUE,
            width = 6,
            DTOutput("bayes_summary_table")
          ),
          box(
            title = "Visualización de Actualización", status = "success", solidHeader = TRUE,
            width = 6,
            plotOutput("bayes_update_visual", height = "300px")
          )
        ),
        
        fluidRow(
          box(
            title = "📖 Interpretación Bayesiana", status = "warning", solidHeader = TRUE,
            width = 12, collapsible = TRUE,
            uiOutput("bayes_interpretation")
          )
        )
      ),
      
      # TAB: MCMC
      tabItem(
        tabName = "mcmc",
        fluidRow(
          box(
            title = "Configuración MCMC", status = "primary", solidHeader = TRUE,
            width = 3, collapsible = TRUE,
            sliderInput("mcmc_iter", "Iteraciones:", 
                        min = 1000, max = 20000, value = 5000, step = 1000),
            sliderInput("mcmc_chains", "Cadenas:", 
                        min = 1, max = 5, value = 3, step = 1),
            sliderInput("mcmc_burnin", "Burn-in:", 
                        min = 100, max = 5000, value = 1000, step = 100),
            hr(),
            actionButton("run_mcmc", "▶️ Ejecutar MCMC", 
                         class = "btn-success btn-block",
                         icon = icon("play")),
            br(),
            downloadButton("download_mcmc", "💾 Descargar Muestras", 
                           class = "btn-warning btn-block"),
            hr(),
            uiOutput("mcmc_diagnostics_ui")
          ),
          
          box(
            title = "Selección de Parámetro", status = "info", solidHeader = TRUE,
            width = 9,
            selectInput("mcmc_param", "Parámetro a visualizar:",
                        choices = c("Intercepto" = "intercept",
                                    "Pendiente (Tendencia)" = "slope",
                                    "Sigma (Error)" = "sigma"),
                        selected = "slope"),
            fluidRow(
              column(6, plotOutput("mcmc_trace", height = "280px")),
              column(6, plotOutput("mcmc_density", height = "280px"))
            )
          )
        ),
        
        fluidRow(
          valueBoxOutput("mcmc_rhat_box", width = 4),
          valueBoxOutput("mcmc_ess_box", width = 4),
          valueBoxOutput("mcmc_accept_box", width = 4)
        ),
        
        fluidRow(
          box(
            title = "Diagnósticos Adicionales", status = "info", solidHeader = TRUE,
            width = 6,
            plotOutput("mcmc_autocorr", height = "300px")
          ),
          box(
            title = "Resumen Estadístico Posterior", status = "success", solidHeader = TRUE,
            width = 6,
            verbatimTextOutput("mcmc_summary")
          )
        ),
        
        fluidRow(
          box(
            title = "📖 Interpretación MCMC", status = "warning", solidHeader = TRUE,
            width = 12, collapsible = TRUE,
            uiOutput("mcmc_interpretation")
          )
        )
      ),
      
      # TAB: COMPARACIÓN
      tabItem(
        tabName = "comparacion",
        fluidRow(
          box(
            title = "Comparación Visual de Métodos", status = "primary", solidHeader = TRUE,
            width = 12,
            plotlyOutput("comparison_plot", height = "450px")
          )
        ),
        
        fluidRow(
          infoBoxOutput("comp_mc_box", width = 4),
          infoBoxOutput("comp_bayes_box", width = 4),
          infoBoxOutput("comp_mcmc_box", width = 4)
        ),
        
        fluidRow(
          box(
            title = "Tabla Comparativa Detallada", status = "info", solidHeader = TRUE,
            width = 12,
            DTOutput("comparison_table"),
            hr(),
            p(class = "text-muted", 
              icon("info-circle"), 
              " Esta comparación muestra las estimaciones de tendencia (pendiente) 
              obtenidas mediante los tres métodos. Los intervalos representan 
              incertidumbre al 90%.")
          )
        ),
        
        fluidRow(
          box(
            title = "📊 Análisis Comparativo", status = "warning", solidHeader = TRUE,
            width = 12, collapsible = TRUE,
            uiOutput("comparison_interpretation")
          )
        )
      ),
      
      # TAB: AYUDA
      tabItem(
        tabName = "ayuda",
        fluidRow(
          box(
            title = "Guía de Uso", status = "primary", solidHeader = TRUE,
            width = 12,
            h3("📚 Bienvenido al Análisis de Incertidumbre Climática"),
            hr(),
            
            h4("1️⃣ Cargar Datos"),
            p("Comience cargando sus datos en la pestaña 'Datos'. Puede usar los datos 
              de ejemplo o subir su propio archivo CSV con columnas 'año' y 'valor'."),
            
            h4("2️⃣ Monte Carlo"),
            p("Método de simulación que genera múltiples escenarios aleatorios para 
              cuantificar la incertidumbre en las predicciones futuras."),
            tags$ul(
              tags$li("Ajuste el número de simulaciones (más = mejor precisión)"),
              tags$li("Configure la desviación estándar del error de medición"),
              tags$li("Seleccione el año para predicción"),
              tags$li("Obtenga intervalos de confianza empíricos")
            ),
            
            h4("3️⃣ Inferencia Bayesiana"),
            p("Método analítico que combina conocimiento previo (prior) con datos 
              observados para obtener distribuciones posteriores."),
            tags$ul(
              tags$li("Defina su distribución a priori (media y desv. estándar)"),
              tags$li("Visualice cómo la evidencia actualiza sus creencias"),
              tags$li("Obtenga intervalos de credibilidad"),
              tags$li("Analice la influencia relativa del prior vs. los datos")
            ),
            
            h4("4️⃣ MCMC (Markov Chain Monte Carlo)"),
            p("Método de muestreo para modelos bayesianos complejos usando el 
              algoritmo Metropolis-Hastings con adaptación."),
            tags$ul(
              tags$li("Configure iteraciones, cadenas y período de burn-in"),
              tags$li("Revise diagnósticos de convergencia (R-hat < 1.1 es bueno)"),
              tags$li("Analice trace plots y densidades posteriores"),
              tags$li("Verifique autocorrelación de las cadenas")
            ),
            
            h4("5️⃣ Comparación"),
            p("Compare los resultados de los tres métodos lado a lado para 
              evaluar consistencia y robustez de las estimaciones."),
            
            hr(),
            h4("📊 Interpretación de Resultados"),
            tags$ul(
              tags$li(tags$b("Intervalos de Confianza (MC):"), 
                      " Rango donde esperamos encontrar el 90% de las predicciones"),
              tags$li(tags$b("Intervalos de Credibilidad (Bayesiano):"), 
                      " Probabilidad del 90% de que el parámetro esté en ese rango"),
              tags$li(tags$b("R-hat (MCMC):"), 
                      " Valores < 1.1 indican buena convergencia"),
              tags$li(tags$b("ESS (MCMC):"), 
                      " Tamaño efectivo de muestra (mayor es mejor)"),
              tags$li(tags$b("Autocorrelación:"), 
                      " Debe decaer rápidamente (indica muestras independientes)")
            ),
            
            hr(),
            p(class = "text-center", 
              tags$b("Desarrollado para análisis científico de cambio climático"),
              br(),
              "Versión 3.0 | 2024")
          )
        )
      )
    )
  )
)

# ============================================================================
# SERVIDOR
# ============================================================================

server <- function(input, output, session) {
  
  # ====== DATOS REACTIVOS ======
  data_reactive <- reactiveVal(NULL)
  mc_results <- reactiveVal(NULL)
  bayes_results <- reactiveVal(NULL)
  mcmc_results <- reactiveVal(NULL)
  
  # ====== CARGAR DATOS ======
  observeEvent(input$load_example, {
    data_reactive(generar_datos_ejemplo())
    showNotification("✅ Datos de ejemplo cargados correctamente", 
                     type = "message", duration = 3)
  })
  
  observeEvent(input$file_input, {
    req(input$file_input)
    tryCatch({
      data <- read.csv(input$file_input$datapath)
      
      if (!all(c("año", "valor") %in% names(data))) {
        showNotification("❌ El archivo debe contener columnas 'año' y 'valor'", 
                         type = "error", duration = 5)
        return(NULL)
      }
      
      data_reactive(data)
      showNotification("✅ Datos cargados correctamente", 
                       type = "message", duration = 3)
    }, error = function(e) {
      showNotification(paste("❌ Error al cargar datos:", e$message), 
                       type = "error", duration = 5)
    })
  })
  
  # ====== VISUALIZACIÓN DE DATOS ======
  output$data_preview <- renderDT({
    req(data_reactive())
    datatable(head(data_reactive(), 10), 
              options = list(pageLength = 5, dom = 't'),
              rownames = FALSE)
  })
  
  output$data_plot <- renderPlotly({
    req(data_reactive())
    data <- data_reactive()
    
    p <- ggplot(data, aes(x = año, y = valor)) +
      geom_line(color = "#3498db", size = 1.2) +
      geom_point(color = "#2c3e50", size = 2.5, alpha = 0.6) +
      geom_smooth(method = "lm", color = "#e74c3c", linetype = "dashed", 
                  se = TRUE, fill = "#e74c3c", alpha = 0.15) +
      labs(title = "Serie Temporal de Datos Climáticos",
           x = "Año", y = "Valor", 
           subtitle = "Línea roja: tendencia lineal con IC 95%") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))
    
    ggplotly(p) %>% layout(hovermode = 'x unified')
  })
  
  output$data_summary <- renderPrint({
    req(data_reactive())
    data <- data_reactive()
    model <- lm(valor ~ año, data = data)
    
    cat("═══════════════════════════════════════════════════════\n")
    cat("               RESUMEN ESTADÍSTICO                     \n")
    cat("═══════════════════════════════════════════════════════\n\n")
    cat(sprintf("📊 N° observaciones:      %d\n", nrow(data)))
    cat(sprintf("📅 Rango temporal:        %d - %d\n", min(data$año), max(data$año)))
    cat(sprintf("📈 Media de valores:      %.4f\n", mean(data$valor)))
    cat(sprintf("📉 Desviación estándar:   %.4f\n", sd(data$valor)))
    cat(sprintf("🔼 Valor máximo:          %.4f\n", max(data$valor)))
    cat(sprintf("🔽 Valor mínimo:          %.4f\n", min(data$valor)))
    cat("\n───────────────────────────────────────────────────────\n")
    cat("                 REGRESIÓN LINEAL                      \n")
    cat("───────────────────────────────────────────────────────\n\n")
    cat(sprintf("Pendiente (tendencia):    %.6f ± %.6f\n", 
                coef(model)[2], summary(model)$coefficients[2, 2]))
    cat(sprintf("Intercepto:               %.4f ± %.4f\n",
                coef(model)[1], summary(model)$coefficients[1, 2]))
    cat(sprintf("R² (ajuste):              %.4f\n", summary(model)$r.squared))
    cat(sprintf("P-valor (pendiente):      %.2e\n", 
                summary(model)$coefficients[2, 4]))
    cat("\n═══════════════════════════════════════════════════════\n")
  })
  
  # ====== MONTE CARLO ======
  observeEvent(input$run_mc, {
    req(data_reactive())
    
    withProgress(message = '🎲 Ejecutando simulaciones Monte Carlo...', {
      tryCatch({
        results <- monte_carlo_sim(
          data_reactive(), 
          n_sim = input$mc_nsim,
          error_sd = input$mc_error_sd,
          año_pred = input$mc_year_pred
        )
        mc_results(results)
        updateTabsetPanel(session, "sidebar", selected = "montecarlo")
        showNotification("✅ Simulación Monte Carlo completada", 
                         type = "message", duration = 3)
      }, error = function(e) {
        showNotification(paste("❌ Error en Monte Carlo:", e$message), 
                         type = "error", duration = 5)
      })
    })
  })
  
  output$mc_histogram <- renderPlotly({
    req(mc_results())
    results <- mc_results()
    df <- data.frame(pred = results$simulations)
    
    p <- ggplot(df, aes(x = pred)) +
      geom_histogram(bins = 50, fill = "#3498db", color = "#2c3e50", alpha = 0.7) +
      geom_vline(xintercept = results$mean_pred, color = "#e74c3c", 
                 linetype = "dashed", size = 1.2) +
      geom_vline(xintercept = results$ci_lower, color = "#f39c12", 
                 linetype = "dotted", size = 1) +
      geom_vline(xintercept = results$ci_upper, color = "#f39c12", 
                 linetype = "dotted", size = 1) +
      labs(title = paste("Histograma de Predicciones", input$mc_year_pred),
           x = "Valor Predicho", y = "Frecuencia") +
      theme_minimal(base_size = 12)
    
    ggplotly(p) %>% layout(hovermode = 'x')
  })
  
  output$mc_density <- renderPlotly({
    req(mc_results())
    results <- mc_results()
    df <- data.frame(pred = results$simulations)
    
    p <- ggplot(df, aes(x = pred)) +
      geom_density(fill = "#27ae60", alpha = 0.5, size = 1.2) +
      geom_vline(xintercept = results$median_pred, color = "#2980b9", 
                 linetype = "dashed", size = 1.2) +
      labs(title = "Densidad de Probabilidad",
           x = "Valor Predicho", y = "Densidad") +
      theme_minimal(base_size = 12)
    
    ggplotly(p) %>% layout(hovermode = 'x')
  })
  
  output$mc_convergence <- renderPlot({
    req(mc_results())
    results <- mc_results()
    
    # Calcular media acumulada
    cumulative_means <- cumsum(results$slopes) / seq_along(results$slopes)
    
    par(mar = c(4, 4, 3, 2))
    plot(cumulative_means, type = "l", col = "#3498db", lwd = 2,
         xlab = "Número de Simulaciones",
         ylab = "Media Acumulada de Pendiente",
         main = "Convergencia de Monte Carlo",
         cex.main = 1.2, cex.lab = 1.1)
    abline(h = results$slope_mean, col = "#e74c3c", lty = 2, lwd = 2)
    legend("topright", 
           legend = c("Media acumulada", "Media final"),
           col = c("#3498db", "#e74c3c"),
           lty = c(1, 2), lwd = 2, bty = "n")
    grid()
  })
  
  output$mc_mean_box <- renderValueBox({
    req(mc_results())
    valueBox(
      format(mc_results()$mean_pred, digits = 4),
      "Media de Predicciones",
      icon = icon("chart-line"),
      color = "blue"
    )
  })
  
  output$mc_median_box <- renderValueBox({
    req(mc_results())
    valueBox(
      format(mc_results()$median_pred, digits = 4),
      "Mediana",
      icon = icon("chart-bar"),
      color = "green"
    )
  })
  
  output$mc_sd_box <- renderValueBox({
    req(mc_results())
    valueBox(
      format(mc_results()$sd_pred, digits = 4),
      "Desviación Estándar",
      icon = icon("arrows-alt-h"),
      color = "orange"
    )
  })
  
  output$mc_ci_box <- renderValueBox({
    req(mc_results())
    results <- mc_results()
    valueBox(
      format(results$ci_upper - results$ci_lower, digits = 4),
      "Amplitud IC 90%",
      icon = icon("expand"),
      color = "red"
    )
  })
  
  output$mc_summary_table <- renderDT({
    req(mc_results())
    results <- mc_results()
    
    df <- data.frame(
      Estadística = c("Media", "Mediana", "Desv. Estándar", 
                      "IC 5%", "IC 95%", "Amplitud IC"),
      Valor = c(
        results$mean_pred,
        results$median_pred,
        results$sd_pred,
        results$ci_lower,
        results$ci_upper,
        results$ci_upper - results$ci_lower
      )
    )
    
    datatable(df, options = list(dom = 't'), rownames = FALSE) %>%
      formatRound('Valor', digits = 6)
  })
  
  output$mc_interpretation <- renderUI({
    req(mc_results())
    results <- mc_results()
    
    cv <- (results$sd_pred / results$mean_pred) * 100
    range_pct <- ((results$ci_upper - results$ci_lower) / results$mean_pred) * 100
    
    tagList(
      h4(icon("chart-bar"), " Interpretación Monte Carlo"),
      hr(),
      
      div(class = "interpretation-box",
          h5(tags$b("🎯 Predicción Central:")),
          p(sprintf(
            "La predicción media para el año %d es %.4f con una mediana de %.4f. 
            La cercanía entre media y mediana (diferencia: %.4f) indica una distribución %s.",
            input$mc_year_pred, results$mean_pred, results$median_pred,
            abs(results$mean_pred - results$median_pred),
            ifelse(abs(results$mean_pred - results$median_pred) < 0.01, 
                   "simétrica y confiable", "con cierta asimetría")
          ))
      ),
      
      div(class = "interpretation-box",
          h5(tags$b("📊 Incertidumbre:")),
          p(sprintf(
            "La desviación estándar de %.4f representa un coeficiente de variación del %.2f%%. 
            Esto indica una incertidumbre %s en las predicciones.",
            results$sd_pred, cv,
            ifelse(cv < 5, "baja", ifelse(cv < 10, "moderada", "alta"))
          ))
      ),
      
      div(class = "interpretation-box",
          h5(tags$b("🎲 Intervalo de Confianza 90%%:")),
          p(sprintf(
            "Con 90%% de confianza, el valor predicho estará entre %.4f y %.4f. 
            La amplitud del intervalo (%.4f, equivalente al %.2f%% del valor medio) refleja 
            la %s precisión de la predicción.",
            results$ci_lower, results$ci_upper, 
            results$ci_upper - results$ci_lower, range_pct,
            ifelse(range_pct < 5, "alta", ifelse(range_pct < 15, "moderada", "baja"))
          ))
      ),
      
      div(class = "interpretation-box",
          h5(tags$b("📈 Tendencia Estimada:")),
          p(sprintf(
            "La pendiente promedio estimada es %.6f ± %.6f (IC 90%%: [%.6f, %.6f]). 
            Esto representa un cambio anual de %.6f unidades.",
            results$slope_mean, results$slope_sd,
            results$slope_ci_lower, results$slope_ci_upper,
            results$slope_mean
          ))
      )
    )
  })
  
  output$download_mc <- downloadHandler(
    filename = function() {
      paste0("monte_carlo_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(mc_results())
      results <- mc_results()
      df <- data.frame(
        simulacion = 1:length(results$simulations),
        prediccion = results$simulations,
        pendiente = results$slopes
      )
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  # ====== BAYESIANO ======
  observeEvent(input$run_bayes, {
    req(data_reactive())
    
    tryCatch({
      results <- bayesian_inference(
        data_reactive(),
        prior_mean = input$bayes_prior_mean,
        prior_sd = input$bayes_prior_sd
      )
      bayes_results(results)
      updateTabsetPanel(session, "sidebar", selected = "bayesiano")
      showNotification("✅ Inferencia Bayesiana completada", 
                       type = "message", duration = 3)
    }, error = function(e) {
      showNotification(paste("❌ Error en inferencia:", e$message), 
                       type = "error", duration = 5)
    })
  })
  
  output$bayes_distributions <- renderPlotly({
    req(bayes_results())
    results <- bayes_results()
    
    all_means <- c(results$prior_mean, results$likelihood_mean, results$posterior_mean)
    all_sds <- c(results$prior_sd, results$likelihood_sd, results$posterior_sd)
    
    x_min <- min(all_means - 4 * all_sds)
    x_max <- max(all_means + 4 * all_sds)
    
    x <- seq(x_min, x_max, length.out = 500)
    
    df <- data.frame(
      x = rep(x, 3),
      density = c(
        dnorm(x, results$prior_mean, results$prior_sd),
        dnorm(x, results$likelihood_mean, results$likelihood_sd),
        dnorm(x, results$posterior_mean, results$posterior_sd)
      ),
      type = factor(rep(c("Prior", "Verosimilitud", "Posterior"), each = length(x)),
                    levels = c("Prior", "Verosimilitud", "Posterior"))
    )
    
    posterior_data <- df[df$type == "Posterior", ]
    ci_data <- posterior_data[posterior_data$x >= results$cred_lower & 
                                posterior_data$x <= results$cred_upper, ]
    
    p <- ggplot() +
      geom_ribbon(data = ci_data, 
                  aes(x = x, ymin = 0, ymax = density),
                  fill = "#e74c3c", alpha = 0.2) +
      geom_line(data = df, 
                aes(x = x, y = density, color = type),
                size = 1.5) +
      geom_vline(xintercept = results$prior_mean, 
                 color = "#3498db", linetype = "dotted", size = 0.8, alpha = 0.6) +
      geom_vline(xintercept = results$likelihood_mean, 
                 color = "#27ae60", linetype = "dashed", size = 0.8, alpha = 0.6) +
      geom_vline(xintercept = results$posterior_mean, 
                 color = "#e74c3c", linetype = "solid", size = 1, alpha = 0.8) +
      geom_vline(xintercept = c(results$cred_lower, results$cred_upper),
                 color = "#c0392b", linetype = "dotted", size = 0.5, alpha = 0.5) +
      scale_color_manual(
        values = c("Prior" = "#3498db", 
                   "Verosimilitud" = "#27ae60", 
                   "Posterior" = "#e74c3c"),
        name = "Distribución"
      ) +
      labs(
        title = "Actualización Bayesiana: Prior → Posterior",
        subtitle = sprintf("IC 90%% Posterior: [%.6f, %.6f]", 
                           results$cred_lower, results$cred_upper),
        x = "Tendencia Climática (pendiente)", 
        y = "Densidad de Probabilidad"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        panel.grid.minor = element_blank()
      )
    
    ggplotly(p, tooltip = c("x", "y", "colour")) %>% 
      layout(hovermode = 'x unified')
  })
  
  output$bayes_update_visual <- renderPlot({
    req(bayes_results())
    results <- bayes_results()
    
    par(mfrow = c(1, 1), mar = c(4, 4, 3, 2))
    
    x_min <- min(results$prior_mean - 4*results$prior_sd, 
                 results$posterior_mean - 4*results$posterior_sd)
    x_max <- max(results$prior_mean + 4*results$prior_sd,
                 results$posterior_mean + 4*results$posterior_sd)
    x <- seq(x_min, x_max, length.out = 500)
    
    plot(x, dnorm(x, results$prior_mean, results$prior_sd),
         type = "l", col = "#3498db", lwd = 3, lty = 2,
         xlab = "Tendencia (pendiente)", 
         ylab = "Densidad",
         main = "Actualización Bayesiana",
         ylim = c(0, max(c(dnorm(x, results$prior_mean, results$prior_sd),
                           dnorm(x, results$posterior_mean, results$posterior_sd)))),
         cex.main = 1.3, cex.lab = 1.1)
    
    lines(x, dnorm(x, results$posterior_mean, results$posterior_sd),
          col = "#e74c3c", lwd = 3)
    
    arrows(x0 = results$prior_mean, 
           y0 = dnorm(results$prior_mean, results$prior_mean, results$prior_sd) * 0.7,
           x1 = results$posterior_mean,
           y1 = dnorm(results$posterior_mean, results$posterior_mean, results$posterior_sd) * 0.7,
           col = "#f39c12", lwd = 2, length = 0.15, code = 2)
    
    legend("topright", 
           legend = c("Prior (creencia inicial)", 
                      "Posterior (actualizada con datos)",
                      "Cambio en la creencia"),
           col = c("#3498db", "#e74c3c", "#f39c12"),
           lty = c(2, 1, 1), lwd = c(3, 3, 2),
           bty = "n", cex = 0.9)
    
    abline(v = results$prior_mean, col = "#3498db", lty = 3, lwd = 1.5)
    abline(v = results$posterior_mean, col = "#e74c3c", lty = 3, lwd = 1.5)
    
    text(results$prior_mean, 
         max(dnorm(x, results$prior_mean, results$prior_sd)) * 0.95,
         sprintf("Prior\n%.6f", results$prior_mean),
         col = "#3498db", cex = 0.8, pos = 2)
    
    text(results$posterior_mean, 
         max(dnorm(x, results$posterior_mean, results$posterior_sd)) * 0.95,
         sprintf("Posterior\n%.6f", results$posterior_mean),
         col = "#e74c3c", cex = 0.8, pos = 4)
  })
  
  output$bayes_prior_box <- renderValueBox({
    req(bayes_results())
    results <- bayes_results()
    valueBox(
      sprintf("%.6f", results$prior_mean),
      "Prior (Media)",
      icon = icon("brain"),
      color = "blue"
    )
  })
  
  output$bayes_posterior_box <- renderValueBox({
    req(bayes_results())
    results <- bayes_results()
    valueBox(
      sprintf("%.6f", results$posterior_mean),
      "Posterior (Media)",
      icon = icon("certificate"),
      color = "red"
    )
  })
  
  output$bayes_credible_box <- renderValueBox({
    req(bayes_results())
    results <- bayes_results()
    valueBox(
      sprintf("%.6f", results$cred_upper - results$cred_lower),
      "Amplitud IC 90%",
      icon = icon("ruler-horizontal"),
      color = "orange"
    )
  })
  
  output$bayes_summary_table <- renderDT({
    req(bayes_results())
    results <- bayes_results()
    
    df <- data.frame(
      Distribución = c("Prior", "Verosimilitud", "Posterior", 
                       "IC Posterior 5%", "IC Posterior 95%"),
      Media = c(results$prior_mean, results$likelihood_mean, 
                results$posterior_mean, results$cred_lower, results$cred_upper),
      Desv_Estándar = c(results$prior_sd, results$likelihood_sd, 
                        results$posterior_sd, NA, NA)
    )
    
    datatable(df, options = list(dom = 't'), rownames = FALSE) %>%
      formatRound(c('Media', 'Desv_Estándar'), digits = 8)
  })
  
  output$bayes_interpretation <- renderUI({
    req(bayes_results())
    results <- bayes_results()
    
    change_pct <- abs((results$posterior_mean - results$prior_mean) / results$prior_mean * 100)
    direction <- if(results$posterior_mean > results$prior_mean) "aumentó" else "disminuyó"
    uncertainty_reduction <- (1 - results$posterior_sd / results$prior_sd) * 100
    
    informativeness <- if(results$bayes_factor > 2) {
      "muy informativos"
    } else if(results$bayes_factor > 1.5) {
      "moderadamente informativos"
    } else {
      "poco informativos"
    }
    
    tagList(
      h4(icon("lightbulb"), " Interpretación Bayesiana"),
      hr(),
      
      div(class = "interpretation-box",
          h5(tags$b("🎯 Actualización de Creencias:")),
          p(sprintf(
            "La estimación posterior de la tendencia es %.6f, lo que representa un cambio del %.2f%% 
            respecto a la creencia inicial (prior: %.6f). Los datos observados %s nuestra estimación.",
            results$posterior_mean, change_pct, results$prior_mean, direction
          ))
      ),
      
      div(class = "interpretation-box",
          h5(tags$b("📊 Reducción de Incertidumbre:")),
          p(sprintf(
            "La desviación estándar se redujo de %.6f (prior) a %.6f (posterior), 
            representando una reducción de incertidumbre del %.2f%%. Los datos observados fueron %s 
            (factor de Bayes: %.2f).",
            results$prior_sd, results$posterior_sd, uncertainty_reduction, 
            informativeness, results$bayes_factor
          ))
      ),
      
      div(class = "interpretation-box",
          h5(tags$b("📈 Intervalo de Credibilidad 90%%:")),
          p(sprintf(
            "Existe un 90%% de probabilidad (interpretación bayesiana) de que la verdadera tendencia 
            climática esté entre %.6f y %.6f (por año). Este intervalo tiene una amplitud de %.6f.",
            results$cred_lower, results$cred_upper, results$cred_upper - results$cred_lower
          ))
      ),
      
      div(class = "interpretation-box",
          h5(tags$b("🔍 Influencia Relativa:")),
          p(sprintf(
            "Precisión del prior: %.2f | Precisión de la verosimilitud: %.2f | Precisión posterior: %.2f. 
            %s tuvo mayor influencia en la estimación final (proporción: %.1f%% vs %.1f%%).",
            results$prior_precision, results$like_precision, results$post_precision,
            if(results$like_precision > results$prior_precision) "Los datos observados" else "El prior",
            (results$like_precision / results$post_precision) * 100,
            (results$prior_precision / results$post_precision) * 100
          ))
      )
    )
  })
  
  output$download_bayes <- downloadHandler(
    filename = function() {
      paste0("bayesian_inference_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(bayes_results())
      results <- bayes_results()
      
      x <- seq(results$posterior_mean - 5*results$posterior_sd,
               results$posterior_mean + 5*results$posterior_sd,
               length.out = 1000)
      
      df <- data.frame(
        tendencia = x,
        prior = dnorm(x, results$prior_mean, results$prior_sd),
        likelihood = dnorm(x, results$likelihood_mean, results$likelihood_sd),
        posterior = dnorm(x, results$posterior_mean, results$posterior_sd)
      )
      
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  # ====== MCMC ======
  observeEvent(input$run_mcmc, {
    req(data_reactive())
    
    withProgress(message = '⛓️ Ejecutando MCMC...', value = 0, {
      tryCatch({
        incProgress(0.3, detail = "Muestreando cadenas...")
        
        results <- mcmc_regression(
          data_reactive(),
          n_iter = input$mcmc_iter,
          n_chains = input$mcmc_chains,
          burn_in = input$mcmc_burnin
        )
        
        incProgress(0.7, detail = "Calculando diagnósticos...")
        mcmc_results(results)
        updateTabsetPanel(session, "sidebar", selected = "mcmc")
        showNotification("✅ MCMC completado exitosamente", 
                         type = "message", duration = 3)
      }, error = function(e) {
        showNotification(paste("❌ Error en MCMC:", e$message), 
                         type = "error", duration = 5)
      })
    })
  })
  
  output$mcmc_diagnostics_ui <- renderUI({
    req(mcmc_results())
    results <- mcmc_results()
    
    tagList(
      h5(tags$b("📊 Diagnósticos de Convergencia")),
      hr(),
      tags$div(
        style = "font-size: 12px;",
        p(tags$b("Tasas de Aceptación:")),
        lapply(1:length(results$acceptance_rates), function(i) {
          rate <- results$acceptance_rates[i] * 100
          color <- if(rate > 20 & rate < 50) "green" else "orange"
          p(style = paste0("color:", color),
            sprintf("Cadena %d: %.1f%%", i, rate))
        }),
        hr(),
        p(tags$b("R-hat (convergencia):")),
        p(style = paste0("color:", if(!is.na(results$r_hat[1]) && results$r_hat[1] < 1.1) "green" else "red"),
          sprintf("Intercepto: %.3f", results$r_hat[1])),
        p(style = paste0("color:", if(!is.na(results$r_hat[2]) && results$r_hat[2] < 1.1) "green" else "red"),
          sprintf("Pendiente: %.3f", results$r_hat[2])),
        p(style = paste0("color:", if(!is.na(results$r_hat[3]) && results$r_hat[3] < 1.1) "green" else "red"),
          sprintf("Sigma: %.3f", results$r_hat[3]))
      )
    )
  })
  
  output$mcmc_trace <- renderPlot({
    req(mcmc_results(), input$mcmc_param)
    
    tryCatch({
      results <- mcmc_results()
      
      param_names <- c("intercept", "slope", "sigma")
      param_labels <- c("Intercepto", "Pendiente", "Sigma")
      param_idx <- which(param_names == input$mcmc_param)
      
      chains_df <- do.call(rbind, lapply(1:length(results$chains), function(i) {
        chain_matrix <- as.matrix(results$chains[[i]])
        data.frame(
          iteration = 1:nrow(chain_matrix),
          value = chain_matrix[, param_idx],
          chain = as.factor(i)
        )
      }))
      
      ggplot(chains_df, aes(x = iteration, y = value, color = chain)) +
        geom_line(alpha = 0.7, size = 0.5) +
        scale_color_brewer(palette = "Set1") +
        labs(title = paste("Trace Plot:", param_labels[param_idx]),
             x = "Iteración", y = "Valor",
             subtitle = "Cadenas MCMC después del burn-in") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "top",
              plot.title = element_text(face = "bold"))
    }, error = function(e) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(1, 1, paste("Error:", e$message), col = "red", cex = 1.2)
    })
  })
  
  output$mcmc_density <- renderPlot({
    req(mcmc_results(), input$mcmc_param)
    
    tryCatch({
      results <- mcmc_results()
      
      param_names <- c("intercept", "slope", "sigma")
      param_labels <- c("Intercepto", "Pendiente", "Sigma")
      param_idx <- which(param_names == input$mcmc_param)
      
      df <- data.frame(value = results$combined[, param_idx])
      
      ggplot(df, aes(x = value)) +
        geom_density(fill = "#9b59b6", alpha = 0.6, size = 1.2) +
        geom_vline(xintercept = mean(df$value), color = "#e74c3c", 
                   linetype = "dashed", size = 1.2) +
        geom_vline(xintercept = quantile(df$value, c(0.05, 0.95)), 
                   color = "#f39c12", linetype = "dotted", size = 1) +
        labs(title = paste("Densidad Posterior:", param_labels[param_idx]),
             x = "Valor", y = "Densidad",
             subtitle = "Línea roja: media | Líneas naranjas: IC 90%") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold"))
    }, error = function(e) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(1, 1, paste("Error:", e$message), col = "red", cex = 1.2)
    })
  })
  
  output$mcmc_autocorr <- renderPlot({
    req(mcmc_results())
    results <- mcmc_results()
    
    if (!is.null(results$autocorr)) {
      par(mar = c(4, 4, 3, 2))
      
      # Graficar autocorrelación de la primera cadena
      acf_vals <- results$autocorr[[1]]
      plot(0:(length(acf_vals)-1), acf_vals, type = "h",
           lwd = 2, col = "#3498db",
           xlab = "Lag", ylab = "Autocorrelación",
           main = "Autocorrelación de la Pendiente (Cadena 1)",
           cex.main = 1.2, cex.lab = 1.1)
      abline(h = 0, lty = 2)
      abline(h = c(-1.96/sqrt(length(results$chains[[1]])), 
                   1.96/sqrt(length(results$chains[[1]]))),
             lty = 2, col = "red")
      grid()
    } else {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(1, 1, "No hay datos de autocorrelación disponibles", cex = 1.2)
    }
  })
  
  output$mcmc_rhat_box <- renderValueBox({
    req(mcmc_results())
    results <- mcmc_results()
    rhat_slope <- results$r_hat[2]
    color <- if(!is.na(rhat_slope) && rhat_slope < 1.1) "green" else "red"
    
    valueBox(
      format(rhat_slope, digits = 3),
      "R-hat (Pendiente)",
      icon = icon("check-circle"),
      color = color
    )
  })
  
  output$mcmc_ess_box <- renderValueBox({
    req(mcmc_results())
    results <- mcmc_results()
    ess_slope <- results$ess[2]
    color <- if(!is.na(ess_slope) && ess_slope > 100) "green" else "orange"
    
    valueBox(
      format(round(ess_slope), big.mark = ","),
      "ESS (Pendiente)",
      icon = icon("layer-group"),
      color = color
    )
  })
  
  output$mcmc_accept_box <- renderValueBox({
    req(mcmc_results())
    results <- mcmc_results()
    avg_accept <- mean(results$acceptance_rates) * 100
    color <- if(avg_accept > 20 && avg_accept < 50) "green" else "orange"
    
    valueBox(
      paste0(format(avg_accept, digits = 3), "%"),
      "Tasa Aceptación Media",
      icon = icon("thumbs-up"),
      color = color
    )
  })
  
  output$mcmc_summary <- renderPrint({
    req(mcmc_results())
    results <- mcmc_results()
    
    cat("═══════════════════════════════════════════════════════\n")
    cat("          RESUMEN DE MUESTRAS POSTERIORES              \n")
    cat("═══════════════════════════════════════════════════════\n\n")
    print(results$summary)
    cat("\n───────────────────────────────────────────────────────\n")
    cat("           DIAGNÓSTICOS DE CONVERGENCIA                \n")
    cat("───────────────────────────────────────────────────────\n\n")
    cat(sprintf("R-hat (Intercepto):  %.4f %s\n", 
                results$r_hat[1],
                if(!is.na(results$r_hat[1]) && results$r_hat[1] < 1.1) "✓" else "✗"))
    cat(sprintf("R-hat (Pendiente):   %.4f %s\n", 
                results$r_hat[2],
                if(!is.na(results$r_hat[2]) && results$r_hat[2] < 1.1) "✓" else "✗"))
    cat(sprintf("R-hat (Sigma):       %.4f %s\n", 
                results$r_hat[3],
                if(!is.na(results$r_hat[3]) && results$r_hat[3] < 1.1) "✓" else "✗"))
    cat("\n")
    cat(sprintf("ESS (Intercepto):    %.0f\n", results$ess[1]))
    cat(sprintf("ESS (Pendiente):     %.0f\n", results$ess[2]))
    cat(sprintf("ESS (Sigma):         %.0f\n", results$ess[3]))
    cat("\n═══════════════════════════════════════════════════════\n")
  })
  
  output$mcmc_interpretation <- renderUI({
    req(mcmc_results())
    results <- mcmc_results()
    
    # Evaluar calidad de convergencia
    convergence_quality <- if(all(!is.na(results$r_hat)) && all(results$r_hat < 1.1)) {
      "excelente"
    } else if(all(!is.na(results$r_hat)) && all(results$r_hat < 1.2)) {
      "buena"
    } else {
      "problemática"
    }
    
    # Evaluar ESS
    ess_quality <- if(all(!is.na(results$ess)) && all(results$ess > 1000)) {
      "muy alto"
    } else if(all(!is.na(results$ess)) && all(results$ess > 400)) {
      "adecuado"
    } else {
      "bajo"
    }
    
    # Evaluar tasa de aceptación
    avg_accept <- mean(results$acceptance_rates) * 100
    accept_quality <- if(avg_accept >= 20 && avg_accept <= 50) {
      "óptima"
    } else if(avg_accept >= 15 && avg_accept <= 60) {
      "aceptable"
    } else {
      "subóptima"
    }
    
    # Estadísticas de la pendiente
    slope_mean <- mean(results$combined[, "slope"])
    slope_sd <- sd(results$combined[, "slope"])
    slope_ci <- quantile(results$combined[, "slope"], c(0.05, 0.95))
    
    tagList(
      h4(icon("link"), " Interpretación MCMC"),
      hr(),
      
      div(class = "interpretation-box",
          h5(tags$b("⛓️ Convergencia de Cadenas:")),
          p(sprintf(
            "La convergencia de las cadenas es %s. Valores de R-hat: Intercepto=%.3f, 
            Pendiente=%.3f, Sigma=%.3f. %s",
            convergence_quality,
            results$r_hat[1], results$r_hat[2], results$r_hat[3],
            if(convergence_quality == "excelente") {
              "Todos los R-hat < 1.1 indican convergencia perfecta."
            } else if(convergence_quality == "buena") {
              "Los R-hat están en rango aceptable. Considere aumentar iteraciones para mayor precisión."
            } else {
              "¡ADVERTENCIA! Algunos R-hat > 1.1. Se recomienda aumentar iteraciones o burn-in."
            }
          ))
      ),
      
      div(class = "interpretation-box",
          h5(tags$b("📊 Tamaño Efectivo de Muestra (ESS):")),
          p(sprintf(
            "El ESS es %s: Intercepto=%.0f, Pendiente=%.0f, Sigma=%.0f. 
            Con %d iteraciones totales (después de burn-in de %d), esto representa una 
            eficiencia del %.1f%% para la pendiente.",
            ess_quality,
            results$ess[1], results$ess[2], results$ess[3],
            (results$n_iter - results$burn_in) * length(results$chains),
            results$burn_in,
            (results$ess[2] / ((results$n_iter - results$burn_in) * length(results$chains))) * 100
          ))
      ),
      
      div(class = "interpretation-box",
          h5(tags$b("🎯 Tasa de Aceptación:")),
          p(sprintf(
            "La tasa de aceptación promedio es %.1f%%, lo cual es %s. Las tasas individuales son: %s. 
            El rango óptimo es 20-50%% para algoritmos de Metropolis-Hastings.",
            avg_accept, accept_quality,
            paste(sprintf("Cadena %d: %.1f%%", 
                          1:length(results$acceptance_rates), 
                          results$acceptance_rates * 100), 
                  collapse = ", ")
          ))
      ),
      
      div(class = "interpretation-box",
          h5(tags$b("📈 Estimación Posterior de Tendencia:")),
          p(sprintf(
            "La media posterior de la pendiente es %.6f ± %.6f. El intervalo de credibilidad 
            90%% es [%.6f, %.6f] con amplitud %.6f. Esta estimación se basa en %d muestras 
            efectivas independientes de la distribución posterior.",
            slope_mean, slope_sd,
            slope_ci[1], slope_ci[2], slope_ci[2] - slope_ci[1],
            round(results$ess[2])
          ))
      ),
      
      div(class = "interpretation-box",
          h5(tags$b("🔍 Calidad General del Muestreo:")),
          p(sprintf(
            "Resumen: Convergencia %s | ESS %s | Aceptación %s. %s",
            convergence_quality, ess_quality, accept_quality,
            if(convergence_quality == "excelente" && ess_quality != "bajo") {
              "✓ Los resultados son confiables y pueden usarse para inferencia."
            } else if(convergence_quality == "buena") {
              "⚠ Los resultados son utilizables, pero considere ejecutar más iteraciones para mayor precisión."
            } else {
              "✗ Se recomienda ajustar parámetros y re-ejecutar antes de usar estos resultados."
            }
          ))
      )
    )
  })
  
  output$download_mcmc <- downloadHandler(
    filename = function() {
      paste0("mcmc_samples_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(mcmc_results())
      results <- mcmc_results()
      df <- as.data.frame(results$combined)
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  # ====== COMPARACIÓN ======
  output$comparison_plot <- renderPlotly({
    req(mc_results(), bayes_results(), mcmc_results())
    
    data <- data_reactive()
    mc <- mc_results()
    bayes <- bayes_results()
    mcmc <- mcmc_results()
    
    df <- data.frame(
      Método = factor(c("Monte Carlo", "Bayesiano", "MCMC"),
                      levels = c("Monte Carlo", "Bayesiano", "MCMC")),
      Estimación = c(mc$slope_mean, bayes$posterior_mean, 
                     mean(mcmc$combined[, "slope"])),
      IC_inferior = c(mc$slope_ci_lower, 
                      bayes$cred_lower,
                      quantile(mcmc$combined[, "slope"], 0.05)),
      IC_superior = c(mc$slope_ci_upper,
                      bayes$cred_upper,
                      quantile(mcmc$combined[, "slope"], 0.95))
    )
    
    p <- ggplot(df, aes(x = Método, y = Estimación, color = Método)) +
      geom_point(size = 5) +
      geom_errorbar(aes(ymin = IC_inferior, ymax = IC_superior), 
                    width = 0.3, size = 1.2) +
      scale_color_manual(values = c("Monte Carlo" = "#3498db",
                                    "Bayesiano" = "#27ae60",
                                    "MCMC" = "#9b59b6")) +
      labs(title = "Comparación de Estimaciones de Tendencia",
           subtitle = "Barras: Intervalos de Incertidumbre 90%",
           y = "Pendiente (tendencia anual)",
           x = "") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", size = 16))
    
    ggplotly(p) %>% layout(hovermode = 'x')
  })
  
  output$comp_mc_box <- renderInfoBox({
    req(mc_results())
    
    infoBox(
      "Monte Carlo",
      format(mc_results()$slope_mean, digits = 6),
      icon = icon("dice"),
      color = "blue",
      fill = TRUE
    )
  })
  
  output$comp_bayes_box <- renderInfoBox({
    req(bayes_results())
    
    infoBox(
      "Bayesiano",
      format(bayes_results()$posterior_mean, digits = 6),
      icon = icon("brain"),
      color = "green",
      fill = TRUE
    )
  })
  
  output$comp_mcmc_box <- renderInfoBox({
    req(mcmc_results())
    
    infoBox(
      "MCMC",
      format(mean(mcmc_results()$combined[, "slope"]), digits = 6),
      icon = icon("link"),
      color = "purple",
      fill = TRUE
    )
  })
  
  output$comparison_table <- renderDT({
    req(mc_results(), bayes_results(), mcmc_results())
    
    mc <- mc_results()
    bayes <- bayes_results()
    mcmc <- mcmc_results()
    
    df <- data.frame(
      Método = c("Monte Carlo", "Bayesiano Analítico", "MCMC"),
      Estimación = c(mc$slope_mean, bayes$posterior_mean, 
                     mean(mcmc$combined[, "slope"])),
      Desv_Estándar = c(mc$slope_sd, bayes$posterior_sd,
                        sd(mcmc$combined[, "slope"])),
      IC_90_inferior = c(mc$slope_ci_lower,
                         bayes$cred_lower,
                         quantile(mcmc$combined[, "slope"], 0.05)),
      IC_90_superior = c(mc$slope_ci_upper,
                         bayes$cred_upper,
                         quantile(mcmc$combined[, "slope"], 0.95))
    )
    
    datatable(df, 
              options = list(dom = 't', pageLength = 3),
              rownames = FALSE) %>%
      formatRound(2:5, digits = 8) %>%
      formatStyle('Método',
                  backgroundColor = styleEqual(
                    c('Monte Carlo', 'Bayesiano Analítico', 'MCMC'),
                    c('#ecf0f1', '#d5f4e6', '#e8daef')
                  ))
  })
  
  output$comparison_interpretation <- renderUI({
    req(mc_results(), bayes_results(), mcmc_results())
    
    mc <- mc_results()
    bayes <- bayes_results()
    mcmc <- mcmc_results()
    
    # Calcular diferencias
    estimates <- c(mc$slope_mean, bayes$posterior_mean, mean(mcmc$combined[, "slope"]))
    mean_estimate <- mean(estimates)
    max_diff <- max(abs(estimates - mean_estimate))
    rel_diff <- (max_diff / mean_estimate) * 100
    
    # Evaluar consistencia
    consistency <- if(rel_diff < 1) {
      "excelente"
    } else if(rel_diff < 5) {
      "buena"
    } else {
      "moderada"
    }
    
    # Amplitudes de intervalos
    widths <- c(
      mc$slope_ci_upper - mc$slope_ci_lower,
      bayes$cred_upper - bayes$cred_lower,
      quantile(mcmc$combined[, "slope"], 0.95) - quantile(mcmc$combined[, "slope"], 0.05)
    )
    
    # Verificar solapamiento
    intervals <- data.frame(
      lower = c(mc$slope_ci_lower, bayes$cred_lower, 
                quantile(mcmc$combined[, "slope"], 0.05)),
      upper = c(mc$slope_ci_upper, bayes$cred_upper,
                quantile(mcmc$combined[, "slope"], 0.95))
    )
    
    overlap_start <- max(intervals$lower)
    overlap_end <- min(intervals$upper)
    has_overlap <- overlap_start < overlap_end
    
    tagList(
      h4(icon("balance-scale"), " Análisis Comparativo"),
      hr(),
      
      div(class = "interpretation-box",
          h5(tags$b("🎯 Consistencia entre Métodos:")),
          p(sprintf(
            "Las tres estimaciones son: MC=%.6f, Bayesiano=%.6f, MCMC=%.6f. 
            La media es %.6f con una desviación máxima de %.6f (%.2f%% relativo). 
            La consistencia es %s.",
            estimates[1], estimates[2], estimates[3],
            mean_estimate, max_diff, rel_diff, consistency
          )),
          p(sprintf(
            "%s",
            if(consistency == "excelente") {
              "✓ Las estimaciones son prácticamente idénticas, lo que indica robustez de los resultados."
            } else if(consistency == "buena") {
              "✓ Las estimaciones son muy similares, demostrando convergencia entre metodologías."
            } else {
              "⚠ Existe cierta variación entre métodos. Esto puede deberse a diferentes supuestos o configuraciones."
            }
          ))
      ),
      
      div(class = "interpretation-box",
          h5(tags$b("📏 Amplitud de Intervalos:")),
          p(sprintf(
            "Amplitudes de IC/Credibilidad 90%%: MC=%.6f, Bayesiano=%.6f, MCMC=%.6f. 
            El método %s produce el intervalo más estrecho (mayor precisión), 
            mientras que %s produce el más amplio.",
            widths[1], widths[2], widths[3],
            c("Monte Carlo", "Bayesiano", "MCMC")[which.min(widths)],
            c("Monte Carlo", "Bayesiano", "MCMC")[which.max(widths)]
          )),
          p(sprintf(
            "La razón entre intervalos más ancho/estrecho es %.2f. %s",
            max(widths) / min(widths),
            if(max(widths) / min(widths) < 1.5) {
              "Los métodos concuerdan en el nivel de incertidumbre."
            } else {
              "Existe discrepancia en la cuantificación de incertidumbre entre métodos."
            }
          ))
      ),
      
      div(class = "interpretation-box",
          h5(tags$b("🔄 Solapamiento de Intervalos:")),
          p(sprintf(
            "%s Los tres intervalos %s. %s",
            if(has_overlap) "✓" else "✗",
            if(has_overlap) {
              sprintf("se solapan en el rango [%.6f, %.6f]", overlap_start, overlap_end)
            } else {
              "NO se solapan completamente"
            },
            if(has_overlap) {
              "Esto confirma que los métodos están identificando la misma región de valores plausibles."
            } else {
              "Esto sugiere posible inconsistencia en supuestos o necesidad de revisar configuraciones."
            }
          ))
      ),
      
      div(class = "interpretation-box",
          h5(tags$b("🏆 Recomendación:")),
          p(sprintf(
            "Para este análisis, se recomienda usar el método %s porque: %s",
            if(consistency == "excelente") {
              "MCMC"
            } else {
              c("Monte Carlo", "Bayesiano", "MCMC")[which.min(widths)]
            },
            if(consistency == "excelente") {
              "proporciona la estimación posterior completa con diagnósticos de convergencia robustos."
            } else {
              "ofrece el balance óptimo entre precisión y confiabilidad en las estimaciones."
            }
          )),
          p(HTML(
            "<b>Nota:</b> La concordancia entre métodos valida las conclusiones. 
            Diferencias menores son esperables debido a distintos marcos conceptuales 
            (frecuentista vs. bayesiano) y aleatorización."
          ))
      )
    )
  })
}

# ============================================================================
# EJECUTAR APLICACIÓN
# ============================================================================

shinyApp(ui = ui, server = server)
