library(shiny)
library(markdown)
library(shinyjs)
library(shinyWidgets)

shinyUI(fluidPage(
  title = "Computer Experiments Sample Size", 
  # tags$head(includeScript("google-analytics.js")),
  titlePanel(title = div("Sample Size Calculations for Computer Experiments",
                         style = "font-weight: 700; line-height: 2; font-size: 120%;
                        color: #000033; background-color: lightgrey;")),
  
  
  useShinyjs(),
  
  
  tags$head(tags$style(
    "body { word-wrap: break-word; }")),
  
  tags$style(HTML("
                  .tabbable > .nav > li > a                  
                  {background-color: black;  color:white}"
  )),
  
  
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  ),
  
  
  tabsetPanel(
    
    tabPanel("Sample Size",
             br(),
             tags$div(
               tags$ul(
                 tags$li(h5(strong("Read the 'About' tab for the definitions of the 
	                           different parameters.")))
               )
             ),
             tags$div(
               tags$ul(
                 tags$li(h5(strong("Use the 'Sample Path Plots' tab to choose appropriate correlation family and hyperparameters.")))
               )
             ),
             tags$div(
               tags$ul(
                 tags$li(h5(strong("Choose your configuration, then press the 'Go' button. 
	                           Scroll down for the results.")))
               )
             ),
             fluidRow(
               
               column(width = 3,
                      br(),     
                      selectInput("corr_family", "Select correlation family", c("Gaussian", "Power Exponential", "Matern")),
                      sliderInput(inputId = "d", 
                                  label = "How many inputs do you have in your computer model?",
                                  value=4, min=1, max=20),
                      br(),
                      radioButtons(inputId = "goal", label = "What would you like to do?",
                                   c("Find the sample size required for a given level of Root Average Unexplained Variation" = "var_to_n",
                                     "Calculate the Root Average Unexplained Variation for a given sample size" = "n_to_var")
                      )
               ), 
               
               column(width = 3,
                      br(),
                      uiOutput("action_thetas")
               ),
               
               column(width = 3,
                      br(),
                      uiOutput("action_p")
               ),
               
               column(width = 3,
                      br(),
                      br(),
                      actionButton("go", "Go",
                                   class = "btn-primary")
               )
               
             ), 
             
             fluidRow(
               column(3,""),
               column(7,
                      br(),
                      textOutput("mytext"), 
                      tags$head(tags$style("#mytext{color: black;
                                 font-size: 20px;
                                 font-style: italic;
                                 }")
                      ),
                      plotOutput("graph")
                      
               )
             )
    ),
    
    
    tabPanel("Sample Path Plots",
             
             fluidRow(
               
               column(2,
                      br(),
                      selectInput("corr_family2", "Select correlation family", c("Gaussian", "Power Exponential", "Matern")),
                      numericInput("samp_theta", HTML("Correlation length (&theta;)"), min=.1, step=.1, value=1),
                      uiOutput("p_or_nu"),
                      numericInput("samp_n", "No. of realizations", min=1, step=1, value=7)
               ),
               
               column(8,
                      plotOutput("samp_path")
               )
             )
             
    ),
    
    
    tabPanel("Robust Sample Size",
             tags$div(
               tags$ul(
                 tags$li(h5(strong("Independent uniform distributions are assigned to the correlation length parameters over the specified range.")))
               )
             ),
             
             tags$div(
               tags$ul(
                 tags$li(h5(strong("The selected sample size is the sample quantile (provided by the user) of the corresponding Monte Carlo sample.")))
               )
             ),
             
             tags$div(
               tags$ul(
                 tags$li(h5(strong("Choose your configuration, then press the 'Go' button. 
                                   Scroll down for the results.")))
               )
             ),
             fluidRow(
               column(3,
                      br(),     
                      selectInput("corr_family3", "Select correlation family", c("Gaussian", "Matern")),
                      sliderInput(inputId = "d2", 
                                  label = "How many inputs do you have in your computer model?",
                                  value=4, min=1, max=20),
                      numericInput("var_target2", label="Target RAUV", 
                                   value = 0.05, min = .005, max = .995, step = .005),
                      br(),
                      uiOutput("action_nu")
               ),
               
               column(3,
                      br(),
                      uiOutput("action_thetas_low")
               ),
               
               column(3,
                      br(),
                      uiOutput("action_thetas_high")
               ),
               
               column(3,
                      br(),
                      numericInput("quantile", label="Cutoff Quantile", 
                                   value = 0.90, min = 0, max = 1, step = .05),
                      numericInput("MC_size", label="Monte Carlo Sample Size", 
                                   value = 5000, min = 1000, max = 10000, step = 500),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      actionButton("go2", "Go", class = "btn-primary")
               )
             ),
             fluidRow(
               column(3),
               column(6,
                      br(),
                      br(),
                      plotOutput("histogram")),
               column(3)
             )
             
    ),
    
    
    tabPanel("About",
             fluidRow(
               includeHTML("About.html")
             )
    )
  )
)
)




#  )
#)




