# Load required libraries
library(shiny) #package for this app (R shiny app)
library(dplyr) #package for processing data table
library(epitools) #odds ratio
library(genpwr) #power of test

OR_fun <- function(N_con, N_Gcon, N_cas, N_Gcas, OR=2, Alpha=0.05){
  Geno <- c(N_con, # No. of control
            N_Gcon, # No. of genotype carrier in control
            N_cas, # No. of case
            N_Gcas) # No. of genotype carrier in case
  # Create contingency table
  cases_exposed <- Geno[4]
  cases_unexposed <- Geno[3]-Geno[4]
  controls_exposed <- Geno[2]
  controls_unexposed <- Geno[1]-Geno[2]

  mymatrix <- matrix(c(cases_exposed, cases_unexposed,
                       controls_exposed, controls_unexposed),
                     nrow = 2,
                     dimnames = list(Exposure = c("Carrier", "Non-carrier"),
                                     Outcome = c("Case", "Control")))

  # Calculate OR using epitools package
  oddsratio_result <- oddsratio(mymatrix, method = "wald")


  MAF <- round(Geno[2]/Geno[1], digits = 4)

  power_result <- genpwr.calc(
    calc = "power",
    model = "logistic",
    N = Geno[1]+Geno[3],                # total sample size (cases + controls)
    Case.Rate = Geno[3]/(Geno[1]+Geno[3]),         # proportion of cases in the sample
    MAF = MAF,            # G84E carrier rate (0.34%)
    OR = OR,                  # expected odds ratio
    Alpha = Alpha,            # significance level
    True.Model = "Dominant", # specify genetic model if known
    Test.Model = "Dominant"  # model used for testing
  )

  ResultTable <- data.frame(
    `OR` = paste(as.character(round(oddsratio_result$measure[2,1], digits = 2)),
                    "(",
                    as.character(round(oddsratio_result$measure[2,2], digits = 2)),
                    ", ",
                    as.character(round(oddsratio_result$measure[2,3], digits = 2)),
                    ")", sep = ""),
    `p-value` = ifelse(oddsratio_result$p.value[2,3]<0.001, "<0.001*",
                       ifelse(oddsratio_result$p.value[2,3]<0.05, paste(round(oddsratio_result$p.value[2,3], digits = 4),
                                                                        "*", sep = ""),
                              round(oddsratio_result$p.value[2,2], digits = 2))),
    `Power` = round(power_result$Power_at_Alpha_0.05, digits = 2)
  )

  colnames(ResultTable) <- c("OR (95%CI)", "p-value", "Power")
  return(
    list(
      TF.Table = as.data.frame(oddsratio_result$data),
      OR.Table = ResultTable
    )
  )
}


# Define UI for the Shiny app
ui <- fluidPage(

  titlePanel("Odds Ratio Calculator"),

  sidebarLayout(

    sidebarPanel(
      textInput("N_con", "No. of control", ""),
      textInput("N_Gcon", "No. of mutation carriers in control", ""),
      textInput("N_cas", "No. of case", ""),
      textInput("N_Gcas", "No. of mutation carriers in case", ""),
      actionButton("Go", "Run")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Table", dataTableOutput("TFtable")),
        tabPanel("Odds Ratio", dataTableOutput("OR"))
      )
    )
  )
)

# Define server logic for the Shiny app
server <- function(input, output, session) {

  inputFile <- reactiveValues()

  observeEvent(input$Go, {
    req(input$N_con,
        input$N_Gcon,
        input$N_cas,
        input$N_Gcas)

    oddsratio_result <- OR_fun(as.numeric(input$N_con),
                               as.numeric(input$N_Gcon),
                               as.numeric(input$N_cas),
                               as.numeric(input$N_Gcas))

    inputFile$data <- oddsratio_result[[1]]
    inputFile$OR <- oddsratio_result[[2]]

  })

  output$TFtable <- renderDataTable(
    inputFile$data
  )

  output$OR <- renderDataTable(
    inputFile$OR
  )

}

# Run the application
shinyApp(ui, server)
