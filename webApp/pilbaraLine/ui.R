library(shiny)
library(leaflet)

# Variables for drop-down menus
#cvars<-c("Purpose" = "Purpose", "Cost" = "cost1") # colouring
ivars<-c("Happens" = TRUE, "Doesn't happen" = FALSE) # future irrigation
ipavars<-c("Happens" = TRUE, "Doesn't happen" = FALSE) # ipa water management
dwellvars<-c("Dwellings managed" = TRUE, "Dwellings not managed" = FALSE) # dwelling management

navbarPage("Toad Barrier", id="nav",
           
  tabPanel("Interactive map",
    div(class="outer",
                        
    tags$head(
      # Include our custom CSS
      includeCSS("styles.css"),
      includeScript("gomap.js")
    ),
    
    # If not using custom CSS, set height of leafletOutput to a number instead of percent
    leafletOutput("map", width="100%", height="100%"),
    
    # Shiny versions prior to 0.11 should use class = "modal" instead.
    absolutePanel(id = "outputs", class = "panel panel-default", fixed = TRUE,
        draggable = TRUE, top = "auto", left = "auto", right = 10, bottom = 20,
        width = 280, height = "auto",
        
        h3("Assumptions"),
        selectInput("fIrr", "Future Irrigation", ivars),
        selectInput("fIPA", "IPA Plan", ipavars),
        selectInput("fDwell", "Dwellings", dwellvars),
        #selectInput("colour", "Colour", cvars),
        hr(),
        h3("This area managed"),
        tableOutput("nWB")
        #plotOutput("scatterCollegeIncome", height = 250)
    ),
    
    tags$div(id="cite",
             'Data from Southwell et al. 2017, ', tags$em('J Appl. Ecol.'), '54:1365.'
    )
    )
  ),
           # 
           # tabPanel("Data explorer",
           #          fluidRow(
           #            column(3,
           #                   selectInput("states", "States", c("All states"="", structure(state.abb, names=state.name), "Washington, DC"="DC"), multiple=TRUE)
           #            ),
           #            column(3,
           #                   conditionalPanel("input.states",
           #                                    selectInput("cities", "Cities", c("All cities"=""), multiple=TRUE)
           #                   )
           #            ),
           #            column(3,
           #                   conditionalPanel("input.states",
           #                                    selectInput("zipcodes", "Zipcodes", c("All zipcodes"=""), multiple=TRUE)
           #                   )
           #            )
           #          ),
           #          fluidRow(
           #            column(1,
           #                   numericInput("minScore", "Min score", min=0, max=100, value=0)
           #            ),
           #            column(1,
           #                   numericInput("maxScore", "Max score", min=0, max=100, value=100)
           #            )
           #          ),
           #          hr(),
           #          DT::dataTableOutput("ziptable")
           # ),
           
           conditionalPanel("false", icon("crosshair"))
)

