#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(readr)
library(dplyr)
library(plotly)
library(shinyWidgets)

# ind.coords.T <- read.table("data/Thompson_indCoords_neutral_cands.txt", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
# head(ind.coords.T)
# centroid.T <- read.table("data/Thompson_centroids_neutral_cands.txt", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
# head(centroid.T)
# percent.T <- read.table("data/Thompson_percent_neutral_cands.txt", header = TRUE, stringsAsFactors = FALSE)
# head(percent.T)
# 
# 
# ind.coords <- read.table("/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/11-DAPC/Coho_DAPC/data/noThompson_indCoords_neutral_cands.txt", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
# head(ind.coords)
# centroid <- read.table("/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/11-DAPC/Coho_DAPC/data/noThompson_centroids_neutral_cands.txt", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
# head(centroid)
# percent <- read.table("/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/11-DAPC/Coho_DAPC/data/noThompson_percent_neutral_cands.txt", header = TRUE, stringsAsFactors = FALSE)
# head(percent)


## using the ind.coords dataframe that is already loaded 

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("DAPC Scatterplot"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            radioButtons("region",
                        "Select Region",
                        choices = c("BC", "Thompson"), selected = "BC"),
            radioButtons("marker",
                         "Select SNP Type", 
                         choices = c("neutral", "GEA", "RDA"), selected = "neutral"),
            actionButton("marker_info", "Click here for SNP info", style = "padding:4px; font-size:80%"),
            
            selectInput("x_axis",
                        "Select First Axis",
                        choices = paste0("Axis", seq(1,6,1)), selected = "Axis1"),
            selectInput("y_axis",
                        "Select Second Axis",
                        choices = paste0("Axis", seq(1,6,1)), selected = "Axis2")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           textOutput("description"),
           plotly::plotlyOutput('plot_dapc_scatter', height = 800, width = 800)
        )
     )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$description <- renderText({
        # generate bins based on input$bins from ui.R
        paste0('Selected Region: ', input$region)
    })
    
    
   ind.coords <- reactive({
     if(input$region == "Thompson") {
       i <- read.table("data/Thompson_indCoords_neutral_cands.txt", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
     } else {
       i <- read.table("data/noThompson_indCoords_neutral_cands.txt", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
     }
     return(i)
   })

   centroid <- reactive({
     if(input$region == "Thompson") {
       j <- read.table("data/Thompson_centroids_neutral_cands.txt", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
     } else {
       j <- read.table("data/noThompson_centroids_neutral_cands.txt", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
     }
     return(j)
   })


   percent <- reactive({
     if(input$region == "Thompson") {
       k <- read.table("data/Thompson_percent_neutral_cands.txt", header = TRUE, stringsAsFactors = FALSE)
     } else {
       k <- read.table("data/noThompson_percent_neutral_cands.txt", header = TRUE, stringsAsFactors = FALSE)
     }
     return(k)
   })
   
    
    
    indcoords_subset <- reactive({
        a <- dplyr::filter(ind.coords(), MARKER == input$marker)
        return(a)
    })
        
    
    centroid_subset <- reactive({
        b <- dplyr::filter(centroid(), MARKER == input$marker)
        return(b)
    })
    
    percent_subset <- reactive({
        p <- dplyr::filter(percent(), MARKER == input$marker)
        return(p)
    })
    
    observeEvent(input$marker_info, {
      showModal(modalDialog(title = "Number of SNPs in each dataset",
                HTML("BC Region: <br><b>Neutral =</b> 22,992 SNPs <br><b>GEA (RDA + LFMM) Outliers =</b> 1,117 SNPs <br><b>RDA Outliers = </b>450 SNPs <br><br> Thompson Region: <br><b>Neutral =</b> 9,869 SNPs <br><b>GEA (RDA + LFMM) Outliers = </b>124 SNPs <br><b>RDA Outliers = </b>97 SNPs")))
    })
    
    output$plot_dapc_scatter <- plotly::renderPlotly({
        plot_ly(type = "scatter", mode = "markers") %>%
            add_trace(data = indcoords_subset(), x = ~get(input$x_axis), y = ~get(input$y_axis),
                      text = indcoords_subset()$pop_cu,
                      hoverinfo = "text",
                      opacity = 0.2,
                      marker = list(symbol = "circle", size = 5, color = indcoords_subset()$color),
                      hoverlabel = list(font = list(size = 5)),
                      showlegend = F) %>%
            add_trace(data = centroid_subset(), x = ~get(input$x_axis), y = ~get(input$y_axis),
                      text = centroid_subset()$pop_cu,
                      hoverinfo = "text",
                      hoverlabel = list(font = list(size = 10)),
                      marker = list(symbol = "diamond",
                                    size = 8,
                                    shape = 23, 
                                    color = centroid_subset()$color,
                                    line = list(color = 'black', width = 1)),
                      showlegend = F) %>%
            layout(xaxis = list(title = paste(input$x_axis, " (", round(percent_subset()$percent[which(percent_subset()$Axis == input$x_axis)], 2), "%)", sep = ""), zerolinewidth = 0.01, zerolinecolor = "darkgrey"),
                   yaxis = list(title = paste(input$y_axis, " (", round(percent_subset()$percent[which(percent_subset()$Axis == input$y_axis)], 2), "%)", sep = ""), zerolinewidth = 0.01, zerolinecolor = "darkgrey")
            )
    })
    
}






# Run the application 
shinyApp(ui = ui, server = server)
