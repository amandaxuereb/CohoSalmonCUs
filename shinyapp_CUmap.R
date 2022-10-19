#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(leaflet)
#library(dplyr)
library(tidyverse)
library(leaflet.extras)
library(htmltools)
library(htmlwidgets)
library(readr)
library(data.table)
library(DT)
library(RColorBrewer)


# ## read in sites 
# sites <- read_tsv("data/ALL_SITES_UPDATED.txt", col_names = FALSE)
# head(sites)
# colnames(sites) <- c("SITE", "FULL", "REGION", "CU", "LATITUDE", "LONGITUDE")
# sites$LONGITUDE
# 
# # for colour coding: 
# pop.cus.col <- read_tsv("data/pop_cus_colours.txt")
# nrow(pop.cus.col)
# pop.cus.col$POP[which(pop.cus.col$POP == "LIO")] <- "LYO"
# pop.cus.col <- dplyr::filter(pop.cus.col, POP %in% sites$SITE)
# nrow(pop.cus.col)
# head(pop.cus.col)
# 
# length(unique(pop.cus.col$color))
# cus.col <- select(pop.cus.col, -pop_cu)
# head(cus.col)
# 
# head(sites)
# colnames(sites)[1] <- "POP"
# sites.cu.tb <- left_join(sites, cus.col, by = "POP")
# nrow(sites.cu.tb)
# head(sites.cu.tb)
# 
# unique(sites.cu.tb$CU.x)
# unique(sites.cu.tb$color)
# 
# pal <- colorFactor(palette = unique(sites.cu.tb$color), levels = unique(sites.cu.tb$CU.x))
# 
# 
# ## add in the CU names to table: 
# cu.names <- read_tsv("data/cu.data.txt", col_names = TRUE)
# head(cu.names)
# 
# colnames(cu.names)[1] <- "CU.x"
# 
# sites.cu.tb.names <- left_join(sites.cu.tb, cu.names, by = "CU.x")
# head(sites.cu.tb.names)
# 
# sites.cu.tb.filt <- sites.cu.tb.names[,c(1:6,8,10)]
# head(sites.cu.tb.filt)
# write.table(sites.cu.tb.filt, file = "./data/FULL_DATA_TABLE.txt", col.names = T, row.names = F, quote = F, sep = "\t")

##### START HERE: 


sites.cu.tb <- read_tsv("./data/FULL_DATA_TABLE.txt")
head(sites.cu.tb)
nrow(sites.cu.tb) ## 146 sampling locations

## Add CU polygons 
# use a SpatialPolygonsDataFrame 
library(rgeos)
library(rgdal)

cu.map <- readOGR("data/Coho_Salmon_CU_simplified/", layer = "cu.simp")
class(cu.map)

cu.data <- read_tsv("data/cu.data.txt")
head(cu.data)

cu.map@data <- cu.data

color.list <- read.table("data/color_list.txt", header = F, stringsAsFactors = F)
color.list
nrow(color.list)
head(color.list)

color.listb <- color.list
color.listb$V1[which(color.listb$V1 == "#FFFFFF")] <- 'black'


cu.pal <- colorFactor(palette = color.list$V1, levels = cu.map@data$Full_CU_IN)

cu.palb <- colorFactor(palette = color.listb$V1, levels = cu.map@data$Full_CU_IN)




#### connect the map and the table in the app: 

### ---- UI 
ui <- fluidPage(
  absolutePanel(bottom = 10, left = 10, actionButton("show_help", "Help")
  ),
  
  leaflet::leafletOutput('x2'),
  DT::DTOutput('x1')
)


### ---- SERVER 

server <- function(input, output, session) {
  
  observeEvent(input$show_help, {
    showModal(modalDialog(title = "Help", HTML("The coloured polygons show the CU boundaries and the points (sampling locations) are coloured according to the CU they are assigned to. Click on the polygon to see a popup that shows the CU name and number. Click on a point to show a popup with information about that sampling location. Hovering your cursor over the points also shows a label with the three-letter code for the sampling location. <br><br> Use the magnifying glass at the top left to search for locations on the map (using the 3-letter code). This willzoom to the searched location. Use the button with the 4 small arrows pointing toward each other (just above the magnifying glass) to return to the original map extent. You can also use the 'layer' button in the top right corner to toggle the CU polygons on and off, or select/deselect each of the two regional groups (BC and Thompson). <br><br> The table displays site information, which can be sorted by clicking on column headings. Clicking on a row will zoom the map to that sampling location. You can also search the table by typing in the Search box.")))
  })
  
  sites.data <- reactive({ 
    sites.cu.tb %>%
      select(-color)
  })
  
  
  output$x1 <- DT::renderDT({
    DT::datatable(sites.data(), colnames = c("SITE", "FULL SITE NAME", "n", "CU", "CU NAME", "LATITUDE", "LONGITUDE", "REGION"), selection = "single", options = list(stateSave = TRUE, scrollY = "200px", pageLength = 146, lengthMenu = 146, lengthChange = FALSE))
  })
  
  
  output$x2 <- leaflet::renderLeaflet({
    sites.data() %>%
      leaflet() %>%
      addProviderTiles("CartoDB") %>%
      fitBounds(lng1 = -135, lat1 = 47, lng2 = -118, lat2 = 58) %>%
      addResetMapButton() %>%
      addCircleMarkers(data = sites.data(),
                       radius = 6, stroke = T, color = "black", weight = 1,
                       opacity = 1, fillColor = ~cu.pal(CU), fillOpacity = 1,
                       popup = ~paste0("<b>", POP, "</b>", "<br/>", FULL, "<br/>", "Region: ", REGION, "<br/>", "CU: ", CU, ", ", CU_name, "<br/>", "n = ", n),
                       label = ~POP,
                       group = ~REGION) %>%
      
      addMarkers(data = sites.data(), lng = ~LONGITUDE, lat = ~LATITUDE, label = ~POP, group = ~POP, icon = makeIcon( 
        iconUrl = "http://leafletjs.com/examples/custom-icons/leaf-green.png",
        iconWidth = 1, iconHeight = 1
      )) %>%
      
      addSearchFeatures(targetGroups = sites.data()$POP, options = searchFeaturesOptions(zoom = 8, openPopup = TRUE)) %>%
      
      addPolygons(data = cu.map, weight = 0.5,
                  fillColor = ~cu.pal(Full_CU_IN),
                  fillOpacity = 0.3,
                  color = ~cu.palb(Full_CU_IN),
                  popup = ~paste0("<b>", Full_CU_IN, "</b>", "<br/>", CU_name, "<br/>"),
                  highlight = highlightOptions(weight = 1,
                                               color = "black",
                                               bringToFront = F, opacity = 1,
                                               fillOpacity = 0.2,
                                               sendToBack = T),
                  group = "CU") %>%
      
      addLayersControl(overlayGroups = c("BC", "Thompson", "CU"))
  })
  
  
  observeEvent(input$x1_rows_selected, {
    row_selected = sites.data()[input$x1_rows_selected,]
    proxy <- leafletProxy('x2')
    #print(row_selected)
    proxy %>%
      setView(lng = row_selected$LONGITUDE,
              lat = row_selected$LATITUDE,
              zoom = 8) %>% 
      addPopups(lng = row_selected$LONGITUDE,
                lat = row_selected$LATITUDE,
                popup = row_selected$POP,
                options = popupOptions(closeButton = TRUE)
      ) %>%
      addCircleMarkers(data = sites.data()[input$x1_rows_selected,],
                       radius = 6, stroke = T, color = "black", weight = 1,
                       opacity = 1, fillColor = ~cu.pal(CU), fillOpacity = 1,
                       popup = ~paste0("<b>", POP, "</b>", "<br/>", FULL, "<br/>", "Region: ", REGION, "<br/>", "CU: ", CU, "<br/>", "n = ", n),
                       label = ~POP,
                       group = ~REGION)
    
  })  
  
}


shinyApp(ui = ui, server = server)

