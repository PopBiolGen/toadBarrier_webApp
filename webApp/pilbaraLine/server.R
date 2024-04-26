library(shiny)
library(leaflet)
library(RColorBrewer)
library(scales)
library(lattice)
library(dplyr)
library(sf)

source("functions.R")

wb<-st_read("mapData/waterbodies.shp")
wb<-select(wb, Property, Area, Name, IPA, Type, Purpose, Conversion, Maintenanc)
wb$Purpose<-factor(wb$Purpose, levels = c(levels(wb$Purpose), "Unmanaged"))
wb$Property<-as.character(wb$Property)
wb$Property[is.na(wb$Property)]<-"Unclassified"
wb.xy<-st_coordinates(wb)
ppt<-st_read("mapData/propertiesWGS.shp")
colrs<-terrain.colors(12, alpha = 0.5)


load("modelData/IPAF.RData")  # probOutIPAF
load("modelData/IPAT.RData")  # probOutIPAT

# Define server logic required to plot our map
function(input, output, session) {
  
  output$map <- renderLeaflet(
    leaflet(data=ppt, options=leafletCRS(crsClass = "L.CRS.EPSG4326")) %>%
      addTiles() %>%
      addScaleBar(position = "topleft") %>%
      setView(lng = 120.943, lat = -19.336, zoom = 7) %>%
      addPolygons(label = ppt$PROPERTY_N)
  )
  
  # returns subset to be managed
  ssFun<-reactive(modify.plb(wb, 
                          future.irrigation = input$fIrr, 
                          IPA = input$fIPA, 
                          Dwelling = input$fDwell)
  )
  # returns appropriate probability dataset
  pDat<-reactive(
    if (input$fIPA) probOutIPAT else probOutIPAF
  )
  
  # A reactive expression that returns the set of wbs that are
  # # in bounds right now
  wbInBounds <- reactive({
    if (is.null(input$map_bounds))
      return(NA)
    bounds <- input$map_bounds
    latRng <- range(bounds$north, bounds$south)
    lngRng <- range(bounds$east, bounds$west)
    cenX<-mean(lngRng) # find center of the map
    cenY<-mean(latRng)
    
    # Find susbet of data represented on the map
    ss<-ssFun()
    currWB<-subset(wb[ss,],
                   wb.xy[ss,"Y"] >= latRng[1] & wb.xy[ss,"Y"] <= latRng[2] &
                     wb.xy[ss,"X"] >= lngRng[1] & wb.xy[ss,"X"] <= lngRng[2])
    nWbTemp<-nrow(currWB)
    
    # Find mamangement point closest to the centre
    tPoints<-test.points(wb.xy[ss,])  # generate Darren's test points
    sqDist<-(tPoints[,"X"]-cenX)^2+(tPoints[,"Y"]-cenY)^2 # squared distance from centre
    minInd<-which(sqDist==min(sqDist))[1] # an index for which prop column to use
    # Use appropriate model results
    pDatList<-pDat()
    xTest<-c(pDatList$nWb[,minInd], 200) # assume 200 waterpoints will always generate 0 breach prob
    yTest<-c(pDatList$pSucc[,minInd], 0)
    # estimate probability of success
    if (nWbTemp > 200) {probSucc <- 1} else {
      probSucc<-1-approx(x = xTest,
                          y = yTest,
                          xout = nWbTemp)$y
    }
    
    # get current waterbodies and report table
    currTab<-data.frame(N = as.integer(nWbTemp),
                   Conv = sum(currWB$Conversion)/1e6, #Maintenance and conversion costs, in $M
                   Maint = sum(currWB$Maintenanc)/1e6,
                   pSucc = probSucc)
    currTab
  })
  
 # This observer is responsible for maintaining the circles and legend,
#  according to the variables the user has chosen to map to color and size.
  observe({
    # plot colours
    colourBy <- "Purpose"
    # # plotted subset
    ss<-ssFun()
   # manwb<-wb[ss,] # managed waterbodies
   # manwb[[colourBy]]<-as.factor(as.character(manwb[[colourBy]])) # re-factor
    #unmanwb<-wb[!ss,] # unmanaged waterbodies
    tPoints<-test.points(wb.xy[ss,])  # generate Darren's test points
    
    colourData <- wb[[colourBy]]
    colourData[!ss]<-"Unmanaged"
    pal <- colorFactor("viridis", levels(colourData))
    
    leafletProxy("map") %>%
      clearGroup(group = "wbs") %>%
      addMarkers(lng = tPoints[,"X"], 
                 lat = tPoints[,"Y"], 
                 label = rownames(tPoints),
                 group = "wbs") %>%
      addCircleMarkers(wb,
                       lng = wb.xy[,1],
                       lat = wb.xy[,2],
                       radius=7,
                       color = pal(colourData),
                       label = wb$Name,
                       group = "wbs") %>%
      addLegend("topleft", pal=pal, values=colourData, title=colourBy,
                layerId="colorLegend")
       
      output$nWB<-renderTable(wbInBounds(), digits=2)
  })
  
}
