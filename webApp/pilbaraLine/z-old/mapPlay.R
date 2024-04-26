library(leaflet)
library(RColorBrewer)
library(scales)
library(lattice)
library(dplyr)
library(sf)


wb<-st_read("mapData/waterbodiesWGS.shp")
wb.xy<-st_coordinates(wb)
ppt<-st_read("mapData/propertiesWGS.shp")





leaflet(data=ppt, options=leafletCRS(crsClass = "L.CRS.EPSG4326")) %>%
  addTiles() %>%
  setView(lng = 120.943, lat = -19.336, zoom = 7) %>%
  # addTiles(
  #   urlTemplate = "//{s}.tiles.mapbox.com/v3/benflips.cjyteb8xi03dd2ntbh8msq19l-2hr6i/{z}/{x}/{y}.png",
  #   attribution = 'Maps by <a href="http://www.mapbox.com/">Mapbox</a>'
  # ) %>%
  addPolygons(label = ppt$PROPERTY_N) %>%
  addCircleMarkers(wb,
                   lng = wb.xy[,1], 
                   lat = wb.xy[,2], 
                   radius=7, 
                   color = as.numeric(wb$FEATTYPE),
                   label = wb$NAME)