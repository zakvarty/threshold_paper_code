
# Create csv earthquake catalogue from JSON file

library("sp")
library("jsonlite")
library("readr")
library("lubridate")
library("matrixStats")
source('../00_src/catalogue_creation/get_catalogue.R')
source('../00_src/catalogue_creation/read_JSON_cat.R')

download_new_json = FALSE
if (download_new_json) write_groningen_json(file_path = "./output/")

JSONpath <- './output/2020-05-14_14-46-48_cat.json'
CSVpath <- paste0(substr(JSONpath,start = 1,stop = 32),'.csv')

## Set up
latlong = "+init=epsg:4326"
google = "+init=epsg:3857"
RD = "+init=epsg:28992"
gfo <- read.table("./Groningen_Field_outline.csv", sep=",", header=TRUE)
pol <- Polygon(coords=as.matrix(gfo[,1:2]), hole=FALSE)
pols <- Polygons(list(pol), ID="1")
SpPols <- SpatialPolygons(list(pols), proj4string=CRS(as.character(NA)))
SPgfo <- SpatialPolygonsDataFrame(SpPols, data=data.frame(name="GFO",from="SJB" ), match.ID = TRUE)
rm(pol,pols,SpPols)

proj4string(SPgfo) <- CRS(RD)

M1 <- data.frame(x=c(243584, 231390, 232874, 240358,243584), y=c(585466, 599000, 578950, 570627,585466))
M2 <- data.frame(x=c(250671, 254639, 234923,231390,243584,250671),y=c(584575, 615566, 612228,599000,585466.5,584575))
M3 <- data.frame(x=c(264250,250671,243584,240358,264250),y=c(565767,584575,585467,570627,565767))
M4 <- data.frame(x=c(271325,254639,250671,264250,271325), y=c(587862,615566,584575,565767,587862))
pol1 <- Polygon(coords=as.matrix(cbind(M1$x,M1$y)), hole=FALSE)
pol2 <- Polygon(coords=as.matrix(cbind(M2$x,M2$y)), hole=FALSE)
pol3 <- Polygon(coords=as.matrix(cbind(M3$x,M3$y)), hole=FALSE)
pol4 <- Polygon(coords=as.matrix(cbind(M4$x,M4$y)), hole=FALSE)
pols <- Polygons(list(pol1,pol2,pol3, pol4), ID=c("1"))
APols <- SpatialPolygons(list(pols), proj4string=CRS(RD))
APareas <- SpatialPolygonsDataFrame(APols, data=data.frame(name="Areas",from="Stijn" ), match.ID = TRUE)
rm(pols,APols, pol1, pol2, pol3, pol4)


AreaNames <- c("All regions", "Eemskanaal", "North West", "East", "South East")
AreaNamesNederlands <- c("Totale veld", "Eemskanaal", "Noord-West", "Oost", "Zuid-Oost")

gfo_proj <- spTransform(SPgfo, CRS("+proj=longlat +datum=WGS84"))
APareas_proj <- spTransform(APareas, CRS("+proj=longlat +datum=WGS84"))

GFObuffer <- read.table(file="./GFObuffer1000.csv", header=TRUE,sep=",")
gfoPBSbuffer <- data.frame(PID=rep(1,nrow(GFObuffer)),SID=rep(1,nrow(GFObuffer)),POS=1:nrow(GFObuffer),X=GFObuffer$X,Y=GFObuffer$Y)

## Make and save catalogue
gf_cat <- read_JSON_cat(Tpath= JSONpath,M1=M1,M2=M2,M3=M3,M4=M4,gfoPBSbuffer=gfoPBSbuffer)
readr::write_csv(gf_cat,path = CSVpath)

## Plot to check all okay
if(FALSE){
  # PLot locations
 plot(SPgfo)
 points(gf_cat$RD_X, gf_cat$RD_Y)
  # Plot times and magnitudes
 plot(as.Date(gf_cat$date),gf_cat$mag)
}

## Clear workspace
rm( "APareas",
    "APareas_proj",
    "AreaNames",
    "AreaNamesNederlands",
    "gf_cat",
    "gfo",
    "gfo_proj",
    "GFObuffer",
    "gfoPBSbuffer",
    "google",
    "latlong",
    "M1",
    "M2",
    "M3",
    "M4",
    "RD",
    "read_JSON_cat",
    "SPgfo",
)

## Unload packages to avoid namespace conflicts
warning("make_catalogue.R: Restart R now to clear loaded packages and avoid namespace conflicts.")
