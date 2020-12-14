read_JSON_cat <- function(Tpath="data/EC.json",M1=M1,M2=M2,M3=M3,M4=M4,gfoPBSbuffer=gfoPBSbuffer) {
  EQ <- fromJSON(Tpath)

  tpoints <- cbind(as.numeric(EQ$events$lon), as.numeric(EQ$events$lat))
  tpoints2 <- SpatialPointsDataFrame(tpoints, proj4string = CRS(latlong), data=data.frame(ID=1:nrow(tpoints)))
  tpoints3 <- spTransform(tpoints2, CRS(RD))

  EQ$events$RD_X <- tpoints3@coords[,1]
  EQ$events$RD_Y <- tpoints3@coords[,2]
  rm(tpoints,tpoints2)

  tpoints <- cbind(as.numeric(EQ$events$lon), as.numeric(EQ$events$lat))
  tpoints2 <- SpatialPointsDataFrame(tpoints, proj4string = CRS(latlong), data=data.frame(ID=1:nrow(tpoints)))
  tpoints3 <- spTransform(tpoints2, CRS(RD))

  EQ$events$RD_X <- tpoints3@coords[,1]
  EQ$events$RD_Y <- tpoints3@coords[,2]
  rm(tpoints,tpoints2, tpoints3)

  EQdat <- EQ$events

  EQ$events$Year <- as.numeric(substr(EQ$events$date,1,4))
  EQdat$Date <- strptime(EQdat$date,format="%Y-%m-%d", tz="CET")
  EQdat$Year <- format(EQdat$Date,"%Y")
  EQdat$Month <- format(EQdat$Date,"%m")
  EQdat$Hour <- as.numeric(substr(EQdat$time,1,2))
  EQdat$Minute <- as.numeric(substr(EQdat$time,4,5))
  EQdat$Second <- as.numeric(substr(EQdat$time,7,8))
  EQdat$Hour2 <- rep(0,nrow(EQdat))
  EQdat$Hour2[EQdat$Hour%in%c(0:22) & EQdat$Minute>30] <- EQdat$Hour[EQdat$Hour%in%c(0:22) & EQdat$Minute>30]+1
  EQdat$Hour2[EQdat$Hour%in%c(23) & EQdat$Minute>30] <- 0
  EQdat$Hour2[EQdat$Minute<=30] <- EQdat$Hour[EQdat$Minute<=30]
  EQdat$WD <- format(EQdat$Date,format="%u")
  Tdate <- ymd(EQdat$date)
  EQdat$decimal <- decimal_date(Tdate)
  EQdat$julian <- julian(EQdat$Date)
  EQdat$julian2 <- as.numeric(EQdat$julian)*24 + EQdat$Hour + EQdat$Minute/(60) + EQdat$Second/(60*60)
  EQdat$dt <- EQdat$julian2/24
  EQdat$mag <- as.numeric(EQdat$mag)

  EQdat$Quarter <- rep(NA,nrow(EQdat))
  EQdat$Quarter[EQdat$Month%in%c("01","02","03")] <- 1
  EQdat$Quarter[EQdat$Month%in%c("04","05","06")] <- 2
  EQdat$Quarter[EQdat$Month%in%c("07","08","09")] <- 3
  EQdat$Quarter[EQdat$Month%in%c("10","11","12")] <- 4

  EQdat$INgfo <- point.in.polygon(EQdat$RD_X, EQdat$RD_Y, gfoPBSbuffer$X, gfoPBSbuffer$Y, mode.checked=FALSE)

  EQdat$South <- point.in.polygon(EQdat$RD_X,EQdat$RD_Y, M1[,1],M1[,2])
  EQdat$SouthWest <- point.in.polygon(EQdat$RD_X,EQdat$RD_Y, M2[,1],M2[,2])
  EQdat$NorthEast <- point.in.polygon(EQdat$RD_X,EQdat$RD_Y, M3[,1],M3[,2])
  EQdat$NorthWest <- point.in.polygon(EQdat$RD_X,EQdat$RD_Y, M4[,1],M4[,2])

  EQdat$Location <- rep(0,nrow(EQdat)); EQdat$Location[EQdat$South==1] <- 1; EQdat$Location[EQdat$SouthWest==1] <- 2; EQdat$Location[EQdat$NorthWest==1] <- 3; EQdat$Location[EQdat$NorthEast==1] <- 4;

  EQIn <- EQdat[EQdat$INgfo==1,c("date","time","RD_X","RD_Y","mag","place","Location","type","dt","lat","lon","julian","julian2","decimal","Hour","Month","Quarter","Year")]
  dimnames(EQIn)[[1]] <- as.character(1:nrow(EQIn))
  return(EQIn)
}
