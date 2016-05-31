########################################################################
#MMM - Oct 6, 2014
#This function appends information from a polygon into a point file.
#The user can select 1 or 2 fields they are interested in, and it will
#determine which polygon each point is located within.  Sometimes
#the polygon layer has several fields you might be interested in -
#this function can append up to 2 unique fields to a point file at once.
#
#this function expects:
#1) a csv containing the points
#2) a shapefile containing the polygons
#3) 1 or 2 fieldnames from the polygon file
#
#Example Usage:
#appending a single polygon field
# points.within.polys("surveypoints","Lobster_LFA_GRID_VDC","LFA",NULL)
#appending two polygon field
# points.within.polys("surveypoints","Lobster_LFA_GRID_VDC","LFA","GRID")
########################################################################

points.within.polys <- function(inputpointfile="",
                              inputpolyfile="",
                              level1Field="",
                              level2Field=""){
  if (is.null(level2Field)){
    hideLevel2=T
  }else{
    hideLevel2=F
  }

  library(maptools)
  bioLibrary( "bio.utilities", "bio.polygons")

  #find and load a point dataset
  #rename the fields to lat and lon (unnecessary)
  pointfile <- read.csv(file.path(find.bio.gis(inputpointfile)))
  pointfile<-rename.df(pointfile, "Lat.Start","lat")
  pointfile<-rename.df(pointfile, "Y","lat")
  pointfile<-rename.df(pointfile, "Long.Start","lon")
  pointfile<-rename.df(pointfile, "X","lon")
  #here we check the min value of the longitudes - if it's more than 0, we
  #multiply all longitudes by -1, since this data is probably in the wrong
  #hemisphere
  if (min(pointfile$lon, na.rm = T)>0){
    pointfile$lon<-pointfile$lon*-1
    print("wrong hemisphere - multiplied lon by -1")
  }
  #remove row if lat or lon is NA
  pointfile <- pointfile[!is.na(pointfile$lat) & !is.na(pointfile$lon),]

  #find and load a polygon dataset
  polygonfile <- readShapePoly(file.path(find.bio.gis(inputpolyfile)), proj4string=CRS("+proj=longlat"))


  level1<-U(polygonfile[,level1Field][[1]])
  this.poly<-NULL
  pointfile$loc<-rep("9999", nrow(pointfile))
  for (i in level1){
    if (hideLevel2){
      #if we're just doing a single level, make sure we loop through
      #all of the polygons for this "i"
        for (k in 1:length(polygonfile[(polygonfile[,level1Field][[1]] == i),]@polygons)){
          this.poly<-polygonfile[(polygonfile[,level1Field][[1]] == i),]@polygons[[k]]@Polygons[[1]]@coords
          pointstatus<-point.in.polygon(pointfile$lon,pointfile$lat,this.poly[,1],this.poly[,2])
          pointfile<-cbind(pointfile,pointstatus)
          pointfile$loc<-ifelse(((pointfile$pointstatus=="1") | (pointfile$pointstatus=="2")| (pointfile$pointstatus=="3")),paste(level1Field,"_",i,sep=""),ifelse(pointfile$loc!="9999",pointfile$loc,"none"))
          pointfile$pointstatus<-NULL
        }
      }else{
       level2<-(U(polygonfile[(polygonfile[,level1Field][[1]] == i),][,level2Field][[1]]))
       for (j in level2){
         this.poly<-polygonfile[(polygonfile[,level1Field][[1]] == i)&(polygonfile[,level2Field][[1]] == j),]@polygons[[1]]@Polygons[[1]]@coords
         pointstatus<-point.in.polygon(pointfile$lon,pointfile$lat,this.poly[,1],this.poly[,2])
         pointfile<-cbind(pointfile,pointstatus)
         pointfile$loc<-ifelse(((pointfile$pointstatus=="1") | (pointfile$pointstatus=="2")| (pointfile$pointstatus=="3")),paste(level1Field,"_",i,"_",level2Field,"_",j,sep=""),ifelse(pointfile$loc!="9999",pointfile$loc,"none"))
         pointfile$pointstatus<-NULL
       }
    }
  }
  #write.csv(pointfile,"pointswithinpolysResults.txt")
  #print(paste("Results saved to ", getwd(),"/pointswithinpolysResults.txt", sep=""))
  return(pointfile)
}
