
#Reading the data
env<-read.csv("enviromentalstat.csv")
View(env)
summary(env)
# loading some packages
library(ggplot2)
library(dplyr);library(viridis)
library(ggExtra)
library(tidyverse) 
library(sf)
library(sp)
library(spdep)
library(rgdal)
library(rgeos)
library(tmap)
library(tmaptools)
library(spgwr)
library(grid)
library(gridExtra)
library(spatial)
library(plot3D)
library(plot3Drgl)

#Part1: Exploring the data
#Box Plot
ggplot(data=env,aes(y=AVR,x=""))+
  geom_boxplot(fill = "lightblue",alpha=0.3,outlier.colour ="Purple",colour="blue")+ 
  ggtitle("Temperature Boxplot") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Histogram
ggplot(data=env, aes(x =AVR)) +
  geom_histogram(binwidth = 1, fill = "lightgreen") +
  ggtitle("Average Temperature (Second Quarter)") +
  xlab("Average Temperature") +
  ylab("Frequency")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#Or for the histogram
hist(env$AVR,main= "Histogram", ylim = c(0,200))


#bubble plot(1)
radius <- sqrt(env$AVR/ pi )
symbols(env$LONG, env$LAT, circles=radius)
symbols(env$LONG, env$LAT, circles=radius, inches=0.25, fg="white", bg="red")

plot<-ggplot()+
  geom_point(data=env, aes(LONG,LAT),col="blue", shape=19, size=2.5)+
  theme_bw()+
  xlab("Longitude")+
  ylab("Latitude")
plot
#bubble plot(2)
# Create a scatter plot with points colored by "average"
ggplot(env, aes(x = LONG, y = LAT, color = AVR)) +
  # Add points with size proportional to a third variable
  geom_point(aes(size = AVR)) +
  # Add a legend for the color scale
  scale_color_gradient(low = "blue", high = "red")

#Scatter Plots
install.packages("plotly")   # install the plotly package
library(plotly)              # load the plotly package

scatter3D(env$LONG,env$LAT,env$AVR, zcol=env$AVR)
scatter3D(env$LONG,env$LAT,env$AVR, zcol=env$AVR,pty="g",ticktype="detailed")
library("plot3Drgl")
plotrgl()

#GWRM as a method of exploration
library(spgwr)
library(ggplot2)

model1 <- lm(env$AVR ~ env$LONG + env$LAT)
summary(model1)
plot(model1)

par(mfrow=c(2,2))
plot(model1)

resids<-residuals(model1)
colours <- c("dark blue", "blue", "red", "dark red")
map.resids <- SpatialPointsDataFrame(data=data.frame(resids), coords=cbind(env$LONG,env$LAT)) 
spplot(map.resids, cuts=quantile(resids), col.regions=colours, cex=1) 


#calculate kernel bandwidth
GWRbandwidth <- gwr.sel(env$AVR ~ env$LAT+env$LONG, data=env, coords=cbind(x,y),adapt=T) 

gwr.model = gwr(env$AVR ~ env$LAT+env$LONG, data=env, coords=cbind(x,y), adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

gwr.model

results<-as.data.frame(gwr.model$SDF)
head(results)

env$coefLONG<-results$LONG
env$coefLAT<-results$LAT
########ERROR
gwr.point1<-ggplot(env, aes(x=x,y=y))+geom_point(aes(colour=data$coefLONG))+scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, space = "rgb", na.value = "grey50", guide = "colourbar", guide_legend(title="Coefs"))
gwr.point3+geom_path(env,aes(Long, Lat, group=id), colour="grey")+coord_equal()


gwr.point3<-ggplot(env, aes(x=x,y=y))+geom_point(aes(colour=data$coefLAT))+scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, space = "rgb", na.value = "grey50", guide = "colourbar", guide_legend(title="Coefs"))


#Part2: Calculating Moran's I

#Global Moran
library(ape)
data.dists<-as.matrix(dist(cbind(env$LAT,env$LONG)))
data.dist.inv<-1/data.dists
diag(data.dist.inv)<-0
Moran.I(env$AVR,data.dist.inv)

#Localized Moran
library(lctools)

l.moran2<-l.moransI(data.dists,20,env$AVR, WType='Bi-square', scatter.plot = TRUE, family = "adaptive")
l.moran2
library(spdep)
library(spatialEco)
library(sp)

morans.plot(env$AVR, coords =coordinates(env))
crs    <- CRS("+init=epsg:28992") 

coords <- env[ , c("LONG", "LAT")]

spdf <- SpatialPointsDataFrame(coords      = coords,
                               data        = env, 
                               proj4string = crs)
spplot(spdf)
library(FNN)
library(ncf)


knn<-knearneigh(spdf, k=20, longlat = NULL)
knn2nb<-knn2nb(knn)
local <- localmoran(x = env$AVR, listw = nb2listw(knn2nb))

mp <- moran.plot(env$AVR, nb2listw(knn2nb))

#lisa plot 
library(ncf)
lisa <- lisa(env$LONG, env$LAT, env$AVR, neigh = 20, latlon=TRUE)
lisa
plot(lisa)


#Part3: Trend Surface

#fitting a trend surface model by least squares
library(spatial)
x<-env$LONG
y<-env$LAT
z<-env$AVR 

fit.sfc1 <- surf.ls(1,x,y,z)
summary(fit.sfc1)
trsurf1 <- trmat(fit.sfc1, min(x), max(x), min(y), max(y), 50)
scatter3D(x,y,z,surf = list(x = trsurf1$x, y = trsurf1$y, z = trsurf1$z,
                            NAcol = "grey", shade = 0.1))
plotrgl()

fit.sfc2 <- surf.ls(2,x,y,z)
summary(fit.sfc2)
fit.sfc2$beta
trsurf2 <- trmat(fit.sfc2, min(x), max(x), min(y), max(y), 50)
scatter3D(x,y,z,surf = list(x = trsurf2$x, y = trsurf2$y, z = trsurf2$z,
                            NAcol = "grey", shade = 0.1))
plotrgl()


fit.sfc3 <- surf.ls(3,x,y,z)
summary(fit.sfc3)
fit.sfc3$beta
trsurf3 <- trmat(fit.sfc3, min(x), max(x), min(y), max(y), 50)
scatter3D(x,y,z,surf = list(x = trsurf3$x, y = trsurf3$y, z = trsurf3$z,
                            NAcol = "grey", shade = 0.1))
plotrgl()

fit.sfc4 <- surf.ls(4,x,y,z)
summary(fit.sfc4)
trsurf4 <- trmat(fit.sfc4, min(x), max(x), min(y), max(y), 50)
scatter3D(x,y,z,surf = list(x = trsurf4$x, y = trsurf4$y, z = trsurf4$z,
                            NAcol = "grey", shade = 0.1))
plotrgl()

contour(trsurf4)
contour(trsurf3)

#Part4: IDW 
#Get the optimal value of p
########
library(spatstat)

obs_window <- owin(xrange = c(-124.2,-114), yrange =c(32.55,41.98) )

ppp_av<-ppp(env$LONG,env$LAT,
            marks=env$AVR,window=obs_window)

idw_av <- idw(ppp_av, power=0.05, at="pixels")

plot(idw_av,
     col=heat.colors(20), 
     main="Interpolated spatial variation in tempreture based on IDW\n (Power = 0.05)") 

idw_points <- idw(ppp_av, power = 0.05, at = "points")
library(Metrics)

Metrics::mse(ppp_av$marks, idw_points)

powers <- seq(0.001, 5, 0.01)
mse_result <- NULL
for(power in powers){
  CV_idw <- idw(ppp_av, power=power, at="points")
  mse_result <- c(mse_result,
                  Metrics::mse(ppp_av$marks,CV_idw))
}
optimal_power <- powers[which.min(mse_result)]
optimal_power

plot(powers, mse_result)

#####Optimal p=3.881 
#To get MSE for p=3.881 and plot it

idw_av1 <- idw(ppp_av, power=3.881, at="pixels")


plot(idw_av1,
     col=heat.colors(20), 
     main="Interpolated spatial variation in tempreture based on IDW\n (Power = 3.881)") 

idw_points <- idw(ppp_av, power = 3.881, at = "points")
library(Metrics)

Metrics::mse(ppp_av$marks, idw_points)

#IDW with p=3.881
library(phylin)
library(sp) 

coords <- env[ , c("LONG", "LAT")]

data1  <- env[ , c("AVR")]

#This Line of the code used for projection
crs    <- CRS("+init=epsg:28992") # proj4string of coords

# make the SpatialPointsDataFrame object
spdf <- SpatialPointsDataFrame(coords      = coords,
                               data        = env, 
                               proj4string = crs)


#Making Interpolated maps

Long <- seq(from=-124,to=-114.16,by=0.1)
Lat <- seq(from=32.55,to=41.98,by=0.1)
grid<-cbind(rep(Long,length(Lat)), rep(Lat,each=length(Long)))
library(phylin)

idw<- phylin::idw(data1, coords, grid)


grid.image(idw, grid, p=3.881, main='IDW interpolation', xlab='Longitude', 
           ylab='Latitude', sclab="Genetic distance to sample s2")
points(coords, cex=data1/50)


mycolors <- colorRampPalette(c("#6a0dad", "#9e36b2", "#d163be", "#f89ece", "#b1dafb"))
grid.image(idw, grid, p=3.881, main='IDW interpolation', xlab='Longitude', 
           ylab='Latitude', sclab="Genetic distance to sample s2", colFUN = mycolors)
points(coords, cex=data1/50)

#Part5: Kriging
library(gstat)

evgm <- variogram( env[, c("AVR")]~1,spdf)
plot(evgm)

#Spherical Variogram and its MSE
fvgm1 <- fit.variogram(evgm,vgm("Sph"))
plot(evgm,model=fvgm1)
attr(fvgm1, "SSErr")

#Gaussian Variogram and its MSE
fvgm2 <- fit.variogram(evgm,vgm("Gau"))
plot(evgm,model=fvgm2)
attr(fvgm2, "SSErr")

#Exponential Variogram and its MSE
fvgm3 <- fit.variogram(evgm,vgm("Exp"))
plot(evgm,model=fvgm3)
attr(fvgm3, "SSErr")

##we'll choose exponential

s.grid <- spsample(spdf, type = "regular", n = 6000)
krig.est <- krige(env[, c("AVR")] ~ 1, spdf, newdata = s.grid, model = fvgm3)
spplot(krig.est)
spplot(krig.est["var1.pred"])
spplot(krig.est["var1.var"])
summary(krig.est)





