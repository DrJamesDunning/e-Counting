#this scipt is called by ECountingSim.V8 to generate an environment matrix that mimics 
#a given fish density profile from a .csv file


#DensData <- read.csv("Ch5DensProfileEdited.csv") 
Depth <- max(abs(DensData$FishDepth))
DistBins  <- Voxel     #size of distance bins m
DepthBins <- Voxel     #size of depth bins in m
WidthBins <- Voxel     #size of width bins in m
FishLayerTop <- min(FishLayer)
FishLayerBottom <- min (max(FishLayer),Depth) #Bottom of layer is equal to depth if it was set to be deeper than depth

DensProfile <- data.frame(FishDepth = 1:(Depth/DepthBins)) #select the columns in .csv to use
DensProfile$FishD_m <- NA
DensProfile$FishD_m[-DensData$FishDepth/DepthBins] <- DensData$SsmoothDens
if (FishLayerTop > 0) {DensProfile$FishD_m[1:(FishLayerTop/DepthBins)] <- 0}
if (FishLayerBottom < Depth) {DensProfile$FishD_m[(FishLayerBottom/DepthBins):(Depth/DepthBins)] <- 0}
DensProfile$FishD_m <- spline(1:(Depth/DepthBins), DensProfile$FishD_m,xout=1:(Depth/DepthBins))$y #interpolate
DensProfile$FishD_m[DensProfile$FishD_m<0] <- 0 #replace negative fish densities from interpolation with 0
DensProfile$FishD_m <- DensProfile$FishD_m*Multiplier



DensProfile$FishD <- DensProfile$FishD_m*DistBins*DepthBins*WidthBins

Ndepth <- Depth/DepthBins                               #number of depth bins
Ndist <-  Distance/DistBins                             #number of distance bins
NextPing_m <- (Speed_kts*1.852*1000/3600) * PingInt     #speed of boat in m/s * Pint interval, distance to next ping in m
NextPing <- round(NextPing_m/DistBins)                  #how many distance bins to the next ping

BeamAngle <- BeamAngle_Deg * pi/180      #beam angle in radians
phi <- BeamAngle       #athwart opening angle influenced by the Max gain compensation parameter in SED, eq.79
theta <- BeamAngle     #along ship opening angle influenced by the Max gain compensation parameter in SED, eq.79
m <- tan((90-(BeamAngle_Deg/2))*pi/180) # slope of left side beam limit
xIntercept <- Ndepth/m                  # intercept on x axis of right side beam (left side would be negative)
xIntercept <- xIntercept*DepthBins/DistBins #convert from depth bins units to distance bins unis

verts <- list(x = c(0, xIntercept), y = c(Ndepth, 0)) #start and end coordinates of right side beam

BeamEdge <- bresenham(verts)                    #Generate integer x,y points between vertices by Bresenham's algorithm
BeamEdge$x <- BeamEdge$x[2:length(BeamEdge$x)]  #remove first x coordinate (would end up being 0,0 in matrix)
BeamEdge$y <- BeamEdge$y[2:length(BeamEdge$y)]  #remove first y coordinate (would end up being 0,0 in matrix)
HalfPing <- max(BeamEdge$x)+1                        #x axis coordinate of first ping vertex
MaxDiametre <- HalfPing*2-1                          #diametre (in pixels) of beam at max depth
Centre <- ((MaxDiametre-1)/2)+1                      #centre of all concentric beam cross-sections in MaxDiametre matrix
Npings <- floor((Ndist-(HalfPing*2-1))/NextPing)+1   #How many pings fit in the transect
if (Graph) {
  plot(verts, type = "l", xlab="Distance bins",ylab="Depth bins")
  points(BeamEdge) 
}

#initialise variables
Fish <- array(rep(NA, Ndepth*Ndist*MaxDiametre),c(Ndepth, Ndist, MaxDiametre)) #rows*columns*sheets
Ping3D <- array(rep(NA, Ndepth*Ndist*MaxDiametre),c(Ndepth, Ndist, MaxDiametre))
Ping2D <- matrix(nrow=Ndepth, ncol=Ndist)
AllPings <- array(rep(0, Ndepth*Ndist*MaxDiametre),c(Ndepth, Ndist, MaxDiametre))
Echogram <- matrix(nrow=Ndepth, ncol=Npings)
EchoTracksID <- matrix(nrow=Ndepth, ncol=Npings)  #this will store echogram with fish tracks and fish IDs
EchoTracks <- matrix(nrow=Ndepth, ncol=Npings)    #this will store echogram with fish tracks, counted in order
FishCounts <- vector()
Densities <- vector()

# create fish in water with layer----
for (q in 1:Ndepth) {
  Fish[q,,] <- t(replicate(MaxDiametre, rbinom(Ndist, size = 1, prob=DensProfile$FishD[q])))  #generate fish in water
  }

FishLayerTop <- min(FishLayer)
FishLayerBottom <- min (max(FishLayer),Depth) #Bottom of layer is equal to depth if it was set to be deeper than depth
DelFish <- rep(TRUE,Ndepth)

if (FishLayerTop <= Depth) {        #if layer of fish is within water column, keep fish in layer. If not, delete all fish
  DelFish[(FishLayerTop/DepthBins):(FishLayerBottom/DepthBins)] <- FALSE #FALSE = keep fish, TRUE = delete fish, 
}

Fish[DelFish,,]  <- 0
FishID <- Fish
FishID[Fish==1] <- 1:sum(Fish) # array where each fish has an ID number


#image3d(AllPings_n_fish, col = c(8,2,4), alpha.power = 0.5)  # 3D summary of simulation. col = pings, undetected fish, detected fish
SimFishD <- apply(Fish, 1, sum)/((Ndist*MaxDiametre*1)*(DistBins*DepthBins*WidthBins))
plot(SimFishD,-(1:Ndepth)*DepthBins, xlab="",ylab="Fish depth (m)", yaxt='n', type ="l", lwd = 2, col = 2)
lines(DensProfile$FishD_m,-DensProfile$FishDepth*DepthBins, type ="l", lwd = 2, col =1)
axis(side=2, las=2, at = seq(-Ndepth*DepthBins,0, by=2), lab=(rev(seq(0,Ndepth*DepthBins, by=2))))
mtext(side=1, line=2.5, expression(Fish~density~(fish/m^"3")) )
legend(8,-70,legend = c("Data","Simulation"), col=c(1,2),lty = 1, seg.len = 0.5, lwd = 5)

