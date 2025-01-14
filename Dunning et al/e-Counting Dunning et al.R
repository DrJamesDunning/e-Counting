# #thought experiment about echo-counting. See if fish detected multiple times affect density estimation
# 
# #install.packages(c("misc3d"))
# #install.packages(c("rgl"))
# #install.packages(c("abind"))
# #
# library(abind)
# library(rgl)
# library(misc3d)
# 
# 
# setwd("R:/CLSM/School of Biological Sciences/FisheriesGroup/PhD students/James Dunning/James Code/Echo Counting Simulation")
# source("bresenham.R")
# source("bresenHCirc.R")
# rotate <- function(x) t(apply(x, 2, rev))

### for development purposes
# Depth <- 10
# Distance <- 10
# Speed_kts <- 5
# PingInt <-  0.5
# BeamAngle_Deg <- 18
# MinTrack <- 3
#
# FishD_m <- 0.6
# FishLayer <- c(0,Depth) #No layer, fish are randomly distributed across full water column
# FishLayer <- c(3,7)     #Upper and lower layer of fish [m]
# Voxel <- 0.1         #length of edge of voxel in m
# Graph <- TRUE


#---- user input (comment out if using CallECS)
# Depth <- 10          #depth in m
# Distance <- 20      #length of transect in m
# Speed_kts <- 6       #speed of boat in knots
# PingInt <-  0.5      #ping interval in seconds
# BeamAngle_Deg <- 18  #beam angle in degrees
# MinTrack <- 2
# FishD_m <- 0.4         #fish density in fish/m^3
# FishLayer <- c(0,Depth) #No layer, fish are randomly distributed across full water column
# FishLayer <- c(3,7)
# 
# Voxel <- 0.1         #length of edge of voxel in m
# Graph <- FALSE

#---- start of program,  uncomment all of above if want to use it as a stand-alone program
DistBins  <- Voxel     #size of distance bins m
DepthBins <- Voxel     #size of depth bins in m
WidthBins <- Voxel     #size of width bins in m
FishD <- FishD_m*DistBins*DepthBins*WidthBins

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
  plot(BeamEdge, type = "p", asp=1, pch = 0) 
  lines(verts, xlab="Distance bins",ylab="Depth bins", col = 2, lwd = 2)
  points(-BeamEdge$x,BeamEdge$y, pch = 0)
  lines(-verts$x, verts$y, xlab="Distance bins",ylab="Depth bins", col = 2, lwd = 2)
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
for (q in 1:MaxDiametre) {
Fish[,,q] <- t(replicate(Ndepth, rbinom(Ndist, size = 1, prob=FishD)))  #generate fish in water
}

FishLayerTop <- min(FishLayer)
FishLayerBottom <- min (max(FishLayer),Depth) #Bottom of layer is equal to depth if it was set to be deeper than depth
DelFish <- rep(TRUE,Ndepth)

if (FishLayerTop <= Depth) {        #if layer of fish is within water column, keep fish in layer. If not, delete all fish
DelFish[(FishLayerTop/DepthBins+1):(FishLayerBottom/DepthBins)] <- FALSE #FALSE = keep fish, TRUE = delete fish, 
}

Fish[DelFish,,]  <- 0
FishID <- Fish
FishID[Fish==1] <- 1:sum(Fish) # array where each fish has an ID number

# create 3D ping ----

#2D ping at max width
BeamEdge$yMatrix <- -(BeamEdge$y-Ndepth)                          #convert y coordinates for matrix indexing
EdgeMatrixDx <- t(rbind(BeamEdge$yMatrix,c(BeamEdge$x+HalfPing)))  
Ping2D[EdgeMatrixDx] <- 0                                          #right edge of ping filled with 0 in matrix
EdgeMatrixSx <- t(rbind(BeamEdge$yMatrix,c(-BeamEdge$x+HalfPing)))
Ping2D[EdgeMatrixSx] <- 0         #left edge of ping filled with 0 in matrix

#3D ping 
for (i in 1:Ndepth) {
  Limits <- which(Ping2D[i,]==0)
  r <- (max(Limits)-min(Limits))/2 #radius, unit = distbins
  Circle <- bresenHCirc(r, Centre, fill = TRUE, plot = FALSE) #Bresenham's circle algorithm
  #image(Circle, asp=1)

  for (k in (Centre-r):(Centre+r)) {Ping3D[i,which(Circle[k,]==1),k] <- 0} #create slice of 3D ping sheet by sheet
}

VolPing <- sum(!is.na(Ping3D))*DistBins*DepthBins*WidthBins  #volume of ping in m^3
Width <- MaxDiametre*WidthBins                               #width of body of water in m
VolWater <- Depth*Distance*Width

R1 <- 0             #upper range of analysis layer [m]
R2 <- Depth         #lower range of analysis layer (= depth seabed or sqrt(Ndepth^2+HalfPing^2) ?) [m]
Rm_math <- ((R1^3+R2^3)/2)^(1/3)     #see SonarX manual, page 213, eq. 68 [m]
Area_math <- pi * Rm_math^2 * tan(phi/2) * tan(theta/2) #area of the SED beam in m^2 p.219 eq.87 [m^2]
Rm <- round(Rm_math/DepthBins)                                       #[voxels]

RadiusDepth_Math <- Depth*tan(theta/2)    #maximum radius arithmetic
RadiusEq <- sqrt(RadiusDepth_Math^2/3)    #radius of beam which defines cylinder with equal volume as cone
Req_math <- RadiusEq/tan(theta/2)         #Range (depth) at which equivalent radius is found (arithmetic)
Req <- round(Req_math/DepthBins)          #range (depht) in voxels

Area <- sum(Ping3D[Rm,,]==0, na.rm=T)*DistBins*WidthBins           #[m^2]
Area_Dpth <- sum(Ping3D[Ndepth,,]==0, na.rm=T)*DistBins*WidthBins  #[m^2]
Area_eq <- sum(Ping3D[Req,,]==0, na.rm=T)*DistBins*WidthBins        #[m^2]
CylAtRm <- Area*Depth                                              #[m^3]
CylatDpth <- Area_Dpth*Depth                                       #[m^3]
CylatEq <- Area_eq*Depth                                           #[m^3]
CylRatio <- CylAtRm/VolPing                       #value should be =1.889882=3/4^(1/3), but here may differ as digitalised from voxel sizes
Cones3 <- CylatDpth/VolPing                       #value should be =3 exact, but here may differ as digitalised from voxel sizes
CylRatio_eq <- CylatEq/VolPing                    #value should be =1 exact, but here may differ as digitalised from voxel sizes

# Another way to estimate CylAtRm (yields better estimate of CylRatio but less )
# VolHalfPing1 <- sum(Ping3D[1:Rm,,]==0, na.rm=T)*DistBins*WidthBins*DepthBins        #[m^3]
# VolHalfPing2 <- sum(Ping3D[Rm:Ndepth,,]==0, na.rm=T)*DistBins*WidthBins*DepthBins   #[m^3]
# VolHalfOutPing1 <- Area*Rm*DistBins - VolHalfPing1                                  #Upper cylinder vol - VolHalfPing1, which is x2 VolHalfPing1
# VolHalfOutPing2 <- VolHalfPing2-(Area*(Depth-Rm*DepthBins))                         #VolHalfPing2 - lower cylinder vol (cylinder is inside)
# CylAtRm <- VolPing*2 - VolHalfOutPing2                              #vol of full cylinder: 3/2 ping=upper cyl, 1/2 ping-VolHalfOutPing2=lower cyl
# CylRatio <- CylAtRm/VolPing  

# loop to simulate passing ping and detect fish in the water----
  for (i in 1:Npings) {  
  Ping <- Ping3D
  if (i>1) {
    NextPingSpace <- array(rep(NA, Ndepth*(NextPing*(i-1))*MaxDiametre),c(Ndepth, NextPing*(i-1), MaxDiametre))
    Ping <- Ping[,-((dim(Ping)[2]-(NextPing*(i-1)-1)):dim(Ping)[2]),]  #delete block at end of array -> moves ping along
    Ping <- abind(NextPingSpace,Ping,along=2)                          #add empty block to left of ping shape
    }

  AllPings[Ping==0] <- 1 #cumulative visual of all pings run so far 

  FishCounts[i] <- sum(Fish+Ping, na.rm=T)
  Densities[i] <- FishCounts[i]/VolPing
  Echogram[,i] <- rowSums(rowSums(Fish+Ping, dims = 2, na.rm=T))
  
  EchoTracksID[,i] <- rowSums(rowSums(FishID+Ping, dims = 2, na.rm=T)) #Echogram with Fish tracks
  }
EchoSingleT <- Echogram
EchoSingleT[EchoSingleT>1] <- 0 #Echogram with only single targets

AllPings_n_fish <- AllPings+Fish*2  #matrix: 0=empty water, 1=at least one ping, 2=undetected fish, 3= detected fish
AreaSurveyD <- sum(AllPings[Ndepth,,])*DistBins*WidthBins #Area surveyed, covered, at max depth (dipendent to pings)
AreaSurveyRm <- sum(AllPings[Rm,,])*DistBins*WidthBins    #Area surveyed, covered, at Rm (dipendent to pings)
#image(AllPings[Rm,,])

#Calculate Density contribution of each layer 1 voxel high, per ping 
SliceVolume <- rowSums(Ping3D+1,dims= 1, na.rm = TRUE)*DistBins*DepthBins*WidthBins #volume of each slice of ping [m^3]
MED_Yul_D <- colSums(Echogram*(1*DepthBins)/SliceVolume)  
#MED_Yul_D <- sum(rowSums(Echogram)*(1*DepthBins)/(SliceVolume*Npings)) this is equal to previous
SED_Yul_D <- colSums(EchoSingleT*(1*DepthBins)/SliceVolume) #Dens contr. of each detected fish (SED Yule method)


## fish track calculations #####################
#SailDistance = Distance sailed + 2*half beam (reach of first and last pings extend more than vertex) [m]
SailDistance <- (NextPing*(Npings-1))*DistBins
#SailDistance <- (NextPing*(Npings-1)+(MaxDiametre-1))*DistBins

# Find SED tracks (no minimum track length here)
EchoTracksID[Echogram>1] <- 0 #remove echotracks generated from multiple fish

k <-1
TrackIDs <- unique(EchoTracksID[EchoTracksID>0]) #unique values (exclude all zeros that are matrix background)


TrackCounts <- length(TrackIDs)
for (i in TrackIDs) { #renumber each track from 1 to k tracks
  EchoTracks[EchoTracksID==i] <- k
  k <- k+1
  }

LongTrackIDs <- as.numeric(names(which(table(EchoTracks)>=MinTrack))) #IDs of tracks with minimum #n of detections (can have gaps)
EchoTracksMTL <- EchoTracks
EchoTracksMTL[!(EchoTracks %in% LongTrackIDs)] <- NA


# calc Density contribution for each fish that has been insonified
if (all(AllPings_n_fish!=3)) {TrackSummary_MED<- matrix(0,nrow=1, ncol=4)  #if there are no fish in beam
} else {
TrackSummary_MED <- matrix(nrow=sum(AllPings_n_fish==3), ncol=4)

TrackSummary_MED[,1] <- FishID[AllPings_n_fish==3]
TrackSummary_MED[,2] <- which(FishID %in% TrackSummary_MED[,1]) #arr.ind=TRUE doesn't work when finding vectors
TrackSummary_MED[,3] <- arrayInd(TrackSummary_MED[,2],dim(FishID))[,1]*DepthBins #depth in m of each fish
TrackSummary_MED[,4] <-  1/(2*tan(theta/2)*TrackSummary_MED[,3])
}
colnames(TrackSummary_MED) <- c("FishID","Linear_indx","Fish_Depth","Dens_Contr[m^-1]")

# calc Density contribution for each SED fish track 
if (all(is.na(EchoTracks))) {TrackSummary_SED <- matrix(0,nrow=1, ncol=2)  #if there are no tracks detected
} else {
  TrackSummary_SED <- matrix(nrow=max(EchoTracks, na.rm = TRUE), ncol=2) 
  for (i in 1:max(EchoTracks, na.rm = TRUE)) {
    TrackSummary_SED[i,1] <- (which(EchoTracks==i, arr.ind=TRUE)[1])*DepthBins #depth in m of fish track
    TrackSummary_SED[i,2] <- 1/(2*tan(theta/2)*TrackSummary_SED[i,1])
  }
}
colnames(TrackSummary_SED) <- c("Track_Depth","Dens_Contr[m^-1]")

# calc Density contribution for each SED fish track with minimum track length
if (length(LongTrackIDs)==0) {TrackSummary_SEDLong <- matrix(0,nrow=1, ncol=4)  #if there are no tracks detected
} else {
  TrackSummary_SEDLong <- matrix(nrow=length(LongTrackIDs), ncol=4) 
  for (i in 1:length(LongTrackIDs)) {
    TrackSummary_SEDLong[i,1] <- (which(EchoTracks==LongTrackIDs[i], arr.ind=TRUE)[1])*DepthBins #depth in m of fish track
    TrackSummary_SEDLong[i,2] <- 1/(2*tan(theta/2)*TrackSummary_SEDLong[i,1])
    TrackSummary_SEDLong[i,3] <- LongTrackIDs[i]
    TrackSummary_SEDLong[i,4] <- table(EchoTracks)[LongTrackIDs[i]]
    
  }
}
colnames(TrackSummary_SEDLong) <- c("Track_Depth","Dens_Contr[m^-1]","TrackID","Track_length")



### Graphics -----
if (Graph) {
image3d(Ping3D, col = 1, alpha =1,  )    # 3D image of 1 ping
#image3d(Fish, col = c(8,2,4), alpha.power = 0.5)  # 3D image of fish
image3d(AllPings_n_fish, col = c("grey80","red2","green3"), alpha.power = 0.5)  # 3D summary of simulation. col = pings, undetected fish, detected fish


par(oma=c(0,0,0,0),mar=c(3,4,1,1))
yax <- seq(-Depth, 0, length = Ndepth)
xax <- seq(0, Distance, length = Ndist)

#Echogram
image(1:Npings, yax, rotate(Echogram),axes = FALSE, xlab="",ylab="") 
abline(v = c(0.5:Npings))
axis(1, at = seq(1, Npings, by = 1))
axis(2, at = seq(-Depth, 0, by = 10), labels = rev(seq(0,Depth, by =10)), las=1)
mtext(side = 1, text = "Ping no.", line = 2, cex=1.5)
mtext(side = 2, text = "Depth (m)", line = 3, cex=1.5)

#Fish Tracks Echogram no MTL
image(1:Npings, yax, rotate(EchoTracks),axes = FALSE, xlab="",ylab="", col = sample(rainbow(length(TrackIDs)))) 
abline(v = c(0.5:Npings))
axis(1, at = seq(1, Npings, by = 1))
axis(2, at = seq(-Depth, 0, by = 10), labels = rev(seq(0,Depth, by =10)), las=1)
mtext(side = 1, text = "Ping no.", line = 2, cex=1.5)
mtext(side = 2, text = "Depth (m)", line = 3, cex=1.5)

#Fish Tracks Echogram with MTL
image(1:Npings, yax, rotate(EchoTracksMTL),axes = FALSE, xlab="",ylab="", col = sample(rainbow(length(LongTrackIDs)))) 
abline(v = c(0.5:Npings))
axis(1, at = seq(1, Npings, by = 1))
axis(2, at = seq(-Depth, 0, by = 10), labels = rev(seq(0,Depth, by =10)), las=1)
mtext(side = 1, text = "Ping no.", line = 2, cex=1.5)
mtext(side = 2, text = "Depth (m)", line = 3, cex=1.5)

#dev.off()
}

## Density calculations ----

MED <- data.frame(TruVolD = c(NA,NA),       #this data frame stores calculations with multiple Echo detections (MED) 
                   TruNFish = c(NA,NA),
                   TruAreaD = c(NA,NA),
                   FSampled = c(NA,NA),
                   FSampled_pc =  c(NA,NA),
                   NTargets = c(NA,NA),
                   VolD = c(NA,NA),
                   NFish = c(NA,NA),
                   AreaD = c(NA,NA),
                   AreaDRm = c(NA,NA),
                   AreaDRmOff = c(NA,NA),
                   AreaDYul = c(NA,NA),
                   AreaDTrack = c(NA,NA),
                   AreaDTrackYul = c(NA,NA),
                   AreaDTrackX = c(NA,NA),
                   AreaDTrackXOff = c(NA,NA))

row.names(MED) <- c("Value","Error")   

MED$TruVolD[1] <- sum(Fish)/VolWater            #True fish volume density [m^3]
MED$TruNFish[1] <- sum(Fish)                    #True fish count
MED$TruAreaD[1] <- sum(Fish)/(Distance*Width)   #True fish area density [m^2]
MED$FSampled[1] <- sum(AllPings_n_fish==3)      #AllPings_n_fish <- AllPings+Fish*2, if ==3 then fish
MED$FSampled_pc[1] <- MED$FSampled[1]/MED$TruNFish[1]*100 #how many fish sampled in percentage
MED$NTargets[1] <- sum(FishCounts)            #how many fish counted (with resamples of same fish)
 
MED$VolD[1] <- sum(Echogram)/(VolPing*Npings)                        #Estimated f1sh density from Echogram [m^3]
MED$VolD[2] <- -(MED$TruVolD[1]-MED$VolD[1])/MED$TruVolD[1]*100   #% Error of estimate

MED$NFish[1] <- MED$VolD[1]*VolWater                                    #Estimated fish abundance from Echogram
MED$NFish[2] <- -(MED$TruNFish[1]-MED$NFish[1])/MED$TruNFish[1]*100   #% Error of estimate

MED$AreaD[1] <- mean(colSums(Echogram))*Depth/VolPing                  #Area density SonarX text above eq. 64
MED$AreaD[2] <- -(MED$TruAreaD[1]-MED$AreaD[1])/MED$TruAreaD[1]*100 #% Error of area density estimate

MED$AreaDRm[1] <- (sum(Echogram)*CylRatio)/(Area*Npings)                    #Area Density [fish m^-2] at Rm
MED$AreaDRm[2] <-  -(MED$TruAreaD[1]-MED$AreaDRm[1])/MED$TruAreaD[1]*100 #% Error of area density estimate

MED$AreaDRmOff[1] <- sum(Echogram)/(Area*Npings)                                 #Area Density [fish m^-2] at Rm without ratio factor
MED$AreaDRmOff[2] <- -(MED$TruAreaD[1]-MED$AreaDRmOff[1])/MED$TruAreaD[1]*100 #% Error of area density estimate

MED$AreaDYul[1] <- mean(MED_Yul_D)   #similar approach to (Tschersich, et al., 2015) but with echo-counts
MED$AreaDYul[2] <- -(MED$TruAreaD[1]-MED$AreaDYul[1])/MED$TruAreaD[1]*100 #% Error of area density estimate

MED$AreaDTrack[1] <- MED$FSampled[1]*Depth/(sum(AllPings)*DistBins*WidthBins*DepthBins)
MED$AreaDTrack[2] <- -(MED$TruAreaD[1]-MED$AreaDTrack[1])/MED$TruAreaD[1]*100 #% Error of area density estimate

MED$AreaDtrackYul[1] <- 1/SailDistance * sum(TrackSummary_MED[,"Dens_Contr[m^-1]"]) #Area Density [fi
MED$AreaDtrackYul[2] <- -(MED$TruAreaD[1]-MED$AreaDtrackYul[1])/MED$TruAreaD[1]*100 #% Error of estimate

MED$AreaDTrackX[1] <- MED$FSampled[1]*CylRatio/AreaSurveyRm                         #Area Density with tracks adapted from SonarX eq. 66
MED$AreaDTrackX[2] <- -(MED$TruAreaD[1]-MED$AreaDTrackX[1])/MED$TruAreaD[1]*100 #% Error of estimate

MED$AreaDTrackXOff[1] <- MED$FSampled[1]/AreaSurveyRm                                     #Area Density with tracks SonarX eq. 66
MED$AreaDTrackXOff[2] <- -(MED$TruAreaD[1]-MED$AreaDTrackXOff[1])/MED$TruAreaD[1]*100 #% Error of estimate


SED <- data.frame(TruVolD = c(NA,NA),           #this data frame stores calculations where SED is applied 
                  TruNFish = c(NA,NA),
                  TruAreaD = c(NA,NA),
                  FSampled = c(NA,NA),
                  NTargets = c(NA,NA),
                  FTracked = c(NA,NA),
                  SED_loss = c(NA,NA),
                  Track_loss = c(NA,NA),
                  VolD = c(NA,NA),
                  NFish = c(NA,NA),
                  AreaD = c(NA,NA),
                  AreaDRm = c(NA,NA),
                  AreaDRmOff = c(NA,NA),
                  AreaDYul = c(NA,NA),
                  AreaDTrack = c(NA,NA),
                  AreaDtrackYul = c(NA,NA),
                  AreaD_L_TrkYul = c(NA,NA),
                  AreaDTrackX = c(NA,NA),
                  AreaDTrackXOff = c(NA,NA))                  
row.names(SED) <- c("Value","Error")   


SED$TruVolD[1] <- sum(Fish)/VolWater            #True fish volume density [m^3]
SED$TruNFish[1] <- sum(Fish)                    #True fish count
SED$TruAreaD[1] <- sum(Fish)/(Distance*Width)   #True fish area density [m^2]

SED$FSampled[1] <- sum(AllPings_n_fish==3)      #AllPings_n_fish <- AllPings+Fish*2, if ==3 then fish
SED$NTargets[1] <- sum(EchoSingleT)            #how many fish counted (with resamples of same fish)
SED$FTracked[1] <- TrackCounts
SED$SED_loss[1] <- (MED$NTargets[1]-SED$NTargets[1])/MED$NTargets[1]*100 #percentage of targets lost due to SED
SED$Track_loss[1] <- (SED$FSampled[1]-SED$FTracked[1])/SED$FSampled[1]*100 #percentage of fish lost due to SED


SED$VolD[1] <- sum(EchoSingleT)/(VolPing*Npings)                #Estimated f1sh density from echogram
SED$VolD[2] <- -(SED$TruVolD[1]-SED$VolD[1])/SED$TruVolD[1]*100 #% Error of estimate

SED$NFish[1] <- SED$VolD[1]*VolWater                                 #Estimated fish abundance from echogram[m^3]
SED$NFish[2] <- -(SED$TruNFish[1]-SED$NFish[1])/SED$TruNFish[1]*100  #% Error of estimate

SED$AreaD[1] <-  mean(colSums(EchoSingleT))*Depth/VolPing           #Area density SonarX text above eq. 64
SED$AreaD[2] <- -(SED$TruAreaD[1]-SED$AreaD[1])/SED$TruAreaD[1]*100 #% Error of estimate

SED$AreaDRm[1] <- (sum(EchoSingleT)*CylRatio)/(Area*Npings)             #Area Density [fish m^-2] at Rm
SED$AreaDRm[2] <- -(SED$TruAreaD[1]-SED$AreaDRm[1])/SED$TruAreaD[1]*100 #% Error of estimate

SED$AreaDRmOff[1] <- sum(EchoSingleT)/(Area*Npings)                           #Area Density [fish m^-2] at Rm without ratio factor
SED$AreaDRmOff[2] <- -(SED$TruAreaD[1]-SED$AreaDRmOff[1])/SED$TruAreaD[1]*100 #% Error of estimate

SED$AreaDYul[1] <- mean(SED_Yul_D)   #similar approach to (Tschersich, et al., 2015) but with echo-counts
SED$AreaDYul[2] <- -(SED$TruAreaD[1]-SED$AreaDYul[1])/SED$TruAreaD[1]*100 #% Error of area density estimate

SED$AreaDTrack[1] <- length(LongTrackIDs)*Depth/(sum(AllPings)*DistBins*WidthBins*DepthBins) #use TrackCounts to disregard MTL
SED$AreaDTrack[2] <- -(SED$TruAreaD[1]-SED$AreaDTrack[1])/SED$TruAreaD[1]*100 #% Error of area density estimate

SED$AreaDtrackYul[1] <- 1/SailDistance * sum(TrackSummary_SED[,"Dens_Contr[m^-1]"]) #Area Density [fish m^-2] (Tschersich, et al., 2015)
SED$AreaDtrackYul[2] <- -(SED$TruAreaD[1]-SED$AreaDtrackYul[1])/SED$TruAreaD[1]*100 #% Error of estimate

SED$AreaD_L_TrkYul[1] <- 1/SailDistance * sum(TrackSummary_SEDLong[,"Dens_Contr[m^-1]"]) #Area Density [fish m^-2] (Tschersich etal 2015)
SED$AreaD_L_TrkYul[2] <- -(SED$TruAreaD[1]-SED$AreaD_L_TrkYul[1])/SED$TruAreaD[1]*100 #% Error of estimate

SED$AreaDTrackX[1] <- TrackCounts*CylRatio/AreaSurveyRm                         #Area Density with tracks adapted from SonarX eq. 66
SED$AreaDTrackX[2] <- -(SED$TruAreaD[1]-SED$AreaDTrackX[1])/SED$TruAreaD[1]*100 #% Error of estimate

SED$AreaDTrackXOff[1] <- TrackCounts/AreaSurveyRm                                     #Area Density with tracks SonarX eq. 66
SED$AreaDTrackXOff[2] <- -(SED$TruAreaD[1]-SED$AreaDTrackXOff[1])/SED$TruAreaD[1]*100 #% Error of estimate

