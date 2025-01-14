setwd("R:/CLSM/School of Biological Sciences/FisheriesGroup/PhD students/James Dunning/James Code/Echo Counting Simulation")
library(beepr)

source("bresenHCirc.R")
rotate <- function(x) t(apply(x, 2, rev))

progress <- function (x, max = 100) {
  percent <- x / max * 100
  cat(sprintf('\r[%-50s] %d%%',paste(rep('=', percent / 2), collapse = ''), floor(percent)))
  if (x == max) {cat('\n')}
}

start_time <- Sys.time()

# ---- user input
Depth <- seq(14.5,15, by=0.1)  # doesn't accept zero [m]
#Depth <- 1
BeamAngle_Deg <- 18         #deg
Speed_kts <- 6              #[kts]
PingInt <- 0.5               #[s]
MTL <- 3                    # Minimum track length [pings]
pixel <- 0.001               #[m]
SaveMethods <- F           #needs only one depth value
SaveResultsRainbow <- F           #needs multiple depths
SaveResultsLine    <- F

DirFig <- paste("R:/CLSM/School of Biological Sciences/FisheriesGroup/PhD students/James Dunning/",
                "Thesis/Ch.7 - Annex/Figures/R export temp/", sep="")

# ---- start of program
BeamAngle <- BeamAngle_Deg * pi/180      #beam angle in radians
decimals <- log10(1/pixel)
NextPing_m <- (Speed_kts*1.852*1000/3600) * PingInt #speed of boat in m/s * Pint interval, distance to next ping in m
NextPing_Px <- round(round(NextPing_m,decimals)/pixel)
DistEdge_m <- NextPing_m * (MTL-1) 


r_max <- max(Depth)*tan(BeamAngle/2)                  #to allocate sizes, similar to r
r_Pix_max <- floor((round(r_max,decimals)/pixel)/2)*2
NY_Pix_max <- r_Pix_max*2+1

Probs <- list(Trues = matrix(nrow=length(Depth), ncol=NY_Pix_max),
              Falses = matrix(nrow=length(Depth), ncol=NY_Pix_max),
              ProbSum = matrix(nrow=length(Depth), ncol=NY_Pix_max),
              AxisDist = matrix(nrow=length(Depth), ncol=NY_Pix_max),
              AngleDist = matrix(nrow=length(Depth), ncol=NY_Pix_max)              )  

for (k in 1:length(Depth)) {
  
  ProgText <- paste("  Performing depth =",Depth[k],"; setting",k,"of",length(Depth),"      ")
  cat(ProgText,'\r')
  progress(k/length(Depth)*100)
  
  r <- Depth[k]*tan(BeamAngle/2)                 #radius of cross-section at depth [m]
  r_Pix <- floor((round(r,decimals)/pixel)/2)*2  #radius of cross-section at depth [pixel]
  
  Centre_Pix <- r_Pix+1
  Diametre_Pix <- r_Pix*2+1
  AnalysPing <- ceiling(Diametre_Pix/NextPing_Px)   #the first ping with a full pattern inside
  AllPings <- AnalysPing*2-1                       #how many pings run by simulation (used not to have -1)
  NY_Pix <- Diametre_Pix                           
  NX_Pix <- (r_Pix+NextPing_Px*(AllPings-1)+r_Pix+1)
  
  Circle <- bresenHCirc(r_Pix, r_Pix+1, fill = TRUE, plot = FALSE)
  Circle[is.na(Circle)] <- 0
  
  
  Sheet <- matrix(nrow=NY_Pix, ncol=NX_Pix,0) 
  Sheet[,1:NY_Pix] <- Sheet[,1:NY_Pix] + Circle
  #image(rotate(Sheet)) # within for loop, first beam (step repeated at different depths)
  
  # Creates footprint matrix
  for (i in 2:AllPings) {
    window <- (NextPing_Px*(i-1)+1):(NextPing_Px*(i-1)+Diametre_Pix)
    Sheet[,window] <- Sheet[,window] + Circle
    #image(rotate(Sheet))                     #addition of all overlapping beams (one by one)
    #readline(prompt="continue next ping")
  }
  #image(rotate(Sheet))  #complete footprint (with degrees) at a given depth
  
  # bitmap with repeating pattern of detect/not-detect
  ProcWindow <- Sheet[,(NextPing_Px*(AnalysPing-1)+1):(NextPing_Px*(AnalysPing-1)+NY_Pix)]
  if(SaveMethods) {PatternProc <- ProcWindow}
  #image(rotate(ProcWindow), asp =1)
  #contour(rotate(ProcWindow), add = TRUE)
  
  ProcWindow[ProcWindow < MTL] <- 0
  ProcWindow[ProcWindow >= MTL] <- 1
  if(SaveMethods) {BitmapProc <- ProcWindow}
  
  #image(rotate(ProcWindow), asp =1)
  #ProcWindow[Circle==0] <- -1
  #Circonference <- bresenHCirc(r_Pix, r_Pix+1, fill = FALSE, plot = FALSE)
  #Circonference[is.na(Circonference)] <- 0
  #image(rotate(ProcWindow+Circonference*5), asp =1)
  
  #cut bitmap into a repeating pattern
  PatSize <- colSums(ProcWindow) #vertical (plot view) size of bitmap at each x-axis
  MaxPat <- max(PatSize)
  MinPat <- min(PatSize[PatSize>0])
  MaxPoints <- which(PatSize == MaxPat)
  MinPoints <- which(PatSize == MinPat)
  
  CutPoints <- MinPoints  #set CutPoints as troughs, if condition below is true set as peaks
  if (min(PatSize)==0 | MaxPoints[1] < MinPoints[1]) {CutPoints <- MaxPoints}
  
  #first and last CutPoints should be correct, but need to adjust for consecutive pixels not forming a point
  Repeats <- table(cumsum(c(1, diff(CutPoints) != 1))) #consecutive repeats of max or min cut points
  if (length(Repeats)==1) {PixAdj <- 0}
  if (length(Repeats)==2) {PixAdj <- Repeats[1]}
  if (length(Repeats)==3) {PixAdj <- as.numeric((Repeats[1]+Repeats[length(Repeats)])-Repeats[2])}
  if (length(Repeats)>3) {CutPoints <- CutPoints[(Repeats[1]+1):(length(CutPoints)-sum(Repeats[c(1,2)])+1)]
                          PixAdj <- 1}
  AllPoints <- all(Repeats == 1) #if all cut points are made of just one pixel then true
  
  #largest repeating pattern available 
  #if (EvenOddPat%%2 == 0) {ProcWindow <- ProcWindow[,CutPoints[EvenOddPat/2]:(CutPoints[length(CutPoints)]-EvenOddPat/2)]}
  if (AllPoints) {ProcWindow <- ProcWindow[,CutPoints[1]:(CutPoints[length(CutPoints)]-1)]
  } else {ProcWindow <- ProcWindow[,(CutPoints[1]):(CutPoints[length(CutPoints)]-PixAdj)]}
  #image(rotate(ProcWindow))
  
  
  Probs$Tot[k] <- sum(ProcWindow)/length(ProcWindow) #non-overlapping bias + MTL bias
  for (i in 1:NY_Pix) {
    Probs$Trues[k,i] <- sum(ProcWindow[i,] == 1)/length(ProcWindow[i,])
    Probs$Falses[k,i] <- sum(ProcWindow[i,] == 0)/length(ProcWindow[i,])
    Probs$ProbSum[k,i] <- Probs$Trues[k,i] + Probs$Falses[k,i]
    if (i<=r_Pix) {Probs$AxisDist[k,i] <- r_Pix-(i-1)}
    if (i==r_Pix+1) {Probs$AxisDist[k,i] <- 0}
    if (i>r_Pix+1) {Probs$AxisDist[k,i] <- (1+r_Pix)-i}
    
    Probs$AngleDist[k,i] <- atan(Probs$AxisDist[k,i]*pixel/Depth[k]) * 180/pi
  }
}

par(mar=c(5, 4, 4, 8), xpd=TRUE)
matplot(t(Probs$AxisDist)*pixel,t(Probs$Trues), xlab="Dist from axis [m]" , ylab="", main = "",
        type = "l", lty = 1, col= 1:length(Depth),  pch = 16)
legend('topright',inset=c(-0.3, 0), title= "Depths", legend = Depth, col=1:length(Depth), pch=16, cex =0.6, )
abline(h=1, lty = 2)

matplot(t(Probs$AngleDist),t(Probs$Trues), xlab="Arthwrd angle [deg]" , ylab="", main = "",
        xaxt="n", type = "l", lty = 1, col= 1:length(Depth),  pch = 16)
axis(1, at = seq(-10, 10, by = 2))
abline(h=1, lty = 2)
legend('topright',inset=c(-0.3, 0), title= "Depths", legend = Depth, col=1:length(Depth), pch=16, cex =0.6, )

#image(rotate(ProcWindow), asp =1)
#contour(rotate(ProcWindow), add = TRUE)
par(xpd=FALSE)
Perc100 <- vector()
Perc100[1] <- -Probs$AxisDist[1,max(which(Probs$Trues[1,]==1))]*pixel #error will show 
Depth_range = (Depth[length(Depth)]-Depth[1])/(length(Depth)-1)
plot(abs(Probs$AxisDist[1,])*pixel,Probs$Trues[1,]*Depth_range-min(Depth)-1,
     type = "l", xlim = c(0,max(Probs$AxisDist,na.rm = TRUE)*pixel), 
     ylim = c(-max(Depth)-1,-min(Depth)), xlab ="distance from axis (m)", ylab = "depth (m)") 
abline(h = c(-min(Depth),-min(Depth)-1), col = "grey76")
  abline(v=0)
for (i in 2:(length(Depth))) {
y = Probs$Trues[i,]*Depth_range-Depth_range*(i-1)-min(Depth)-1
lines(abs(Probs$AxisDist[i,])*pixel,y, type = "l")
abline(h = min(y, na.rm = TRUE), col = "grey76")

Perc100[i] <- -Probs$AxisDist[i,max(which(Probs$Trues[i,]==1))]*pixel
}
lines(Perc100,-Depth, col = 2)



#plot total probability at depth (when assuming 3D wedge sampling volume)
names(Probs$Tot) <- Depth
plot(Depth,Probs$Tot*100, type="l", ylim = c(0,100))
abline(h=c(0,100), lty=2)


matrix(round(Probs$Tot*100), dimnames = list(Depth,"Prob")) #total probability by depth
MaxProb <- apply(Probs$Trues,1,max, na.rm=T)
names(MaxProb) <- Depth
FirstDect <- min(as.numeric(names(MaxProb > 0)))


if (SaveMethods){
  
  XPixelPos <- seq(0,1,length.out = dim(Sheet)[2])
  XaxPos <- XPixelPos[seq(r_Pix+1, by = NextPing_Px,  length.out = AllPings)]
  XaxLab <- round(seq(0,by = NextPing_m, length.out = AllPings),1)  
  YUnit <- 1/dim(Sheet)[1]
  YPixelPos <- seq(0,1,length.out = dim(Sheet)[1])
  YaxLab <- -(round(r)):round(r)
  YaxPos <-  (r_Pix*pixel+YaxLab)*0.5/(r_Pix*pixel)
  
  png(paste(DirFig,"FullPattern.png",sep=""), width = 15, height = 6.25, units = "cm", res = 300)
par(mar = c(3,3.5,0.5,0))
image(rotate(Sheet), xaxt ='n', yaxt = 'n', xlab = "", ylab = "")  

axis(side = 1, at = XaxPos, labels = XaxLab)
axis(side = 1, at = XaxPos[seq(1,length(XaxPos),by=2)], labels = XaxLab[seq(1,length(XaxPos),by=2)])
axis(side = 1, at = XaxPos[seq(2,length(XaxPos-1),by=2)], labels = rep(NA,length((XaxPos-1)/2)))
axis(side = 2, at = YaxPos, labels = YaxLab, las = 2) 
axis(side = 2, lwd.tick=0, labels=FALSE)
mtext(side = 1, line = 2, "Distance travelled (m)")
mtext(side = 2, line = 2.7, "Athwartship distance (m)")

# axis(side = 1, line = 2, at = seq(0,1,length.out = dim(Sheet)[2]), labels =1:dim(Sheet)[2]) 
# axis(side = 2, line = 2, at = seq(0,1,length.out = dim(Sheet)[1]), labels =1:dim(Sheet)[1]) 
# abline(v = seq(0,1,length.out = dim(Sheet)[2]))
# abline(h = seq(0,1,length.out = dim(Sheet)[1]))
dev.off()  


png(paste(DirFig,"PatternProc.png",sep=""), width = 15, height = 6.25, units = "cm", res = 300)
par(mar = c(3,3.5,0.5,0))
image(rotate(PatternProc), asp =1, xaxt ='n', yaxt = 'n', xlab = "", ylab = "")
dev.off()

png(paste(DirFig,"BitmapProc.png",sep=""), width = 15, height = 6.25, units = "cm", res = 300)
par(mar = c(3,3.5,0.5,0))
image(rotate(BitmapProc), asp =1, xaxt ='n', yaxt = 'n', xlab = "", ylab = "")
dev.off()

{
png(paste(DirFig,"ProbEg.png",sep=""), width = 15, height = 5, units = "cm", res = 300)
par(mar=c(3, 3.5, 0.5, 0), oma=c(0,0,0,0))
plot(t(Probs$AxisDist)*pixel,t(Probs$Trues)*100, yaxt = 'n', xlab="" , ylab="",type = "l")
#abline(h=c(0,100), lty = 2, col = "grey60")
#lines(t(Probs$AxisDist)*pixel,t(Probs$Trues)*100, yaxt = 'n', xlab="" , ylab="",type = "l")

axis(side = 2, las =2)
mtext(side = 1, line = 2, "Athwartship distance (m)")
mtext(side = 2, line = 2.5, "Probability (%)")

dev.off()
}

}

if (SaveResultsRainbow) {
  
  ## rainbow plot
  png(paste(DirFig,"Res_Rainbow.png",sep=""), width = 15.92, height = 12, units = "cm", res = 300)
  par(mar=c(3.5, 3.4, 0.2, 3.5), xpd=T)
  par(oma=c(0,0,0,0))
  par(mfrow = c(2,1))
  matplot(t(Probs$AxisDist)*pixel,t(Probs$Trues)*100, xlab="" ,xaxt="n", yaxt ='n', ylab="", main = "",
          type = "l", lty = 1, lwd = 1.7, col= 1:length(Depth))
  axis(side=1, at = seq(-8,8, by = 2))
  axis(side=2, las=2)
  mtext(side=1, line = 2, "Athwartship distance (m)")
  mtext(side =2, line = 2.5, "Probability (%)")
  
  legend('topright',inset=c(-0.14, 0), title= "Range (m)", legend = Depth, 
         col=1:length(Depth), lty=1, lwd = 1.7, seg.len = 0.8, cex =0.735 )
  #abline(h=1, lty = 2)
  
  par(mar=c(3, 3.4, 0.7, 3.5), xpd=T)
  matplot(t(Probs$AngleDist),t(Probs$Trues)*100, xlab="" , yaxt ='n', ylab="", main = "",
          xaxt="n", type = "l", lty = 1, lwd = 1.7, col= 1:length(Depth))
  axis(1, at = seq(-8, 8, by = 2))
  axis(side=2, las=2)
  mtext(side=1, line = 2, "Athwartship angle (deg)")
  mtext(side =2, line = 2.5, "Probability (%)")
  #abline(h=1, lty = 2)
  legend('topright',inset=c(-0.14, 0), title= "Range (m)", legend = Depth, 
         col=1:length(Depth), lty=1, lwd = 1.7, seg.len = 0.9, cex =0.735 )
  dev.off()
}
  

  if (SaveResultsLine) {
  #line plot
  #plot total probability at depth (when assuming 3D wedge sampling volume)
  png(paste(DirFig,"Res_Line.png",sep=""), width = 15.92, height = 9, units = "cm", res = 300)
  par(mar=c(3, 3.5, 0.2, 0))
  
  plot(Depth,Probs$Tot*100, type="l", xaxt ="n", yaxt ="n", xlab = "", ylab = "", ylim = c(0,100),
       lwd = 2)
  segments(x0=0, y0=0, x1 = 5, y1 = 0,lwd = 2)
  abline(h=c(0,100), lty=2)
  axis(side = 1)
  axis(side = 2, las= 2)
  mtext(side = 1, line = 2, "Range (m)")
  mtext(side = 2, line = 2.5, "Simulation / wedge (%)")
  dev.off()
  
}


  MaxProb  #all maximum probabilities for each depth
cat("The first possible detection within the dataset is at depth", FirstDect,
    "m with a max probability of", MaxProb[as.character(FirstDect)],"m\n")

cat("The first depth with 100% detection is", 
    Depth[which(Perc100==min(Perc100, na.rm = TRUE))],"m\n")

# plot(Perc100_3,-Depth, type = "l", xlab = "distance from axis (m)", ylab = "depth (m)")
# lines(Perc100_4,-Depth)
# lines(Perc100_5,-Depth)



end_time <- Sys.time()
end_time - start_time

beep()
