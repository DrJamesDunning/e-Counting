
source("bresenHCirc.R")
wd <- ("C:/Users/james/Documents/PhD/Teunis PhD followup/Code/Echo Counting Simulation/For Github/Dunning et al")
          
setwd(wd)
load("Dunning et al_Data.RData")
yax <- seq(-Depth, 0, length = Ndepth)


# horizontal cross sections (circles)
png("10 m Circles.png", width = 9, height = 8.9, units = "cm", res = 300)
par(mar = c(3,4.1,0,0.25))  #default mar=c(5.1, 4.1, 4.1, 2.1)

Maxradius <- 16          # = 16 if 10 m deep, 18? beam angle
Circle <- bresenHCirc(Maxradius, Centre, fill = TRUE, plot = F)
image(Circle, asp=1, col =1, xaxt ="n", yaxt ="n", xlim = c(-0.015,1.015), ylim = c(-0.015,1.015))

for (i in 1:Maxradius) {
Circle <- bresenHCirc(Maxradius-i, Centre, fill = TRUE, plot = F)
if (((Maxradius-i) %% 2) == 1) {image(Circle, asp=1, col ="grey55", add = T) #if odd
} else                         {image(Circle, asp=1, col =1, add = T) }    #if even
}

#axis(side=1, at = seq(0,1,by=1/32*4), lab = seq(-16,16, by=4))
#axis(side=2, at = seq(0,1,by=1/32*4), lab = seq(-16,16, by=4), las=2)
#mtext(side=1, line = 2, "Radius (dm)", cex = )

dev.off()


#Survey matrix
#image3d(AllPings_n_fish, col = c("grey80","red2","#006400"), alpha.power = 0.5)  # 3D summary of simulation. col = pings, undetected fish, detected fish


{# MED, SED, Track Echograms
png("Echograms.png", width = 15.92, height = 8, units = "cm", res = 300)
par(mar = c(4,1.5,0.2,0.25))  #default mar=c(5.1, 4.1, 4.1, 2.1)
par(oma = c(0,2.1,0,0))
par(mfrow=c(1,3))

image(1:Npings, yax, rotate(Echogram),axes = FALSE, xlab="",ylab="") 
abline(v = c(0.5:Npings))
#abline(h = seq(-Depth, 0, by = 1))
axis(1, at = seq(1, Npings, by = 1), cex.axis = 1.2)
axis(2, at = seq(-Depth, 0, by = 2), labels = rev(seq(0,Depth, by =2)), las=1)
mtext(side = 1, text = "Ping no.", line = 2.7, cex=1.5)
mtext(side = 2, text = "Depth (m)", line = 2, cex=1.5)

#SED Echogram
image(1:Npings, yax, rotate(EchoSingleT),axes = FALSE, xlab="",ylab="", col = c("#FFFFC8","#F5A100")) 
abline(v = c(0.5:Npings))
axis(1, at = seq(1, Npings, by = 1))
axis(2, at = seq(-Depth, 0, by = 2), labels = rev(seq(0,Depth, by =2)), las=1)
mtext(side = 1, text = "Ping no.", line = 2.7, cex=1.5)
#mtext(side = 2, text = "Depth (m)", line = 2.5, cex=1.5)

#Track echogram
#Fish Tracks Echogram
image(1:Npings, yax, rotate(EchoTracks),axes = FALSE, xlab="",ylab="", col = sample(rainbow(length(TrackIDs)))) 
abline(v = c(0.5:Npings))
axis(1, at = seq(1, Npings, by = 1))
axis(2, at = seq(-Depth, 0, by = 2), labels = rev(seq(0,Depth, by =2)), las=1)
mtext(side = 1, text = "Ping no.", line = 2.7, cex=1.5)
#mtext(side = 2, text = "Depth (m)", line = 2.5, cex=1.5)

dev.off()}
par(mfrow=c(1,1))




## Long simulation results: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  #####
# LongDataDir <- paste("R:/CLSM/School of Biological Sciences/FisheriesGroup/PhD students/James Dunning/",
#                        "James Code/Echo Counting Simulation/", sep="")
# options(scipen = 1)
# 
# #save relevant Layer Data
# DataLayer <- "6Depths_6FishD_Lay10-40_500_v6v7.RData"  #with fish layer
# 
# load(paste(LongDataDir,DataLayer, sep=""))
# 
# ResL <- Res
# TimeLayer <- end_time - start_time
# 
# #save relevant FullWater Data
# DataFullWater <- "6Depths_6FishD_500_v7v8_Ella.RData"   # TrkDens-line with MTL
# load(paste(LongDataDir,DataFullWater, sep=""))
# Res2 <- Res
# 
# DataFullWater <- "6Depths_6FishD_500_v6v7.RData"        #TrkDens-line without MTL
# load(paste(LongDataDir,DataFullWater, sep=""))
# Res$ArD_SED_L_trk <- Res2$ArD_SED_trk
# 
# 
# #Double Standard errors and make them large enough for minimum graphical appearance
# Res$ArD_MED[,,c(3,4)]        <- Res$ArD_MED[,,c(3,4)]        *1.96
# Res$ArD_MED_RmOff[,,c(3,4)]  <- Res$ArD_MED_RmOff[,,c(3,4)]  *1.96
# Res$ArD_MED_trk[,,c(3,4)]    <- Res$ArD_MED_trk[,,c(3,4)]    *1.96
# Res$ArD_MED_trkYul[,,c(3,4)] <- Res$ArD_MED_trkYul[,,c(3,4)] *1.96
# Res$ArD_SED[,,c(3,4)]        <- Res$ArD_SED[,,c(3,4)]        *1.96
# Res$ArD_SED_RmOff[,,c(3,4)]  <- Res$ArD_SED_RmOff[,,c(3,4)]  *1.96
# Res$ArD_SED_trk[,,c(3,4)]    <- Res$ArD_SED_trk[,,c(3,4)]    *1.96
# Res$ArD_SED_trkYul[,,c(3,4)] <- Res$ArD_SED_trkYul[,,c(3,4)] *1.96
# Res$ArD_SED_L_trk[,,c(3,4)]  <- Res$ArD_SED_L_trk[,,c(3,4)]  *1.96
# Res$ArD_SED_L_trkYul[,,c(3,4)] <- Res$ArD_SED_L_trkYul[,,c(3,4)]  *1.96
# 
# ResL$ArD_MED[,,c(3,4)]       <- ResL$ArD_MED[,,c(3,4)]       *1.96
# ResL$ArD_MED_Yul[,,c(3,4)]   <- ResL$ArD_MED_Yul[,,c(3,4)]   *1.96
# Res$SED_loss[,,c(3,4)]       <- Res$SED_loss[,,c(3,4)]       *1.96
# 
# #delete layer data at 10 m bottom depth (only last voxel layer had fish and at 9.9 m depth, not within 10-40 m)
# ResL$ArD_MED[2,,] <- NA
# ResL$ArD_MED_Yul[2,,] <- NA
# ResL$ArD_MED_trk[2,,] <- NA
# ResL$ArD_MED_trkYul[2,,] <- NA
# save.image(file=paste(LongDataDir,"Ch2_Data.RData", sep=""))

#LongDataDir <- ("C:/Users/james/Documents/PhD/Teunis PhD followup/Code/Echo Counting Simulation/For Github/Dunning et al/")
#load(paste(LongDataDir,"Dunning et al_Data.RData", sep=""))

### Analysis ---------------------------------------------------------
options(scipen = 1)

#useful numbers
end_time - start_time + TimeLayer    #Time to run both simulations
mean(Res$TruArD[,,1]/Res$ArD_MED_RmOff[,,1]) #mean ratio between true density and EchoDens-Rm
1/(mean(Res$ArD_MED_RmOff[,,2])/100+1)       #same as before but calculated from %error
mean(Res$ArD_MED_RmOff[,,2])                 #mean of EchoDens-Rm %error

abs(Res$ArD_SED[,,2]) - abs(Res$ArD_SED_trk[,,2])   #Difference in %errors between EchoDens-point anTrkDens-line (+ve = smaller track error)
max(abs(Res$ArD_SED[,,2]) - abs(Res$ArD_SED_trk[,,2])) #max of previous line
Res$ArD_SED_L_trk[,,2]                                 #column and row names

mean(Res$ArD_MED_trkYul[1,,2])  #mean bias at 5m  bottom depth of TrkDens_Yule (MED)
mean(Res$ArD_MED_trkYul[6,,2])  #mean bias at 80m bottom depth of TrkDens_Yule (MED)

#significance of % Errors: (rows = depth (5,10,20,40,60,80) , columns = density (0.0005, 0.001, 0.005, 0.01,0.1,1))
sum((Res$ArD_MED[,,2]-Res$ArD_MED[,,4]) <= 0 & (Res$ArD_MED[,,2]+Res$ArD_MED[,,4]) >= 0)                               #EchoDens-point (MED)
sum((Res$ArD_MED_RmOff[,,2]-Res$ArD_MED_RmOff[,,4]) <= 0 & (Res$ArD_MED_RmOff[,,2]+Res$ArD_MED_RmOff[,,4]) >= 0)       #EchoDens-Rm    (MED)
sum((Res$ArD_MED_trk[,,2]-Res$ArD_MED_trk[,,4]) <= 0 & (Res$ArD_MED_trk[,,2]+Res$ArD_MED_trk[,,4]) >= 0)               #TrkDens-Point   (MED)
sum((Res$ArD_MED_trkYul[,,2]-Res$ArD_MED_trkYul[,,4]) <= 0 & (Res$ArD_MED_trkYul[,,2]+Res$ArD_MED_trkYul[,,4]) >= 0)   #TrkDens-Yule   (MED)
sum((Res$ArD_MED_Yul[,,2]-Res$ArD_MED_Yul[,,4]) <= 0 & (Res$ArD_MED_Yul[,,2]+Res$ArD_MED_Yul[,,4]) >= 0)               #EchoDens-Layer (MED)

sum((Res$ArD_SED[,,2]-Res$ArD_SED[,,4]) <= 0 & (Res$ArD_SED[,,2]+Res$ArD_SED[,,4]) >= 0)                               #EchoDens-point (SED)
sum((Res$ArD_SED_RmOff[,,2]-Res$ArD_SED_RmOff[,,4]) <= 0 & (Res$ArD_SED_RmOff[,,2]+Res$ArD_SED_RmOff[,,4]) >= 0)       #EchoDens-Rm    (SED)
sum((Res$ArD_SED_trk[,,2]-Res$ArD_SED_trk[,,4]) <= 0 & (Res$ArD_SED_trk[,,2]+Res$ArD_SED_trk[,,4]) >= 0)               #TrkDens-Point   (SED)
sum((Res$ArD_SED_trkYul[,,2]-Res$ArD_SED_trkYul[,,4]) <= 0 & (Res$ArD_SED_trkYul[,,2]+Res$ArD_SED_trkYul[,,4]) >= 0)   #TrkDens-Yule   (SED)
sum((Res$ArD_SED_L_trk[,,2]-Res$ArD_SED_L_trk[,,4]) <= 0 & (Res$ArD_SED_L_trk[,,2]+Res$ArD_SED_L_trk[,,4]) >= 0)             #TrkDens-Point   (SED+MTL)
sum((Res$ArD_SED_L_trkYul[,,2]-Res$ArD_SED_L_trkYul[,,4]) <= 0 & (Res$ArD_SED_L_trkYul[,,2]+Res$ArD_SED_L_trkYul[,,4]) >= 0) #TrkDens-Yule   (SED+MTL)

sum((ResL$ArD_MED[,,2]-ResL$ArD_MED[,,4]) <= 0 & (ResL$ArD_MED[,,2]+ResL$ArD_MED[,,4]) >= 0,na.rm=T)                         #EchoDens-point (MED w Layer)
sum((ResL$ArD_MED_Yul[,,2]-ResL$ArD_MED_Yul[,,4]) <= 0 & (ResL$ArD_MED_Yul[,,2]+ResL$ArD_MED_Yul[,,4]) >= 0,na.rm=T) #EchoDens-Layer (MED w Layer)
sum((Res$SED_loss[,,1]-Res$SED_loss[,,3]) <= 0 & (Res$SED_loss[,,1]+Res$SED_loss[,,3]) >= 0)                         #SED loss




# MED comparison 
{png("MED Comparisons.png", width = 15.92, height = 10, units = "cm", res = 300)
par(mfrow=c(2,2))
par(mar = c(1,1,0,0.1))  #default mar=c(5.1, 4.1, 4.1, 2.1)
par(oma = c(2.3,2.7,3.5,0))

Xtitle <- "Bottom depth (m)"
YLim <- c(-100,100)
Ytitle <- "Mean %Error"
TitleX <- 21

{
Title <- expression("A) EchoDens"["Point"])
matplot(MovingVar1,Res$ArD_MED[,,2], ylim = YLim, xaxt ='n', xlab='', yaxt='n', ylab='',
        type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
abline(a=0,b=0, lty=2)
arrows(MovingVar1, Res$ArD_MED[,,2]-Res$ArD_MED[,,4], MovingVar1, Res$ArD_MED[,,2]+Res$ArD_MED[,,4],
       length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
mtext(side = 2, line = 2.7,  Ytitle)
mtext(side = 3, line = -1.1, at = TitleX-3, Title, cex = 0.8)
axis(side = 1,  at = c(5,10,20,40,60,80), labels = F)
axis(side = 2, las = 2)
  

Title <- expression("B) EchoDens"["Rm"])
matplot(MovingVar1,Res$ArD_MED_RmOff[,,2], ylim = YLim, xaxt ='n', yaxt ='n', xlab='', ylab='',
        type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
abline(a=0,b=0, lty=2)
arrows(MovingVar1, Res$ArD_MED_RmOff[,,2]-Res$ArD_MED_RmOff[,,4], MovingVar1, Res$ArD_MED_RmOff[,,2]+Res$ArD_MED_RmOff[,,4],
       length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1))) #Add SE
mtext(side = 3, line = -1.1, at = TitleX-4, Title, cex = 0.8)
axis(side = 1,  at = c(5,10,20,40,60,80), labels =F)
axis(side = 2, labels = F)


Title <- expression("C) TrkDens"["Point"])
matplot(MovingVar1,Res$ArD_MED_trk[,,2], ylim = YLim ,xaxt ='n', xlab='', yaxt ='n', ylab='',
        type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
abline(a=0,b=0, lty=2)
arrows(MovingVar1, Res$ArD_MED_trk[,,2]-Res$ArD_MED_trk[,,4], MovingVar1, Res$ArD_MED_trk[,,2]+Res$ArD_MED_trk[,,4],
       length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
mtext(side = 1, line = 2.2, Xtitle)
mtext(side = 2, line = 2.7, Ytitle)
mtext(side = 3, line = -1.1, at = TitleX-5, Title, cex = 0.8)
axis(side = 1,  at = c(5,10,20,40,60,80))
axis(side = 2,las = 2)


Title <- expression("D) TrkDens"["Yule"])
matplot(MovingVar1,Res$ArD_MED_trkYul[,,2], ylim = YLim ,xaxt ='n', xlab='', yaxt ='n', ylab='',
        type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
abline(a=0,b=0, lty=2)
arrows(MovingVar1, Res$ArD_MED_trkYul[,,2]-Res$ArD_MED_trkYul[,,4], MovingVar1, Res$ArD_MED_trkYul[,,2]+Res$ArD_MED_trkYul[,,4],
       length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
mtext(side = 1, line = 2.2, Xtitle)
mtext(side = 3, line = -1.1, at = TitleX-5, Title, cex = 0.8)
axis(side = 1,  at = c(5,10,20,40,60,80))
axis(side = 2, labels =F)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
text(0.1, 1.01, "MED data", cex = 1.5)
legend(-0.37,1.035,horiz = T,title= "", xpd = TRUE,
       legend = MovingVar2, col=1:length(MovingVar2), 
       pch=16, lty = 1, seg.len = 0.8, lwd=1.5,  cex =1, x.intersp = 0.5, text.width = 0.2,
       bty='n')
text(-0.65,0.857, expression(Dens[Input]~"(Fish/m"^3*"):"))
rect(xleft = -0.868, xright =1.075, ybottom = 0.79, ytop=0.925 )

}
dev.off()
}




# SED comparison 
{png("SED Comparisons.png", width = 15.92, height = 14, units = "cm", res = 300)
  par(mfrow=c(3,2))
  par(mar = c(1,1,0,0.1))  #default mar=c(5.1, 4.1, 4.1, 2.1)
  par(oma = c(2.3,2.7,3.5,0))
  
  YLim <- c(-100,100)
  Ytitle <- "Mean %Error"
  TitleX <- 21
  
  {
    Title <- expression("A) EchoDens"["Point"])
    matplot(MovingVar1,Res$ArD_SED[,,2], ylim = YLim, xaxt ='n', xlab='', yaxt='n', ylab='',
            type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
    abline(a=0,b=0, lty=2)
    arrows(MovingVar1, Res$ArD_SED[,,2]-Res$ArD_SED[,,4], MovingVar1, Res$ArD_SED[,,2]+Res$ArD_SED[,,4],
           length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
    mtext(side = 2, line = 2.5,  Ytitle)
    mtext(side = 3, line = -1.3, at = TitleX-3, Title, cex = 0.8)
    axis(side = 1,  at = c(5,10,20,40,60,80), labels = F)
    axis(side = 2, las = 2)
    
    
    Title <- expression("B) EchoDens"["Rm"])
    matplot(MovingVar1,Res$ArD_SED_RmOff[,,2], ylim = YLim, xaxt ='n', yaxt ='n', xlab='', ylab='',
            type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
    abline(a=0,b=0, lty=2)
    arrows(MovingVar1, Res$ArD_SED_RmOff[,,2]-Res$ArD_SED_RmOff[,,4], MovingVar1, Res$ArD_SED_RmOff[,,2]+Res$ArD_SED_RmOff[,,4],
           length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1))) #Add SE
    mtext(side = 3, line = -1.3, at = TitleX-4, Title, cex = 0.8)
    axis(side = 1,  at = c(5,10,20,40,60,80), labels =F)
    axis(side = 2, labels = F)
    
    
    Title <- expression("C) TrkDens"["Point"])
    matplot(MovingVar1,Res$ArD_SED_trk[,,2], ylim = YLim ,xaxt ='n', xlab='', yaxt ='n', ylab='',
            type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)

    abline(a=0,b=0, lty=2)
    arrows(MovingVar1, Res$ArD_SED_trk[,,2]-Res$ArD_SED_trk[,,4], MovingVar1, Res$ArD_SED_trk[,,2]+Res$ArD_SED_trk[,,4],
           length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
    mtext(side = 2, line = 2.5, Ytitle)
    mtext(side = 3, line = -1.3, at = TitleX-5, Title, cex = 0.8)
    axis(side = 2,las = 2)
    
    
    Title <- expression("D) TrkDens"["Yule"])
    matplot(MovingVar1,Res$ArD_SED_trkYul[,,2], ylim = YLim ,xaxt ='n', xlab='', yaxt ='n', ylab='',
            type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
    abline(a=0,b=0, lty=2)
    arrows(MovingVar1, Res$ArD_SED_trkYul[,,2]-Res$ArD_SED_trkYul[,,4], MovingVar1, Res$ArD_SED_trkYul[,,2]+Res$ArD_SED_trkYul[,,4],
           length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
    mtext(side = 3, line = -1.3, at = TitleX-5, Title, cex = 0.8)
    axis(side = 2, labels =F)
    
    
    Title <- expression("E) TrkDens"["Point"] + MTL)
    matplot(MovingVar1,Res$ArD_SED_L_trk[,,2], ylim = YLim ,xaxt ='n', xlab='', yaxt ='n', ylab='',
            type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
    
    abline(a=0,b=0, lty=2)
    arrows(MovingVar1, Res$ArD_SED_L_trk[,,2]-Res$ArD_SED_L_trk[,,4], MovingVar1, Res$ArD_SED_L_trk[,,2]+Res$ArD_SED_L_trk[,,4],
           length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
    mtext(side = 1, line = 2.2, Xtitle)
    mtext(side = 2, line = 2.5, Ytitle)
    mtext(side = 3, line = -1.3, at = TitleX-0, Title, cex = 0.8)
    axis(side = 1,  at = c(5,10,20,40,60,80))
    axis(side = 2,las = 2)
    
    
    Title <- expression("F) TrkDens"["Yule"] + MTL)
    matplot(MovingVar1,Res$ArD_SED_L_trkYul[,,2], ylim = YLim ,xaxt ='n', xlab='', yaxt ='n', ylab='',
            type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
    abline(a=0,b=0, lty=2)
    arrows(MovingVar1, Res$ArD_SED_L_trkYul[,,2]-Res$ArD_SED_L_trkYul[,,4], MovingVar1, Res$ArD_SED_L_trkYul[,,2]+Res$ArD_SED_L_trkYul[,,4],
           length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
    mtext(side = 1, line = 2.2, Xtitle)
    mtext(side = 3, line = -1.3, at = TitleX-0, Title, cex = 0.8)
    axis(side = 1,  at = c(5,10,20,40,60,80))
    axis(side = 2, labels =F)
    
    
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    text(0.1, 1.05, "SED data", cex = 1.5)
    legend(-0.37,1.042,horiz = T,title= "", xpd = TRUE,
           legend = MovingVar2, col=1:length(MovingVar2), 
           pch=16, lty = 1, seg.len = 0.8, lwd=1.5,  cex =1, x.intersp = 0.5, text.width = 0.2,
           bty='n')
    text(-0.65,0.94, expression(Dens[Input]~"(Fish/m"^3*"):"))
    rect(xleft = -0.912, xright =1.076, ybottom = 0.9, ytop=0.978 )
    
  }
  dev.off()
}



# Layer comparison
{png("Layer Comparisons.png", width = 15.92, height = 10, units = "cm", res = 300)
  par(mfrow=c(2,2))
  par(mar = c(1,1,0,0.1))  #default mar=c(5.1, 4.1, 4.1, 2.1)
  par(oma = c(2.3,2.7,3.5,0))

YLim <- c(-100,100)
Ytitle <- "Mean %Error"
TitleX <- 21
  
  
Title <- expression("A) EchoDens"["Point"])
matplot(MovingVar1,ResL$ArD_MED[,,2], ylim = YLim, xaxt ='n', xlab='', yaxt='n', ylab='',
        type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
abline(a=0,b=0, lty=2)
arrows(MovingVar1, ResL$ArD_MED[,,2]-ResL$ArD_MED[,,4], MovingVar1, ResL$ArD_MED[,,2]+ResL$ArD_MED[,,4],
       length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
mtext(side = 2, line = 2.7,  Ytitle)
mtext(side = 3, line = -1.1, at = TitleX+45, Title, cex = 0.8)
axis(side = 2, las = 2)


Title <- expression("B) EchoDens"["Layer"])
matplot(MovingVar1,ResL$ArD_MED_Yul[,,2], ylim = YLim, xaxt ='n', yaxt ='n', xlab='', ylab='',
        type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
abline(a=0,b=0, lty=2)
arrows(MovingVar1, ResL$ArD_MED_Yul[,,2]-ResL$ArD_MED_Yul[,,4], MovingVar1, ResL$ArD_MED_Yul[,,2]+ResL$ArD_MED_Yul[,,4],
       length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1))) #Add SE
mtext(side = 3, line = -1.2, at = TitleX+45, Title, cex = 0.8)
axis(side = 2, labels = F)


Title <- expression("C) TrkDens"["Point"])
matplot(MovingVar1,ResL$ArD_MED_trk[,,2], ylim = YLim, xaxt ='n', xlab='', yaxt='n', ylab='',
        type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
abline(a=0,b=0, lty=2)
arrows(MovingVar1, ResL$ArD_MED_trk[,,2]-ResL$ArD_MED_trk[,,4], MovingVar1, ResL$ArD_MED_trk[,,2]+ResL$ArD_MED_trk[,,4],
       length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
mtext(side = 1, line = 2.2, Xtitle)
mtext(side = 2, line = 2.7,  Ytitle)
mtext(side = 3, line = -1.1, at = TitleX+45, Title, cex = 0.8)
axis(side = 1,  at = c(20,40,60,80))
axis(side = 2, las = 2)


Title <- expression("D) TrkDens"["Yule"])
matplot(MovingVar1,ResL$ArD_MED_trkYul[,,2], ylim = YLim, xaxt ='n', yaxt ='n', xlab='', ylab='',
        type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
abline(a=0,b=0, lty=2)
arrows(MovingVar1, ResL$ArD_MED_trkYul[,,2]-ResL$ArD_MED_trkYul[,,4], MovingVar1, ResL$ArD_MED_trkYul[,,2]+ResL$ArD_MED_trkYul[,,4],
       length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1))) #Add SE
mtext(side = 1, line = 2.2, Xtitle)
mtext(side = 3, line = -1.1, at = TitleX+45, Title, cex = 0.8)
axis(side = 1,  at = c(20,40,60,80))
axis(side = 2, labels = F)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
text(0.1, 1.01, "MED w/ Fish Layer", cex = 1.5)
legend(-0.37,1.035,horiz = T,title= "", xpd = TRUE,
       legend = MovingVar2, col=1:length(MovingVar2), 
       pch=16, lty = 1, seg.len = 0.8, lwd=1.5,  cex =1, x.intersp = 0.5, text.width = 0.2,
       bty='n')
text(-0.65,0.857, expression(Dens[Input]~"(Fish/m"^3*"):"))
rect(xleft = -0.868, xright =1.075, ybottom = 0.79, ytop=0.925 )

dev.off()
  }

