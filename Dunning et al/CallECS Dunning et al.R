#This script calls and loops the program ECountingSim.V3.R
#you can make one and only one of the user input variables as a vector, the others need to be constant. 
#Then paste the name of the vector variable on the right side of "MovingVar1 <-!

#install.packages(c("abind"))
#install.packages(c("rgl"))
#install.packages(c("misc3d"))

# Clean workspace
rm(list=ls(all=TRUE))
#dev.off()

start_time <- Sys.time()

library(abind)
library(rgl)
library(misc3d)


#setwd("F:/Greenland2 October/Echo-counting 4 Fletcher")
setwd ("R:/CLSM/School of Biological Sciences/FisheriesGroup/PhD students/James Dunning/James Code/Echo Counting Simulation")
source("bresenham.R")
source("bresenHCirc.R")

rotate <- function(x) t(apply(x, 2, rev))

progress <- function (x, max = 100) {
  percent <- x / max * 100
  cat(sprintf('\r[%-50s] %d%%',paste(rep('=', percent / 2), collapse = ''), floor(percent)))
  if (x == max) {cat('\n')}
}

se <- function(x) {
  out <- vector()
  for (indx in 1:dim(x)[2]) {
    out[indx] <- sd(x[,indx],na.rm = TRUE)/sqrt(sum(!is.na(x[,indx])))
  }
  return(out)
}


#---- user input
Var_Depth <-     c(5,10,20,40,60,80)       #depth in m
#Var_Depth <-     c(5,80)
Var_Distance <-  seq(50,50, by = 1)       #length of transect in m
Var_Speed_kts <- seq(6,6, by = 0.1)       #speed of boat in knots
Var_PingInt <-   seq(0.5,0.5, by = 0.1)       #ping interval in seconds
Var_Angle <-     seq(18,18, by = 1)       #beam angle in degrees
Var_MinTrack <-  seq(3,3, by = 1)        #minimum number of detections in a track (gaps allowed)
Var_FishD_m <-   c(0.0005, 0.001, 0.005, 0.01, 0.1, 1)    #fish density in fish/m^3
#Var_FishD_m <-   c(0.1, 1)

AddLayer  <- FALSE     #TRUE= use layer defined below, FALSE= disregard layer, have fish in all water column
FishLayer <- c(10,40)  #upper and lower limits of fish layer in m (the rest is empty water)

Var_Voxel <-     seq(0.1,0.1, by = 0.01)  #length of edge of voxel in m

Iterations <- 500         #how many times to repeat simulation for a given setting
MovingVar1 <- Var_Depth   #paste here one of two vector variables from list above
MovingVar2 <- Var_FishD_m #paste here one of two vector variables from list above

WriteReports <- TRUE       #Do you want to write reports? (useful for long processes)
ReportFreq <- 2            #how often to write a report (hrs)

Xtitle <- "Depth (m)"     #for making graphs
legTitle <- "(fish m^-3)"  #for making graphs
subtitle <- paste("Dist=",Var_Distance,"m; Speed=",Var_Speed_kts,"kts; Ping=",Var_PingInt,
                  "s; Angle=",Var_Angle,"°; voxel=",Var_Voxel,"m^3; AddLayer=",AddLayer,
                  "; Layer=",FishLayer[1],"-",FishLayer[2],"m; Iterations=",Iterations,sep="")

Graph <- FALSE           #if true, graphics for EACH simulation will be produced, much slower to run!


{
  
  #---- Check requirements, stop with error message if not met ########################
  
  if (min(Var_Distance) < (max(Var_Depth)*tan(((max(Var_Angle)*pi/180)/2)))*2 ) {
    stop(paste("Minimum distance too short to include greatest ping. Min distance =",min(Var_Distance), 
               "m ; Greatest ping =", (max(Var_Depth)*tan(((max(Var_Angle)*pi/180)/2)))*2)," m")
  }
  
  #---- start of program ################################################################
  NextReport <- ReportFreq
  
  cols <- length(MovingVar2)
  rows <- length(MovingVar1)
  DimNamesR <- list(MovingVar1,MovingVar2,c("Val","%Er","SEval","SE%er"))
  DimNamesT <- list(1:Iterations,c("Val","%Er")) #for itTable
  
  Res <- list(Npings = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              F_Sampled = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR), 
              F_Sampled_pc = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              MED_NTargets = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              SED_NTargets = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR), 
              SED_loss = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              FTracked = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              Track_loss = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              TruVolD = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),  
              VolD_MED = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              VolD_SED = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              AbsNTrue = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              AbsN_MED = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              AbsN_SED = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              TruArD = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_MED = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_MED_Rm = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_MED_RmOff = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_MED_Yul = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_SED = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_SED_Rm = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_SED_RmOff = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_SED_Yul = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_MED_trk = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_MED_trkYul = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_MED_trkX = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_MED_trkXOff = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_SED_trk = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_SED_trkYul = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_SED_L_trkYul = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_SED_trkX = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR),
              ArD_SED_trkXOff = array(rep(NA, rows*cols*4),c(rows, cols, 4), dimnames = DimNamesR)  )
  
  # for (i in 1:length(Res)) {
  #   Res[[i]] <- matrix(ncol = length(MovingVar2), nrow = length(MovingVar1))
  #   rownames(Res[[i]]) <- MovingVar1
  #   colnames(Res[[i]]) <- MovingVar2
  # }
  
  
  
  #Start of iteration loops ----
  
  for (v2 in 1:length(MovingVar2)) {
    if (setequal(MovingVar2,Var_Depth))     {Depth <- Var_Depth[v2]}         else {Depth <- Var_Depth}
    if (setequal(MovingVar2,Var_Distance))  {Distance <- Var_Distance[v2]}   else {Distance <- Var_Distance}
    if (setequal(MovingVar2,Var_Speed_kts)) {Speed_kts <- Var_Speed_kts[v2]} else {Speed_kts <- Var_Speed_kts}
    if (setequal(MovingVar2,Var_PingInt))   {PingInt <- Var_PingInt[v2]}     else {PingInt <- Var_PingInt}
    if (setequal(MovingVar2,Var_Angle))     {BeamAngle_Deg <- Var_Angle[v2]} else {BeamAngle_Deg <- Var_Angle}
    if (setequal(MovingVar2,Var_MinTrack))  {MinTrack <- Var_MinTrack[v2]}   else {MinTrack <- Var_MinTrack}
    if (setequal(MovingVar2,Var_FishD_m))   {FishD_m <- Var_FishD_m[v2]}     else {FishD_m <- Var_FishD_m}
    if (setequal(MovingVar2,Var_Voxel))     {Voxel <- Var_Voxel[v2]}         else {Voxel <- Var_Voxel}
    
    
    for (v1 in 1:length(MovingVar1)) {
      
      if (setequal(MovingVar1,Var_Depth))     {Depth <- Var_Depth[v1]}
      if (setequal(MovingVar1,Var_Distance))  {Distance <- Var_Distance[v1]}   
      if (setequal(MovingVar1,Var_Speed_kts)) {Speed_kts <- Var_Speed_kts[v1]} 
      if (setequal(MovingVar1,Var_PingInt))   {PingInt <- Var_PingInt[v1]}     
      if (setequal(MovingVar1,Var_Angle))     {BeamAngle_Deg <- Var_Angle[v1]} 
      if (setequal(MovingVar1,Var_MinTrack))  {MinTrack <- Var_MinTrack[v1]}   
      if (setequal(MovingVar1,Var_FishD_m))   {FishD_m <- Var_FishD_m[v1]}     
      if (setequal(MovingVar1,Var_Voxel))     {Voxel <- Var_Voxel[v1]}         
      
      if(AddLayer == FALSE) {
        FishLayer <- c(0,Depth)   #No layer, fish are randomly distributed across full water column
      }
      
      itTable <- list(F_Sampled = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      F_Sampled_pc = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      MED_NTargets = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      SED_NTargets = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      SED_loss = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      FTracked = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      Track_loss = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      TruVolD = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),  
                      VolD_MED = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      VolD_SED = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      AbsNTrue = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      AbsN_MED = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      AbsN_SED = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      TruArD = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_MED = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_MED_Rm = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_MED_RmOff = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_MED_Yul = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_SED = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_SED_Rm = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_SED_RmOff = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_SED_Yul = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_MED_trk = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_MED_trkYul = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_MED_trkX = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_MED_trkXOff = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),  
                      ArD_SED_trk = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_SED_trkYul = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_SED_L_trkYul = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_SED_trkX = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT),
                      ArD_SED_trkXOff = matrix(ncol = 2, nrow = Iterations, dimnames = DimNamesT)  )
      
      for (it in 1:Iterations) {
        
        ProgText <- paste("  Performing iteration",it,"of",Iterations, "for setting number",(v2-1)*length(MovingVar1)+v1,
                          "of",length(MovingVar1)*length(MovingVar2),"      ")
        cat(ProgText,'\r')
        Perc <- (Iterations*(v2-1)*length(MovingVar1)+Iterations*(v1-1)+it)/(Iterations*length(MovingVar1)*length(MovingVar2))*100
        progress(Perc)
        
        if (WriteReports & Sys.time() >= (start_time+NextReport*3600)) {
          writeLines(paste(ProgText,'\r',round(Perc),"% reached at ",Sys.time(),sep=""),".ProgressReport_Ella.txt")
          NextReport <- NextReport + ReportFreq
        }
        
        
        
        source("e-Counting Dunning et al.R")
        
        itTable$F_Sampled[it,]    <- MED$FSampled    # col 1 = estimate value, col 2 = % error (deviation from true value)
        itTable$F_Sampled_pc[it,] <- MED$FSampled_pc
        itTable$MED_NTargets[it,] <- MED$NTargets
        itTable$SED_NTargets[it,] <- SED$NTargets
        
        itTable$SED_loss[it,]   <- SED$SED_loss
        itTable$FTracked[it,]   <- SED$FTracked
        itTable$Track_loss[it,] <- SED$Track_loss
        
        itTable$TruVolD[it,]  <- MED$TruVolD
        itTable$VolD_MED[it,] <- MED$VolD
        itTable$VolD_SED[it,] <- SED$VolD
        
        itTable$AbsNTrue[it,] <- MED$TruNFish
        itTable$AbsN_MED[it,] <- MED$NFish
        itTable$AbsN_SED[it,] <- SED$NFish
        
        itTable$TruArD[it,]        <- MED$TruAreaD
        itTable$ArD_MED[it,]       <- MED$AreaD
        itTable$ArD_MED_Rm[it,]    <- MED$AreaDRm
        itTable$ArD_MED_RmOff[it,] <- MED$AreaDRmOff
        itTable$ArD_MED_Yul[it,]   <- MED$AreaDYul
        itTable$ArD_SED[it,]       <- SED$AreaD
        itTable$ArD_SED_Rm[it,]    <- SED$AreaDRm
        itTable$ArD_SED_RmOff[it,] <- SED$AreaDRmOff
        itTable$ArD_SED_Yul[it,]   <- SED$AreaDYul
        
        itTable$ArD_MED_trk[it,]     <- MED$AreaDTrack
        itTable$ArD_MED_trkYul[it,]  <- MED$AreaDtrackYul
        itTable$ArD_MED_trkX[it,]    <- MED$AreaDTrackX
        itTable$ArD_MED_trkXOff[it,] <- MED$AreaDTrackXOff
        itTable$ArD_SED_trk[it,]     <- SED$AreaDTrack
        itTable$ArD_SED_trkYul[it,]  <- SED$AreaDtrackYul
        itTable$ArD_SED_L_trkYul[it,]<- SED$AreaD_L_TrkYul
        itTable$ArD_SED_trkX[it,]    <- SED$AreaDTrackX
        itTable$ArD_SED_trkXOff[it,] <- SED$AreaDTrackXOff
        
      }
      
      # mean values and %errors
      Res$Npings[v1,v2,1]         <- Npings #page1=mean value, page2=mean %error, page3=SEv, page4=SEer
      Res$F_Sampled[v1,v2,1:2]    <- colMeans(itTable$F_Sampled, na.rm = TRUE)
      Res$F_Sampled_pc[v1,v2,1:2] <- colMeans(itTable$F_Sampled_pc, na.rm = TRUE)
      Res$MED_NTargets[v1,v2,1:2] <- colMeans(itTable$MED_NTargets, na.rm = TRUE)
      Res$SED_NTargets[v1,v2,1:2] <- colMeans(itTable$SED_NTargets, na.rm = TRUE)
      
      Res$SED_loss[v1,v2,1:2]   <- colMeans(itTable$SED_loss, na.rm = TRUE)
      Res$FTracked[v1,v2,1:2]   <- colMeans(itTable$FTracked, na.rm = TRUE)
      Res$Track_loss[v1,v2,1:2] <- colMeans(itTable$Track_loss, na.rm = TRUE)
      
      Res$TruVolD[v1,v2,1:2]  <- colMeans(itTable$TruVolD, na.rm = TRUE)
      Res$VolD_MED[v1,v2,1:2] <- colMeans(itTable$VolD_MED, na.rm = TRUE)
      Res$VolD_SED[v1,v2,1:2] <- colMeans(itTable$VolD_SED, na.rm = TRUE)
      
      Res$AbsNTrue[v1,v2,1:2] <- colMeans(itTable$AbsNTrue, na.rm = TRUE)
      Res$AbsN_MED[v1,v2,1:2] <- colMeans(itTable$AbsN_MED, na.rm = TRUE)
      Res$AbsN_SED[v1,v2,1:2] <- colMeans(itTable$AbsN_SED, na.rm = TRUE)
      
      Res$TruArD[v1,v2,1:2]        <- colMeans(itTable$TruArD, na.rm = TRUE)
      Res$ArD_MED[v1,v2,1:2]       <- colMeans(itTable$ArD_MED, na.rm = TRUE)
      Res$ArD_MED_Rm[v1,v2,1:2]    <- colMeans(itTable$ArD_MED_Rm, na.rm = TRUE)
      Res$ArD_MED_RmOff[v1,v2,1:2] <- colMeans(itTable$ArD_MED_RmOff, na.rm = TRUE)
      Res$ArD_MED_Yul[v1,v2,1:2]   <- colMeans(itTable$ArD_MED_Yul, na.rm = TRUE)
      Res$ArD_SED[v1,v2,1:2]       <- colMeans(itTable$ArD_SED, na.rm = TRUE)
      Res$ArD_SED_Rm[v1,v2,1:2]    <- colMeans(itTable$ArD_SED_Rm, na.rm = TRUE)
      Res$ArD_SED_RmOff[v1,v2,1:2] <- colMeans(itTable$ArD_SED_RmOff, na.rm = TRUE)
      Res$ArD_SED_Yul[v1,v2,1:2]   <- colMeans(itTable$ArD_SED_Yul, na.rm = TRUE)
      
      Res$ArD_MED_trk[v1,v2,1:2]     <- colMeans(itTable$ArD_MED_trk, na.rm = TRUE)
      Res$ArD_MED_trkYul[v1,v2,1:2]  <- colMeans(itTable$ArD_MED_trkYul, na.rm = TRUE)
      Res$ArD_MED_trkX[v1,v2,1:2]    <- colMeans(itTable$ArD_MED_trkX, na.rm = TRUE)
      Res$ArD_MED_trkXOff[v1,v2,1:2] <- colMeans(itTable$ArD_MED_trkXOff, na.rm = TRUE)
      Res$ArD_SED_trk[v1,v2,1:2]     <- colMeans(itTable$ArD_SED_trk, na.rm = TRUE)
      Res$ArD_SED_trkYul[v1,v2,1:2]  <- colMeans(itTable$ArD_SED_trkYul, na.rm = TRUE)
      Res$ArD_SED_L_trkYul[v1,v2,1:2]<- colMeans(itTable$ArD_SED_L_trkYul, na.rm = TRUE)
      Res$ArD_SED_trkX[v1,v2,1:2]    <- colMeans(itTable$ArD_SED_trkX, na.rm = TRUE)
      Res$ArD_SED_trkXOff[v1,v2,1:2] <- colMeans(itTable$ArD_SED_trkXOff, na.rm = TRUE)
      
      
      #Standard errors 
      Res$F_Sampled[v1,v2,3:4]    <- se(itTable$F_Sampled)    #page1=mean value, page2=mean %error, page3=SEv, page4=SEer
      Res$F_Sampled_pc[v1,v2,3:4] <- se(itTable$F_Sampled_pc)
      Res$MED_NTargets[v1,v2,3:4] <- se(itTable$MED_NTargets)
      Res$SED_NTargets[v1,v2,3:4] <- se(itTable$SED_NTargets)
      
      Res$SED_loss[v1,v2,3:4]   <- se(itTable$SED_loss)
      Res$FTracked[v1,v2,3:4]   <- se(itTable$FTracked)
      Res$Track_loss[v1,v2,3:4] <- se(itTable$Track_loss)
      
      Res$TruVolD[v1,v2,3:4]  <- se(itTable$TruVolD)
      Res$VolD_MED[v1,v2,3:4] <- se(itTable$VolD_MED)
      Res$VolD_SED[v1,v2,3:4] <- se(itTable$VolD_SED)
      
      Res$AbsNTrue[v1,v2,3:4] <- se(itTable$AbsNTrue)
      Res$AbsN_MED[v1,v2,3:4] <- se(itTable$AbsN_MED)
      Res$AbsN_SED[v1,v2,3:4] <- se(itTable$AbsN_SED)
      
      Res$TruArD[v1,v2,3:4]        <- se(itTable$TruArD)
      Res$ArD_MED[v1,v2,3:4]       <- se(itTable$ArD_MED)
      Res$ArD_MED_Rm[v1,v2,3:4]    <- se(itTable$ArD_MED_Rm)
      Res$ArD_MED_RmOff[v1,v2,3:4] <- se(itTable$ArD_MED_RmOff)
      Res$ArD_MED_Yul[v1,v2,3:4]   <- se(itTable$ArD_MED_Yul)
      Res$ArD_SED[v1,v2,3:4]       <- se(itTable$ArD_SED)
      Res$ArD_SED_Rm[v1,v2,3:4]    <- se(itTable$ArD_SED_Rm)
      Res$ArD_SED_RmOff[v1,v2,3:4] <- se(itTable$ArD_SED_RmOff)
      Res$ArD_SED_Yul[v1,v2,3:4]   <- se(itTable$ArD_SED_Yul)
      
      Res$ArD_MED_trk[v1,v2,3:4]     <- se(itTable$ArD_MED_trk)
      Res$ArD_MED_trkYul[v1,v2,3:4]  <- se(itTable$ArD_MED_trkYul)
      Res$ArD_MED_trkX[v1,v2,3:4]    <- se(itTable$ArD_MED_trkX)
      Res$ArD_MED_trkXOff[v1,v2,3:4] <- se(itTable$ArD_MED_trkXOff)
      Res$ArD_SED_trk[v1,v2,3:4]     <- se(itTable$ArD_SED_trk)
      Res$ArD_SED_trkYul[v1,v2,3:4]  <- se(itTable$ArD_SED_trkYul)
      Res$ArD_SED_L_trkYul[v1,v2,3:4]<- se(itTable$ArD_SED_L_trkYul)
      Res$ArD_SED_trkX[v1,v2,3:4]    <- se(itTable$ArD_SED_trkX)
      Res$ArD_SED_trkXOff[v1,v2,3:4] <- se(itTable$ArD_SED_trkXOff)
      
      
      
    }
    if (WriteReports & v1 == length(MovingVar1)) {
      CycleTime <- Sys.time()-start_time
      writeLines(paste(v2,"full cycles out of",length(MovingVar2),"completed in",round(CycleTime,1),units(CycleTime),'\r',
                       "Average time for one cycle", round((CycleTime/v2),1), units(CycleTime),'\r',
                       "ETA:", start_time+(CycleTime/v2)*length(MovingVar2)),".ETA_Processing_Ella.txt")
      writeLines(paste(ProgText,'\r',round(Perc),"% reached at ",Sys.time(),sep=""),".ProgressReport_Ella.txt")
    }
  }
  
  #Graphics ----
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  
  Title <- "Number of pings in transect"
  Ytitle <- "N#"
  matplot(MovingVar1,Res$Npings[,,1], xlab=Xtitle , ylab=Ytitle, main = Title,
          xlim = c(0,max(Var_Depth)), ylim = c(0,max(Res$Npings[,,1])),
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "% Fish Sampled"
  Ytitle <- "Mean % +/- SE"
  matplot(MovingVar1,Res$F_Sampled_pc[,,1], ylim = c(0,100), xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  arrows(MovingVar1, Res$F_Sampled_pc[,,1]-Res$F_Sampled_pc[,,3], MovingVar1, Res$F_Sampled_pc[,,1]+Res$F_Sampled_pc[,,3],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "% Targets lost due to SED"
  Ytitle <- "Mean % +/- SE"
  matplot(MovingVar1,Res$SED_loss[,,1], ylim = c(0,100), xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  abline(a=0,b=0, lty=2)
  arrows(MovingVar1, Res$SED_loss[,,1]-Res$SED_loss[,,3], MovingVar1, Res$SED_loss[,,1]+Res$SED_loss[,,3],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "% Fish(tracks) lost due to SED"
  Ytitle <- "Mean % +/- SE"
  matplot(MovingVar1,Res$Track_loss[,,1], ylim = c(0,100), xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  abline(a=0,b=0, lty=2)
  arrows(MovingVar1, Res$Track_loss[,,1]-Res$Track_loss[,,3], MovingVar1, Res$Track_loss[,,1]+Res$Track_loss[,,3],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  # Title <- "%Error Volume density FullPing" 
  # Ytitle <- "Mean %Error +/- SE"
  # matplot(MovingVar1,Res$VolD_MED[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
  #         type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  # legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  # abline(a=0,b=0, lty=2)
  # arrows(MovingVar1, Res$VolD_MED[,,2]-Res$VolD_MED[,,4], MovingVar1, Res$VolD_MED[,,2]+Res$VolD_MED[,,4],
  #        length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  # mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "%Error Volume density only SED FullPing"
  Ytitle <- "Mean %Error +/- SE"
  matplot(MovingVar1,Res$VolD_SED[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  abline(a=0,b=0, lty=2)
  arrows(MovingVar1, Res$VolD_SED[,,2]-Res$VolD_SED[,,4], MovingVar1, Res$VolD_SED[,,2]+Res$VolD_SED[,,4],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "%Error Area density MED Full ping"
  Ytitle <- "Mean %Error +/- SE"
  matplot(MovingVar1,Res$ArD_MED[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  abline(a=0,b=0, lty=2)
  arrows(MovingVar1, Res$ArD_MED[,,2]-Res$ArD_MED[,,4], MovingVar1, Res$ArD_MED[,,2]+Res$ArD_MED_Rm[,,4],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "%Error Area density MED Echocount Yule"
  Ytitle <- "Mean %Error +/- SE"
  matplot(MovingVar1,Res$ArD_MED_Yul[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  abline(a=0,b=0, lty=2)
  arrows(MovingVar1, Res$ArD_MED_Yul[,,2]-Res$ArD_MED_Yul[,,4], MovingVar1, Res$ArD_MED_Yul[,,2]+Res$ArD_MED_Yul[,,4],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  # Title <- "%Error Area density MED at Rm"
  # Ytitle <- "Mean %Error +/- SE"
  # matplot(MovingVar1,Res$ArD_MED_Rm[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
  #         type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  # legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  # abline(a=0,b=0, lty=2)
  # arrows(MovingVar1, Res$ArD_MED_Rm[,,2]-Res$ArD_MED_Rm[,,4], MovingVar1, Res$ArD_MED_Rm[,,2]+Res$ArD_MED_Rm[,,4],
  #        length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  # mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  # 
  # 
  # Title <- "%Error Area density MED at Rm no psi"
  # Ytitle <- "Mean %Error +/- SE"
  # matplot(MovingVar1,Res$ArD_MED_RmOff[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
  #         type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  # legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  # abline(a=0,b=0, lty=2)
  # arrows(MovingVar1, Res$ArD_MED_RmOff[,,2]-Res$ArD_MED_RmOff[,,4], MovingVar1, Res$ArD_MED_RmOff[,,2]+Res$ArD_MED_RmOff[,,4],
  #        length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  # mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "%Error Area density MED Fish Tracks Yule 2000"
  Ytitle <- "Mean %Error +/- SE"
  matplot(MovingVar1,Res$ArD_MED_trkYul[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  abline(a=0,b=0, lty=2)
  arrows(MovingVar1, Res$ArD_MED_trkYul[,,2]-Res$ArD_MED_trkYul[,,4], MovingVar1, Res$ArD_MED_trkYul[,,2]+Res$ArD_MED_trkYul[,,4],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "%Error Area density MED Fish Tracks Volumetric"
  Ytitle <- "Mean %Error +/- SE"
  matplot(MovingVar1,Res$ArD_MED_trk[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  abline(a=0,b=0, lty=2)
  arrows(MovingVar1, Res$ArD_MED_trk[,,2]-Res$ArD_MED_trk[,,4], MovingVar1, Res$ArD_MED_trk[,,2]+Res$ArD_MED_trk[,,4],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  # Title <- "%Error Area density MED Fish Tracks SonarX with psi"
  # Ytitle <- "Mean %Error +/- SE"
  # matplot(MovingVar1,Res$ArD_MED_trkX[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
  #         type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  # legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  # abline(a=0,b=0, lty=2)
  # arrows(MovingVar1, Res$ArD_MED_trkX[,,2]-Res$ArD_MED_trkX[,,4], MovingVar1, Res$ArD_MED_trkX[,,2]+Res$ArD_MED_trkX[,,4],
  #        length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  # mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  # 
  # 
  # Title <- "%Error Area density MED Fish Tracks SonarX no psi"
  # Ytitle <- "Mean %Error +/- SE"
  # matplot(MovingVar1,Res$ArD_MED_trkXOff[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
  #         type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  # legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  # abline(a=0,b=0, lty=2)
  # arrows(MovingVar1, Res$ArD_MED_trkXOff[,,2]-Res$ArD_MED_trkXOff[,,4], MovingVar1, Res$ArD_MED_trkXOff[,,2]+Res$ArD_MED_trkXOff[,,4],
  #        length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  # mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "%Error Area density SED FullPing"
  Ytitle <- "Mean %Error +/- SE"
  matplot(MovingVar1,Res$ArD_SED[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  abline(a=0,b=0, lty=2)
  arrows(MovingVar1, Res$ArD_SED[,,2]-Res$ArD_SED[,,4], MovingVar1, Res$ArD_SED[,,2]+Res$ArD_SED[,,4],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "%Error Area density SED Echocount Yule"
  Ytitle <- "Mean %Error +/- SE"
  matplot(MovingVar1,Res$ArD_SED_Yul[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  abline(a=0,b=0, lty=2)
  arrows(MovingVar1, Res$ArD_SED_Yul[,,2]-Res$ArD_SED_Yul[,,4], MovingVar1, Res$ArD_SED_Yul[,,2]+Res$ArD_SED_Yul[,,4],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  # Title <- "%Error Area density at Rm SED"
  # Ytitle <- "Mean %Error +/- SE"
  # matplot(MovingVar1,Res$ArD_SED_Rm[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
  #         type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  # legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  # abline(a=0,b=0, lty=2)
  # arrows(MovingVar1, Res$ArD_SED_Rm[,,2]-Res$ArD_SED_Rm[,,4], MovingVar1, Res$ArD_SED_Rm[,,2]+Res$ArD_SED_Rm[,,4],
  #        length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  # mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  # 
  # 
  # Title <- "%Error Area density at Rm SED no psi"
  # Ytitle <- "Mean %Error +/- SE"
  # matplot(MovingVar1,Res$ArD_SED_RmOff[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
  #         type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  # legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  # abline(a=0,b=0, lty=2)
  # arrows(MovingVar1, Res$ArD_SED_RmOff[,,2]-Res$ArD_SED_RmOff[,,4], MovingVar1, Res$ArD_SED_RmOff[,,2]+Res$ArD_SED_RmOff[,,4],
  #        length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  # mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "%Error Area density SED Fish Tracks Yule 2000"
  Ytitle <- "Mean %Error +/- SE"
  matplot(MovingVar1,Res$ArD_SED_trkYul[,,2], ylim = c(-100,300), xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  abline(a=0,b=0, lty=2)
  arrows(MovingVar1, Res$ArD_SED_trkYul[,,2]-Res$ArD_SED_trkYul[,,4], MovingVar1, Res$ArD_SED_trkYul[,,2]+Res$ArD_SED_trkYul[,,4],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "%Error Area density SED Fish Tracks min length Yule 2000"
  Ytitle <- "Mean %Error +/- SE"
  matplot(MovingVar1,Res$ArD_SED_L_trkYul[,,2], ylim = c(-100,300), xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  abline(a=0,b=0, lty=2)
  arrows(MovingVar1, Res$ArD_SED_L_trkYul[,,2]-Res$ArD_SED_L_trkYul[,,4], MovingVar1, Res$ArD_SED_L_trkYul[,,2]+Res$ArD_SED_L_trkYul[,,4],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, paste(subtitle,"; MinTrk=",MinTrack,sep=""))
  
  
  
  Title <- "%Error Area density SED Fish Tracks Volumetric"
  Ytitle <- "Mean %Error +/- SE"
  matplot(MovingVar1,Res$ArD_SED_trk[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  abline(a=0,b=0, lty=2)
  arrows(MovingVar1, Res$ArD_SED_trk[,,2]-Res$ArD_SED_trk[,,4], MovingVar1, Res$ArD_SED_trk[,,2]+Res$ArD_SED_trk[,,4],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  # Title <- "%Error Area density SED Fish Tracks SonarX with psi"
  # Ytitle <- "Mean %Error +/- SE"
  # matplot(MovingVar1,Res$ArD_SED_trkX[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
  #         type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  # legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  # abline(a=0,b=0, lty=2)
  # arrows(MovingVar1, Res$ArD_SED_trkX[,,2]-Res$ArD_SED_trkX[,,4], MovingVar1, Res$ArD_SED_trkX[,,2]+Res$ArD_SED_trkX[,,4],
  #        length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  # mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  # 
  # 
  # 
  # Title <- "%Error Area density SED Fish Tracks SonarX no psi"
  # Ytitle <- "Mean %Error +/- SE"
  # matplot(MovingVar1,Res$ArD_SED_trkXOff[,,2], ylim = c(-100,100), xlab=Xtitle , ylab=Ytitle, main = Title,
  #         type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  # legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  # abline(a=0,b=0, lty=2)
  # arrows(MovingVar1, Res$ArD_SED_trkXOff[,,2]-Res$ArD_SED_trkXOff[,,4], MovingVar1, Res$ArD_SED_trkXOff[,,2]+Res$ArD_SED_trkXOff[,,4],
  #        length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  # mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "True Area density"
  Ytitle <- "Density (fish m^-2) +/- SE"
  matplot(MovingVar1,Res$TruArD[,,1], xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  arrows(MovingVar1, Res$TruArD[,,1]-Res$TruArD[,,3], MovingVar1, Res$TruArD[,,1]+Res$TruArD[,,3],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "Area density Estimate FullPing"
  Ytitle <- "Density (fish m^-2) +/- SE"
  matplot(MovingVar1,Res$ArD_MED[,,1], xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  arrows(MovingVar1, Res$ArD_MED[,,1]-Res$ArD_MED[,,3], MovingVar1, Res$ArD_MED[,,1]+Res$ArD_MED[,,3],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "Area density Estimate Yule Echocount"
  Ytitle <- "Density (fish m^-2) +/- SE"
  matplot(MovingVar1,Res$ArD_MED_Yul[,,1], xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  arrows(MovingVar1, Res$ArD_MED_Yul[,,1]-Res$ArD_MED_Yul[,,3], MovingVar1, Res$ArD_MED_Yul[,,1]+Res$ArD_MED_Yul[,,3],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "Area density Estimate Tracks Yule 2000"
  Ytitle <- "Density (fish m^-2) +/- SE"
  matplot(MovingVar1,Res$ArD_MED_trkYul[,,1], xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  arrows(MovingVar1, Res$ArD_MED_trkYul[,,1]-Res$ArD_MED_trkYul[,,3], MovingVar1, Res$ArD_MED_trkYul[,,1]+Res$ArD_MED_trkYul[,,3],
         length=0.05, angle=90, code=3, col=rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "True Volume Density"
  Ytitle <- "Density (fish m^-3) +/- SE"
  matplot(MovingVar1,Res$TruVolD[,,1], xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  arrows(MovingVar1, Res$TruVolD[,,1]-Res$TruVolD[,,3], MovingVar1, Res$TruVolD[,,1]+Res$TruVolD[,,3],
         length=0.05, angle=90, code=3, col= rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  Title <- "Volume density Estimate FullPing"
  Ytitle <- "Density (fish m^-3) +/- SE"
  matplot(MovingVar1,Res$VolD_MED[,,1], xlab=Xtitle , ylab=Ytitle, main = Title,
          type = "b", lty = 1, col= 1:length(MovingVar2),  pch = 16)
  legend('topright',inset=c(-0.36,0),title= legTitle, legend = MovingVar2, col=1:length(MovingVar2), pch=16, cex =0.5, )
  arrows(MovingVar1, Res$VolD_MED[,,1]-Res$VolD_MED[,,3], MovingVar1, Res$VolD_MED[,,1]+Res$VolD_MED[,,3],
         length=0.05, angle=90, code=3, col= rep(1:length(MovingVar2), each=length(MovingVar1)))
  mtext(side=3, line=0, at= min(MovingVar1), adj=-0.05, cex=1, subtitle)
  
  
  end_time <- Sys.time()
  
  # filename <- 'LastRun_CallECS.RData'
  # save.image(file=filename)
  # 
  filename <- '6Depths_6FishD_500.RData'
  setwd ("R:/CLSM/School of Biological Sciences/FisheriesGroup/PhD students/James Dunning/James Code/Echo Counting Simulation")
  # save.image(file=filename)
  #load('5Depths_4FishD_200_FTracks.v2.RData')
  
  cat("End of processing. Cheerio!\n") 
  if (WriteReports){
    writeLines(paste("End of processing reached at ",Sys.time(),'\r',"Cheerio!", sep=""),".End_of_processing_Ella.txt")
  }
}
end_time - start_time
