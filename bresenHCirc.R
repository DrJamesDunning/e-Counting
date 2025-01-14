


bresenHCirc <- function(r, Centre, fill = TRUE, plot = TRUE)
{
CX <- vector()
CY <- vector()
MSize <- (Centre-1)*2+1
Circle <- matrix(nrow=MSize, ncol=MSize)
AddCircle <- t(rbind(c(Centre+r,Centre,Centre-r,Centre),c(Centre,Centre+r,Centre,Centre-r)))
Circle[AddCircle] <- 1

CX[1] <- r
CY[1] <- 0
i <- 1

while (CX[i] > CY[i]) {
  RE1 <- abs((CX[i]^2 - 2*CX[i] + 1) + (CY[i]^2 + 2*CY[i] + 1) - r^2)
  RE2 <- abs(CX[i]^2 + (CY[i]^2 + 2*CY[i] + 1) - r^2)
  if (RE1 < RE2) {
    CX[i+1] <- CX[i] - 1
    CY[i+1] <- CY[i] + 1
  }  else {
    CX[i+1] <- CX[i]
    CY[i+1] <- CY[i] + 1
  }
  
  AddCircle <- t(rbind(c(CX,CY,-CY,-CX,-CX,-CY,CY,CX),c(CY,CX,CX,CY,-CY,-CX,-CX,-CY)))
  AddCircle <- AddCircle+Centre
  Circle[AddCircle] <- 1
  
  i <- i+1
}

if (fill) {
  for (i in (Centre-r):(Centre+r)) {
    border <- which(Circle[i,]==1)
    Circle[i,border[1]:max(border)] <- 1
  }

}

if (plot) {image(rotate(Circle),asp=1)}

return(Circle)
}
