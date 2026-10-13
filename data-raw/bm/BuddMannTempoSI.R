#Supplementary code for Budd & Mann (2024) A covariant model of molecular and diversification patterns through time and the history of large clades

library(deSolve)
library(ReacTran)
library(parallel)

#NB: By far the most computationally intensive part of this script is calculating Q: the diversity probability distribution. 
#As such, an option is included below to skip this calculation if running the code multiple times [assuming Q is already calculated and in the environment]

#Set calcQ = F if Q matrix has already been calculated for multiple runs, as this is the slowest element (days)
calcQ=F


#Set up global parameters for model and PDE solvers
lx=101 #number of points on the log-tempo grid

s=1
theta = 0.01 #parameters of the diversification 
r=0.01
mu=0.5
lambda = mu+r

#Set up a grid of log-tempos for the PDE solver
Grid=setup.grid.1D(x.up=-10, x.down=10, N=lx)
x = Grid$x.mid
ex = setup.prop.1D(grid=Grid, func=function(x) exp(x))
x0idx = which.min(abs(x))

##Define required PDE functions - refer to documentation of ReacTran package for details of notation.

pdeG <- function(t, C, parms){ #PDE for the diversity generating function
  A=1
  tran1 <- tran.1D(C=C, A=A, D=1, v=0, dx = Grid)
  tran2 <- tran.1D(C=C, A=A, D=0, v=1, dx = Grid)
  
  dC = (tran1$dC*parms[3]*parms[4]^2 + tran2$dC*parms[3]*Grid$x.mid + (parms[1]*C-parms[2])*(C-1))*ex$mid
  
  list(dC=dC)
}

pdeU <- function(t, C, parms){ #PDE for for U=G-mu/lambda
  A=1
  tran1 <- tran.1D(C=C, A=A, D=1, v=0, dx = Grid)
  tran2 <- tran.1D(C=C, A=A, D=0, v=1, dx = Grid)
  
  dC = (tran1$dC*parms[3]*parms[4]^2 + tran2$dC*parms[3]*Grid$x.mid + parms[1]*C*(C-1+parms[2]/parms[1]))*ex$mid
  
  list(dC=dC)
}

pdeM <- function(t, C, parms){ #PDE for mean diversity (unconditioned)
  A=1
  tran1 <- tran.1D(C=C, A=A, D=1, v=0, dx = Grid)
  tran2 <- tran.1D(C=C, A=A, D=0, v=1, dx = Grid)
  
  dC = (tran1$dC*parms[2]*parms[3]^2 + tran2$dC*parms[2]*Grid$x.mid + parms[1]*C)*ex$mid
  list(dC=dC)
}

pdeMM2 <- function(t, CC2, parms){ #joint PDEs for mean diversity and mean-square diversity (unconditioned)
  A=1
  lambda = parms[1]
  mu=parms[2]
  theta=parms[3]
  s=parms[4]
  r=lambda-mu
  C=CC2[1:lx]
  C2=CC2[lx+(1:lx)]
  tran1 <- tran.1D(C=C, A=A, D=1, v=0, dx = Grid)
  tran2 <- tran.1D(C=C, A=A, D=0, v=1, dx = Grid)
  
  dC = (tran1$dC*theta*s^2 + tran2$dC*theta*Grid$x.mid + r*C)*ex$mid
  
  tran1 <- tran.1D(C=C2, A=A, D=1, v=0, dx = Grid)
  tran2 <- tran.1D(C=C2, A=A, D=0, v=1, dx = Grid)
  
  dC2 = (tran1$dC*theta*s^2 + tran2$dC*theta*Grid$x.mid + r*C2 +2*lambda*C^2)*ex$mid
  dC2 = dC2+dC
  dCC2 = c(dC, dC2)
  list(dCC2=dCC2)
}

pdeC <- function(t, CC2, parms){ #PDE for correction term C for conditioning on survival
  A=1
  lambda = parms[1]
  mu=parms[2]
  theta=parms[3]
  s=parms[4]
  C=CC2[1:lx]
  C2=CC2[lx+(1:lx)]
  tran1 <- tran.1D(C=C, A=A, D=1, v=0, dx = Grid)
  tran2 <- tran.1D(C=C, A=A, D=0, v=1, dx = Grid)
  dC = (tran1$dC*theta*s^2 + tran2$dC*theta*Grid$x.mid + ((2*C2-1)*lambda+mu)*C)*ex$mid
  tran1 <- tran.1D(C=C2, A=A, D=1, v=0, dx = Grid)
  tran2 <- tran.1D(C=C2, A=A, D=0, v=1, dx = Grid)
  dC2 =  (tran1$dC*theta*s^2 + tran2$dC*theta*Grid$x.mid + lambda*C2*(C2-1))*ex$mid
  
  list(dC=c(dC, dC2))
}

evalQ <- function(MM=2^21){
  #Function to evaluate the diversity probability distribution. MM specifies the number of evaluations and thus indirectly to what number of species to distribution is calculated. 
  #with r= theta = 0.01, s = 1, MM needs to be at least 2^20 to reach high enough number of species (higher MM gives more precise answers)
  
  mm = 0:(MM-1)
  z = matrix(NA, lx,MM)
  idxs = splitIndices(MM, MM/(2^10)) #Split up the jobs into 2^10 part chunks to avoid maxing out the RAM
  #If you have more than 8GB of RAM you could increase the chunk sizes but it doesn't make much difference
  for (j in 1:length(idxs)){
    print(j)
    g <- function(z0){zvode(rep(z0-mu/lambda, lx), c(0, 500), pdeU, parms=c(lambda, mu, theta, s), atol=1e-18)[2, 2:(lx+1)] + mu/lambda}
    #z = sapply(exp(2*pi*1i*mm/MM), g) #Use this line instead of the one below to remove parallel processing
    zz = mcmapply(g, exp(2*pi*1i*mm[idxs[[j]]]/MM), mc.cores=6, mc.preschedule = T) #Using 6 cores is optimal on my computer for avoiding RAM problems - can be increased for more powerful machines 
    z[, idxs[[j]]] = zz
    gc() #Garbage collection to avoid RAM issues
  }
  Q = matrix(NA, dim(z)[1], dim(z)[2])
  for (i in 1:dim(z)[1]){
    Q[i,] = abs(fft(z[i, ]))/MM
  }
return(Q)
}

evalM <- function(r, theta, s){
  #Calculate matrix M_t,y,x - expected diversity of species with log-tempo y at time t, starting from log-tempo x
  M = array(NA, c(501, lx, lx))
  for (i in 1:lx){
    m0 = rep(0, lx)
    m0[i] = 1 #Set up correct initia condition
    tmp  = vode(m0,seq(0, 500, 1), pdeM,  parms=c(r, theta, s), atol = 1e-15)
    M[, i, ] = tmp[, 2:(lx+1)]
  }
  return(M)
}

##----------------------------------------Calculate some fundamental quantities that are used across analyses


if (calcQ){Q = evalQ(MM=2^21)} #This line takes a long time, so it is a good idea to save Q and set calcQ=F for multiple runs with the same params

M = evalM(r, theta, s) #M matrix with standard parameters
Mhighr = evalM(2*r, theta, s) #M matrix with higher r 
Mlowr = evalM(0, theta, s) #M matrix with r=0

#Calculate extinction (D) and survival (S) probabilities for different starting tempos over time
D = vode(rep(0, lx), seq(0, 500, 1), pdeG, parms=c(lambda, mu, theta, s), atol=1e-15)[, 2:(lx+1)]
S=1-D

#Calculate the mean diversity and mean-square diversity over time for different starting tempos
meanNfromx = vode(rep(1, lx),seq(0, 500, 1), pdeM,  parms=c(r, theta, s), atol = 1e-15)[, 2:(lx+1)]
meanNN2fromx = vode(rep(1, 2*lx),seq(0, 500, 1), pdeMM2,  parms=c(lambda, mu, theta, s), atol = 1e-15)[, 2:(2*lx+1)]

#Get mean clade size and mean experienced clade size among clades that survive 500 MYrs
meanclade = meanNN2fromx[501, x0idx]/(1-D[501, x0idx])
meanexclade = (meanNN2fromx[501, x0idx+lx]/(1-D[501, x0idx]))/meanclade


##----------------------------------------Figure 1 analysis


#Calculate probability that one sister group takes a given proportion of total diversity in crown group of 500Myrs
#in increments of dw between 0 and 1.

#Define increments
dw = 0.05
ws = seq(0, 0.5-dw, dw)
#Get prob. distribution of surviving clade sizes starting at log-tempo x=0 
p = Q[x0idx, 2:2^21]
p = p/sum(p)
P = 1-cumsum(p)
n = 1:length(p)

Pw = rep(NA, length(ws))
Pw[1] = 0
for (i in 2:length(ws)){
  w = ws[i]
  n2c = ceiling((1-w)*n/w)
  Pw[i] = sum(p*P[n2c], na.rm=T)
}

ws = c(ws, 0.5, rev(1-ws))
Pw = c(Pw, 0.5, 1-rev(Pw))

#Make figure 1
png("~/Dropbox/R/tempo/tempofigs/tempofig1.png", width = 18, height=9, units="cm", res=300)
layout(matrix(1:2, 1,2))
par(mar = c(4.1,4.1, 1.3,1.3))
plot(1:round(2*meanexclade), Q[x0idx, 2:(1+round(2*meanexclade))], type='l', bty='l', log='y', xlab='Species', ylab='Probability')
abline(v=meanclade, col='blue')
abline(v=meanexclade,col='red')
p_geom = 1/meanclade
x_geom= seq(1, 2*meanexclade, 10000)
Q_geom = dgeom(x_geom, p_geom)
lines(x_geom, Q_geom, lty=2)
text(meanclade+50, 1e-5, 'Mean\nclade', pos=4, font=3, cex=0.7, col='blue')
text(meanexclade+50, 1e-5, 'Mean experienced\nclade', pos=4, font=3, cex=0.7, col='red')
text(meanclade*2, Q[x0idx, 2], 'A', font=2, cex=1.5)
b=barplot(diff(Pw), col='white', border='black', space=0, ylab='Probability', xlab='Proportion in one sister group')
lines(seq(0,20, 1),rep(1/20, 21), lty=2)
axis(1, at = seq(min(b)-mean(diff(b)/2), max(b)+mean(diff(b)/2), length.out=5), labels=signif(seq(0,1, length.out=5), 2))
text(b[3], max(diff(Pw))*0.96, 'B', font=2, cex=1.5)
dev.off()
layout(1)

##----------------------------------------Figure 2 analysis

#Apply conditioning on survival for expected diversification plot
mean_species_surviving = meanNfromx/S
mean_lineages=matrix(NA, 501, lx)
for (i in 1:lx){
m = M[, , i]
Sp = 1-D[rev(1:501), ]
mean_lineages[, i] = rowSums(m*Sp)/S[501, i]
}

meanNfromx_correction=matrix(NA, 501, lx)
meanNfromx_conditional = matrix(NA, 501, lx)
for (i in 1:501){
m0 = c(D[501-(i-1), ], D[501-(i-1), ]-mu/lambda)
meanNfromx_correction[i, ]  = vode(m0,c(500-(i-1), 500), pdeC,  parms=c(lambda, mu, theta, s), atol = 1e-15)[2, 2:(lx+1)]
meanNfromx_conditional[i,] = (meanNfromx[i, ]-meanNfromx_correction[i,])/(1-D[501, ])
}

#Make figure 2
png("~/Dropbox/R/tempo/tempofigs/tempofig2.png", width = 18, height=9, units="cm", res=300)

layout(matrix(1:2, 1,2))
par(mar = c(4.1,4.1, 1.3,1.3))
plot(0:500, 0:500, type='n', xaxt='n', xlab = 'Time (Ma)', ylab ='Species', bty='l', log='y', ylim=c(1, max(meanNfromx_conditional)), main='All clades')
axis(1, at =seq(0, 500, 100), labels=seq(500, 0, -100))
idx = which.min(abs(Grid$x.mid))
lines(0:500, meanNfromx[, idx], col='black', lwd=2)
idx = which.min(abs(Grid$x.mid+2))
lines(0:500, meanNfromx[, idx], col='blue', lwd=2)
idx = which.min(abs(Grid$x.mid-2))
lines(0:500, meanNfromx[, idx], col='red', lwd=2)
legend('topright', lty=1, col=c('red', 'black', 'blue'), legend = c(expression(x[0] == 2), 
                                                                    expression(x[0] == 0),
                                                                    expression(x[0] == -2)), cex=0.75)
text(10, max(meanNfromx_conditional), 'A', font=2, cex=1.5)
plot(0:500, 0:500, type='n', xaxt='n', xlab = 'Time (Ma)', ylab ='Species', bty='l', log='y', ylim=c(1, max(meanNfromx_conditional)), main='Surviving clades')
axis(1, at =seq(0, 500, 100), labels=seq(500, 0, -100))
idx = which.min(abs(Grid$x.mid))
lines(0:500, meanNfromx_conditional[, idx], col='black', lwd=2)
lines(0:500, mean_lineages[, idx], col='black', lwd=2, lty=2)
idx = which.min(abs(Grid$x.mid+2))
lines(0:500, meanNfromx_conditional[, idx], col='blue', lwd=2)
lines(0:500, mean_lineages[, idx], col='blue', lwd=2, lty=2)
idx = which.min(abs(Grid$x.mid-2))
lines(0:500, meanNfromx_conditional[, idx], col='red', lwd=2)
lines(0:500, mean_lineages[, idx], col='red', lwd=2, lty=2)
text(10, max(meanNfromx_conditional), 'B', font=2, cex=1.5)
dev.off()
layout(1)



##----------------------------------------Figure 3 analysis

evalmx <- function(M=M, x0=c(0,-2, 2)){
#From the M matrix, evaluate the mean and sd log-tempo over time for a process starting at log-tempos x0
mx = matrix(NA, 501, 3)
mx2 = mx
for (it in 1:501){
  idx = which.min(abs(Grid$x.mid)-x0[1])
  p = M[it, , idx]
  p= p/sum(p)
  mx[it, 1] = sum(p*Grid$x.mid)
  mx2[it, 1] = sum(p*Grid$x.mid^2)
  
  idx = which.min(abs(Grid$x.mid-x0[2]))
  p = M[it, , idx]
  p= p/sum(p)
  mx[it, 2] = sum(p*Grid$x.mid)
  mx2[it, 2] = sum(p*Grid$x.mid^2)
  
  idx = which.min(abs(Grid$x.mid-x0[3]))
  p = M[it, , idx]
  p= p/sum(p)
  mx[it, 3] = sum(p*Grid$x.mid)
  mx2[it, 3] = sum(p*Grid$x.mid^2)
}
sdx = sqrt(mx2-mx^2)
return(list(mx = mx, sdx=sdx))
}

#mx matrix for processes with standard parameters, starting from x = 0,-2 and 2
tmp = evalmx(M=M, x0=c(0,-2, 2))
mx = tmp$mx
sdx = tmp$sdx

#mx matrix for processes with high r, starting from x = 0,-2 and 2
tmp = evalmx(M=Mhighr, x0=c(0,-2, 2))
mxhr = tmp$mx
sdxhr = tmp$sdx

#mx matrix for processes with r=0, starting from x = 0,-2 and 2
tmp = evalmx(M=Mlowr, x0=c(0,-2, 2))
mxlr = tmp$mx
sdxlr = tmp$sdx

#Make figure 3
png("~/Dropbox/R/tempo/tempofigs/tempofig3.png", width = 18, height=9, units="cm", res=300)

layout(matrix(1:2, 1,2))
par(mar = c(4.1,4.1, 1.3,1.3))
cols = c('black', 'blue', 'red')
plot(0:500, mx[, 1], ylim = c(-2.5, 2.5), bty='l', xlab='Time (Ma)', ylab='Log tempo (x)', type='n', xaxt='n')
axis(1, at =seq(0, 500, 100), labels=seq(500, 0, -100))
for (i in 1:3){
  polygon(c(0:500, rev(0:500)), c(mx[,i]-sdx[, i], rev(mx[, i]+sdx[, i])), border=NA, col=scales::alpha(cols[i], 0.2))
}
for (i in 1:3){
  lines(0:500, mx[, i], col=cols[i])
}
legend('bottomright', lty=1, col=c('red', 'black', 'blue'), legend = c(expression(x[0] == 2), 
                                                                       expression(x[0] == 0),
                                                                       expression(x[0] == -2)), cex=0.75)
text(5, 2.5, 'A', font=2, cex=1.5)

cols = c('black', 'darkgreen', 'purple')
plot(0:500, mx[, 1], ylim = c(-2.5, 2.5), bty='l', xlab='Time (Ma)', ylab='Log tempo (x)', type='n', xaxt='n')
axis(1, at =seq(0, 500, 100), labels=seq(500, 0, -100))
polygon(c(0:500, rev(0:500)), c(mx[,1]-sdx[, 1], rev(mx[, 1]+sdx[, 1])), border=NA, col=scales::alpha(cols[1], 0.2))
polygon(c(0:500, rev(0:500)), c(mxlr[,1]-sdxlr[, 1], rev(mxlr[, 1]+sdxlr[, 1])), border=NA, col=scales::alpha(cols[2], 0.2))
polygon(c(0:500, rev(0:500)), c(mxhr[,1]-sdxhr[, 1], rev(mxhr[, 1]+sdxhr[, 1])), border=NA, col=scales::alpha(cols[3], 0.2))

lines(0:500, mx[, 1], col=cols[1])
lines(0:500, mxlr[, 1], col=cols[2])
lines(0:500, mxhr[, 1], col=cols[3])
legend('bottomleft', legend = c(expression(paste('r = 2', theta)), expression(paste('r = ', theta)), 'r = 0'), col=cols[c(3,1,2)], lty=1, cex=0.75)
text(5, 2.5, 'B', font=2, cex=1.5)
dev.off()


##----------------------------------------Figure 4 analysis

#Calculate the long term log-tempo distribution as the final distribution from the highest starting tempo
Px = M[501, ,lx]
Px=Px/sum(Px) #Normalize

#Taking this stable distribution as the distribution from which the starting tempo might be drawm, caluclate marginal distributions and means:

#Marginal mean and mean-square diversity over time, integrated over possible starting tempos
margmean =sum(meanNfromx[501,]*Px/(1-Q[, 1]))
margexmean = (sum(meanNN2fromx[501,102:202]*Px/(1-Q[, 1])))/margmean

#Marginal diversity distribution for surviving clades, integrated over possible starting tempos
Qmarg=rep(0, dim(Q)[2])
for (i in 1:lx)
{
  Qmarg = Qmarg+Px[i]*Q[i, ]
}
Qmargp = Qmarg[2:dim(Q)[2]]
Qmargp = Qmargp/sum(Qmargp)

cladesizes = round(c(margexmean))
ExtnT=rep(NA, 501)
ExsqtnT = ExtnT
cladesize = round(c(margexmean))

#Calculate the expected historical tempo of members of a mean experienced sized clade
PnTx0 = Q[, cladesize+1]
PnT = sum(PnTx0*Px)
PxtnT= matrix(NA, 501, lx)
for (i in 1:501){
  Pxtx0 = M[i, ,]
    for (j in 1:lx){
      Pxtx0[, j] = Pxtx0[, j]/sum(Pxtx0[, j])
    }
    for (j in 1:lx){
      PxtnT[i, j] = sum(Pxtx0[j, ]*PnTx0*Px)/PnT
    }
  
    ExtnT[i] = sum(x*PxtnT[i, ])
    ExsqtnT[i] = sum(x^2*PxtnT[i, ])
  }

sdxtnT = sqrt(ExsqtnT-ExtnT^2)
UxtnT = ExtnT+sdxtnT
LxtnT = ExtnT-sdxtnT

#Calculate the historical tempo of ancestors of a randomly chosen modern species
Px2 = matrix(NA, 501, lx)
Ex = rep(NA, 501)
Lx = Ex
Ux = Ex
Exsq = Ex
for (i in 1:501){
  Px2[i,] = Px*meanNfromx[i, ]
  Px2[i, ] =Px2[i, ]/sum(Px2[i, ])
  Ex[i] = sum(Px2[i,]*Grid$x.mid)
  Exsq[i] = sum(Px2[i,]*Grid$x.mid^2)
}
sdx=sqrt(Exsq-Ex^2)
Ex = rev(Ex)
sdx = rev(sdx)
Ux = rev(Ux)
Lx = rev(Lx)
Ux = Ex+sdx
Lx = Ex-sdx

#Make figure 4
png("~/Dropbox/R/tempo/tempofigs/tempofig4.png", width = 18, height=9, units="cm", res=300)
par(mar = c(4.1,4.1, 1.3,1.3))
layout(matrix(1:2, 1,2))
rangemin = min(c(min(Lx), min(LxtnT)))
rangemax = max(c(max(Ux), max(UxtnT)))
plot(0:500, Ex, type='n', bty='l', xlab='Time (Ma)', ylab='Log tempo (x)', ylim=c(rangemin, rangemax), yaxt='n', xaxt='n')
axis(1, at=seq(0, 500, 100), labels=seq(500, 0, -100))
axis(2, at=seq(-10, 10, 0.5), labels = signif((seq(-10, 10, 0.5)), 2))
polygon(c(0:500, rev(0:500)), c(Ux, rev(Lx)), border = NA, col=scales::alpha('black', 0.2))

lines(0:500, Ex,lwd=2)
text(30, rangemax, 'A', font=2, cex=1.5)

plot(0:500, ExtnT, type='n', bty='l', xlab='Time (Ma)', ylab='Log tempo (x)', xaxt='n', ylim=c(rangemin, rangemax))
axis(1, at=seq(0, 500, 100), labels = rev(seq(0, 500, 100)))

polygon(c(0:500, rev(0:500)), c(UxtnT, rev(LxtnT)), border = NA, col=scales::alpha('black', 0.2))
lines(0:500, ExtnT, col='black', lwd=2)

text(30, rangemax, 'B', font=2, cex=1.5)
dev.off()
layout(1)


##----------------------------------------Figure 5 analysis

ts = seq(0, 20, 0.01) #times to evaluate branch length probability at

tmp=vode(ex$mid*(lambda+mu),seq(0, 20, 0.01), pdeM,  parms=c(-(lambda+mu), theta, s), atol = 1e-15)[, 2:(lx+1)]

ys1 = matrix(NA, 3, length(ts))
ys2 = ys1
idx = which.min(abs(Grid$x.mid+2))
ys1[1, ]=  tmp[, idx]
ys1[1, ] = 1.3*ys1[1, ]/max(ys1[1, ])
ys2[1, ]=  dexp(ts, lambda+mu)
ys2[1, ] = 1.3*ys2[1, ]/max(ys2[1, ])
idx = which.min(abs(Grid$x.mid))
ys1[2, ]= tmp[, idx]
ys1[2, ] = 1.3*ys1[2, ]/max(ys1[2, ])
ys2[2, ]= dexp(ts, lambda+mu)
ys2[2, ] = 1.3*ys2[2, ]/max(ys2[2, ])
idx = which.min(abs(Grid$x.mid-2))
ys1[3, ]= tmp[, idx]
ys1[3, ] = 1.3*ys1[3, ]/max(ys1[3, ])
ys2[3, ]= dexp(ts, lambda+mu)
ys2[3, ] = 1.3*ys2[3, ]/max(ys2[3, ])

#Make figure 5
png("~/Dropbox/R/tempo/tempofigs/tempofig5.png", width = 18, height=9, units="cm", res=300)
par(mar = c(4.1,3, 1.3,1.3))
layout(matrix(1:2, 1, 2))

cols = c('blue', 'black', 'red')
plot.new()
plot.window(xlim = c(0,15), ylim = c(0.9,  4.6), bty='l', xaxs='i')
axis(1)
axis(2, at=c(0, 10))
title(xlab='Branch duration (MYrs)', line =2.2)
title(ylab='Density', line =1)
for (j in 1:3){
  polygon(c(ts, rev(ts)), c(rep(4-j, length(ts)), 4-j+rev(ys1[j,])), border = cols[j], col=scales::alpha(cols[j], 0.2))
  abline(h=4-j, col = cols[j])
}
text(1.5, 4.5, 'A', font=2, cex=1.5)
plot.new()
plot.window(xlim = c(0,15), ylim = c(0.9,  4.6), bty='l', xaxs='i')
axis(1)
axis(2, at=c(0, 10))
title(xlab='Molecular change (MYrs-equiv.)', line =2.2)
title(ylab='Density', line =1)
for (j in 1:3){
  polygon(c(ts, rev(ts)), c(rep(4-j, length(ts)), 4-j+rev(ys2[j,])), border = cols[j], col=scales::alpha(cols[j], 0.2))
  abline(h=4-j, col = cols[j])
}
text(1.5, 4.5, 'B', font=2, cex=1.5)
dev.off()



