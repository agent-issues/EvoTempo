#Script to estimate prior distributions for relaxed clock models
#first get variance between clock and nonclock tree branch lengths to use for IGR and TK02 model
#then simulate branch lengths under CPP model to get combination of rate multiplier and lambda
#finally check for autocorrelation in evolutionary rates

library(phangorn)
library(phyloch)
library(gplots)

#-------------------------------------------------------------------------------
#functions
#-------------------------------------------------------------------------------
#compares the branch lengths of two trees, including only branches that occur
#in both topologies
compare.branches<-function(tree1,tree2,outgroup=""){
  if (length(tree1$edge)<length(tree2$edge)){
    reftree<-tree1
    testtree<-tree2
  } else {
    reftree<-tree2
    testtree<-tree1
  }
  maxNbBr<-length(reftree$edge)
  BrlenCompMatrix<-matrix(rep(NA,maxNbBr*2),nrow=maxNbBr,ncol=2)
 
  childRef<-reftree$edge[,2]
  childTest<-testtree$edge[,2]
 
  for(i in 1:length(childRef)){
      #don´t use branch if leading to outgroup
      if(!childRef[i]==which(reftree$tip.label==outgroup)){
           refVector<-list.descendents.of(reftree,childRef[i])
           if(!is.na(refVector)) {              
              refNames<-reftree$tip.label[refVector]
              refNames<-sort(refNames)
           
              for(j in 1:length(childTest)){
                testVector<-list.descendents.of(testtree,childTest[j])
                if(!is.na(testVector)){
                  testNames<-testtree$tip.label[testVector]
                  testNames<-sort(testNames)
                  if(all(testNames==refNames)){
                    BrlenCompMatrix[i,1]<-reftree$edge.length[i]
                    BrlenCompMatrix[i,2]<-testtree$edge.length[j]
                    j=length(testtree$edge)
                  }
                }
              }
           }          
      }        
  } 
  #return matrix WITH na´s
  return(BrlenCompMatrix)     
}
#-------------------------------------------------------------------------------
#returns a vector containing the indices of the tips which are descendents of "nodeIndex"
list.descendents.of<-function(tree,nodeIndex){
  if(is.na(nodeIndex)) return(NA)
  if(nodeIndex<1) return(NA)
  
  parent<-tree$edge[,1]
  child<-tree$edge[,2]
  root<-as.integer(parent[!match(parent, child, 0)][1]) #match gives the positions at which the entries of the parent vector occur in the child vector. If =0, then no match and =root node
  tips<-child[is.na(!match(child,parent))]
  
  #if it is a tip
  if(length(which(tips==nodeIndex))>0){
    return(nodeIndex)
  }
  
  is.desc<-c(rep(FALSE,length(tips)))

  for (i in 1:length(tips)){
    index<-tips[i]
    while(index!=root){
      index<-parent[child==index]
      if(index==nodeIndex){
        is.desc[i]=TRUE
        index<-root
      }
    }
  }
  return(tips[is.desc])
}

#-------------------------------------------------------------------------------
#simulating a tree with tk02 branch lengths
simulate_tk02Tree<-function(nbTips=8,treeshape="balanced",rootRate=1.0,tk02Var=0.5,usertree=NULL){

  if(is.null(usertree)) {
    simtree<-stree(nbTips,treeshape)
    simtree$edge.length<-c(rep(1/log2(nbTips),nbTips*2-2))
  } else {
    simtree<-usertree
    nbTips<-length(usertree$tip.label)
    if(is.null(simtree$edge.length)){
        simtree$edge.length<-c(rep(1,length(usertree$edge.length)))
    }
  }

  tk02.node.rate <- c(rep(NA,length(simtree$edge.length)+1))
  tk02.edge.rate <- c(rep(NA,length(simtree$edge.length)))
    
  parent <- simtree$edge[,1]
  child <- simtree$edge[,2]
  root <- as.integer(parent[!match(parent, child, 0)][1]) #match gives the positions at which the entries of the parent vector occur in the child vector. If =0, then no match and =root node

  #to get rates for each node, start with root node -> child node rate sampled from lognorm with mean= ancestral rate and var=brlens*autocorrFactor (="bmvar" in MB) 
  brTodo <- c(rep(1,length(simtree$edge.length)))
  tk02.node.rate[root] <- rootRate
  
  while(sum(brTodo) > 0){
    for (i in 1:length(brTodo)) {

      if(brTodo[i]){  #if node is still to do
        if(!is.na(tk02.node.rate[parent[i]])) {  #if rate already calculated for parent node
        
          meanLinear <-  tk02.node.rate[parent[i]]
          sdLinear <- sqrt(tk02Var*simtree$edge.length[i])
          meanLogNorm <- log(meanLinear) - tk02Var*tk02Var/2
          sdLogNorm <- sqrt(log((sdLinear/exp(2*log(meanLinear)))+1))          
          tk02.node.rate[child[i]] <- rlnorm(1, meanlog = meanLogNorm, sdlog = sdLogNorm)  
          
          brTodo[i] <- 0      
        }
      }
    }
  }

  for (i in 1:length(tk02.edge.rate)) {
    tk02.edge.rate[i] <- mean(c(tk02.node.rate[simtree$edge[i,1]],tk02.node.rate[simtree$edge[i,2]]))
    simtree$edge.length[i] <- simtree$edge.length[i] * tk02.edge.rate[i]
  }

  return(simtree)
}
#-----------------------------------------------------
#simulating a tree with cpp branch lengths
simulate_cppTree<-function(nbTips=8,treeshape="balanced",rootRate=1.0,lambda=1.0,stdevLog=1.0,usertree=NULL){
  if(is.null(usertree)) {
    simtree<-stree(nbTips,treeshape)
    simtree$edge.length<-c(rep(1/log2(nbTips),nbTips*2-2))
  } else {
    simtree<-usertree
    nbTips<-length(usertree$tip.label)
    if(is.null(simtree$edge.length)){
        simtree$edge.length<-c(rep(1,length(usertree$edge.length)))
    }
  }
  
  nbBr<-length(simtree$edge.length)
  nbCPP<-c(rep(0,nbBr))
  waitTimeMatrix<-matrix(nrow=nbBr,ncol=1000)  #not very elegant... dynamic memory alloc in R?
  effBrlens<-c(rep(0.0,nbBr))

  #draw nb of cpp events and waiting times for each branch
  for (i in 1:nbBr){
    restLength<-simtree$edge.length[i]
    while (restLength > 0){
      t<--(1/lambda)*log(runif(1))
      if(t<=restLength) {
        nbCPP[i]<-nbCPP[i]+1
        restLength<-restLength-t
        waitTimeMatrix[i,nbCPP[i]]<-t
      } else {
        waitTimeMatrix[i,nbCPP[i]+1]<-restLength
        restLength<-0
      }
    }
  }

  parent<-simtree$edge[,1]
  child<-simtree$edge[,2]
  root <- as.integer(parent[!match(parent, child, 0)][1]) #match gives the positions at which the entries of the parent vector occur in the child vector. If =0, then no match and =root node

  brRates<-c(rep(NA,nbBr))

  #go through branches, update branchlength if parent branch has been done already
  brTodo<-c(rep(1,nbBr))
  while(sum(brTodo)>0){
    for (i in 1:nbBr) {

      if(brTodo[i]){  #if branch is still to do

        if(parent[i]==root){  #if parent=root, we can begin here
          effBrlens[i]<-0
          for (k in 1:(nbCPP[i]+1)){
            if(k==1) {
              Rate<-rootRate
            } else Rate<-Rate*rlnorm(1,meanlog = 0, sdlog = stdevLog)
            effBrlens[i]<-effBrlens[i]+Rate*waitTimeMatrix[i,k]
          }
          brRates[i]<-Rate #store final rate of branch for its children
          brTodo[i]=0

        } else if(!brTodo[child==parent[i]]){       #or if ancestral has been done
            for (k in 1:(nbCPP[i]+1)){
              if(k==1) {
                Rate<-brRates[child==parent[i]]
              } else Rate<-Rate*rlnorm(1,meanlog = 0, sdlog = stdevLog)
              effBrlens[i]<-effBrlens[i]+Rate*waitTimeMatrix[i,k]
            }
            brRates[i]<-Rate #store final rate of branch for its children
            brTodo[i]=0
        }
      }
    }
  }

  cpp.tree<-simtree
  cpp.tree$edge.length<-effBrlens

  #windows(10,10)
  #plot(cpp.tree)
  
  return(cpp.tree)
}
#-------------------------------------------------------------------------------

#*******************************************************************************
#main commands
#*******************************************************************************

tree1<-read.nexus("nonclock.con.tre")
tree2<-read.nexus("strictclock.con.tre")

BrlensMatrix<-compare.branches(tree1,tree2,outgroup="Orthoptera")
BrlensMatrix2<-na.omit(BrlensMatrix)

#Plotting
windows(15,8)
layout(matrix(1:2,1,2))
hist(BrlensMatrix2[,1],main="non-clock branch lengths")
hist(BrlensMatrix2[,2],main="strict-clock branch lengths")

windows(10,10)
plot(BrlensMatrix2[,2],BrlensMatrix2[,1],main="comparison of clock versus non-clock branch lengths",xlab="strict-clock brlens",ylab="non-clock brlens",xlim=c(0,max(BrlensMatrix2)),ylim=c(0,max(BrlensMatrix2)))
regression<-lm(BrlensMatrix2[,1]~BrlensMatrix2[,2]-1)
abline(coef=c(0,regression$coefficients))
abline(coef=c(0,1),col="red")

#calculate variance between clock and non-clock branches
obsVariance<-c(rep(NA,length(BrlensMatrix2[,2])))
for (i in 1:length(BrlensMatrix2[,2])){
    obsVariance[i]<-BrlensMatrix2[i,1]-BrlensMatrix2[i,2]
    obsVariance[i]<-obsVariance[i]^2
}
windows(10,10)
plot(BrlensMatrix2[,2],obsVariance,main="variance in branchlengths",xlab="t[i]",ylab="(b[i]-t[i])^2")
regression2<-lm(obsVariance~BrlensMatrix2[,2]-1)
abline(coef=c(0,regression2$coefficients))
regression2$coefficients

#take slope as median -> calculate mean (=rate) for exp distr of IBR model
as.numeric(regression2$coefficients)/log(2)

a<-list.descendents.of(tree1,tree1$edge[which(obsVariance==max(obsVariance))+1,2])  #+1 because root branch was not sampled by function ´compare.branches´!
tree1$tip.label[a]

#-------------------------------------------------------------------------------
#simulate different values for the variance parameter in the TK02 model to find 
#an optimal value to use as a prior
expVar <- 0.01867348
nreps <- 100 #nreps per tk02var parameter

tk02varVector <- seq(from=0.005, to=0.05, by=0.005)
nVar <- length(tk02varVector)

Outputmatrix <- matrix(nrow=nreps,ncol=nVar)

for (i in 1:nVar){
    for (k in 1:nreps){
        #simulate tree with tk02 branchlengths    
        tk02tree <- simulate_tk02Tree(rootRate=1.0,tk02Var=tk02varVector[i],usertree=tree2)
      
        #compare branches of bm tree with strict clock tree
        BrlensMatrix <- compare.branches(tk02tree,tree2,outgroup="Orthoptera")
        BrlensMatrix <- na.omit(BrlensMatrix)

        obsVariance <- c(rep(NA,length(BrlensMatrix[,2])))
        for (m in 1:length(obsVariance)){
            obsVariance[m] <- BrlensMatrix[m,1]-BrlensMatrix[m,2]
            obsVariance[m] <- obsVariance[m]^2
        }
        regr <- lm(obsVariance~BrlensMatrix[,2]-1)
        Outputmatrix[k,i] <- as.numeric(regr$coefficients)
         #reset matrix
        for(n in 1:length(BrlensMatrix[,2])){
            BrlensMatrix[n,1]<-NA
            BrlensMatrix[n,2]<-NA
        }
    }
}
dimnames(Outputmatrix)<-list(c(1:nreps),round(tk02varVector, digits=6))

boxplot(Outputmatrix, style="quantile", las=3,xlab="tk02var values", ylab="obtained variance in brlens")
abline(h=expVar, col="red")


#-------------------------------------------------------------------------------
#simulate different combinations of cpp brlens to find ideal lambda and ratemultiplier (cppratepr and cppmultdevpr)

stricttree<-read.nexus("strictclock.con.tre")
expVar<-0.01867348

nreps<-100 #nreps per lambda/rmult combination
LambdaVector<-seq(from=0.25, to=2.5, by=0.25)
nbLambda<-length(LambdaVector)
minRmult<-(-1.5) 
maxRmult<-0.75     
stepsizeRmult<-0.25
RmultVector<-exp(minRmult+(((1:10)-1)*stepsizeRmult))
nbRmult <- length(RmultVector)
Outputmatrix<-matrix(nrow=nbLambda,ncol=nbRmult)

for (i in 1:nbLambda){
  for (j in 1:nbRmult){
    Lvalue<-LambdaVector[i]
    sdLog<-RmultVector[j]
    BrlensVar<-c(rep(NA,nreps))
    
    for (k in 1:nreps){
      cpptree<-simulate_cppTree(nbTips=8,treeshape="balanced",rootRate=1.0,lambda=Lvalue,stdevLog=sdLog,usertree=stricttree)

      if(max(cpptree$edge.length<1000)){

        #new way: compare branches of cpp with strict clock tree
        BrlensMatrix<-compare.branches(cpptree,stricttree,outgroup="Orthoptera")
        BrlensMatrix<-na.omit(BrlensMatrix)

        obsVariance<-c(rep(NA,length(BrlensMatrix[,2])))
        for (m in 1:length(obsVariance)){
            obsVariance[m]<-BrlensMatrix[m,1]-BrlensMatrix[m,2]
            obsVariance[m]<-obsVariance[m]^2
        }
        regr<-lm(obsVariance~BrlensMatrix[,2]-1)
        BrlensVar[k]<-as.numeric(regr$coefficients)
        write.table(BrlensVar[k],file="LowNcppSimuls_Values.txt",append=TRUE,row.names=F,col.names=F)

        #reset matrix
        for(n in 1:length(BrlensMatrix[,2])){
            BrlensMatrix[n,1]<-NA
            BrlensMatrix[n,2]<-NA
        }
      } else {
        BrlensVar[k]<-NA
      }
    }
    
    Outputmatrix[i,j]<-median(BrlensVar)-expVar
    if (Outputmatrix[i,j]<0) Outputmatrix[i,j]<-(-Outputmatrix[i,j])
  }
  write.table(Outputmatrix[i,],file="Simuls_Matrix1.txt",append=TRUE,row.names=F)

}

#plotting from file after saving (make new files from saved values in Simuls_Matrix1.txt in table format including row and column labels
Outputmatrix2<-read.table(file="Lambda-5-55_Rmult--2-5-0-5.txt",header=T)
attach(Outputmatrix2)
for (i in 1:dim(Outputmatrix2)[1]){
  for (j in 1:dim(Outputmatrix2)[2]){
    if (Outputmatrix2[i,j]>1.0) Outputmatrix2[i,j]=1.0
  }   
}
#Outputmatrix2<-log(Outputmatrix2)
Outputmatrix2<-Outputmatrix2[1:11,]
nbBreaks<-50
myBreaks<-c(0,(quantile(as.matrix(Outputmatrix2),probs=c((1:nbBreaks)/nbBreaks),na.rm=T)))
image(1:dim(Outputmatrix2)[1], 1:dim(Outputmatrix2)[2], as.matrix(Outputmatrix2), col = rainbow(nbBreaks,end=(4/6)),breaks=myBreaks,col.axis="white",xlab="rMult", ylab="Lambda",main="difference to targeted variance")
axis(1,at=c(1:dim(Outputmatrix2)[1]),labels=rownames(Outputmatrix2),cex.axis=0.8,las=2)
axis(2,at=c(1:dim(Outputmatrix2)[2]),labels=names(Outputmatrix2),cex.axis=0.8,las=2)



#-------------------------------------------------------------------------------
#check for autocorrelation of rates: get log(b[i]/t[i] / b[anc]/t[anc]) for each branch
AutocorVec<-c(rep(NA,length(tree1$edge[,1])))

for(i in 1:(length(tree1$edge[,1]))){   #ohne root
    parentBranch<-which(tree1$edge[,2]==tree1$edge[i,1])
    if(length(parentBranch)==1){#if there is a parent branch 
         AutocorVec[i]<-(BrlensMatrix[i,1]/BrlensMatrix[i,2])/(BrlensMatrix[parentBranch,1]/BrlensMatrix[parentBranch,2]) #should be 1 for complete autocorrelation
         if(log(AutocorVec[i])>0) {
            AutocorVec[i]<-log(AutocorVec[i])
         } else AutocorVec[i]<-(-log(AutocorVec[i]))
    }
}
AutocorVec<-na.omit(AutocorVec)

#and the same for random pairs
nreps <- 1000
HistBreaks <- c(0.25*(0:20))
HistMatrix <- matrix(nrow = nreps, ncol = length(HistBreaks)-1)

for(a in 1:nreps){
    UncorVec <- c(rep(NA,length(AutocorVec)))

    for(i in 1:length(AutocorVec)){
      branch1<-ceiling(runif(1)*length(tree1$edge[,1]))
      branch2<-branch1
      while(branch2==branch1){
        branch2<-ceiling(runif(1)*length(tree1$edge[,1]))
        #exclude parent-descendent pairs
        if(tree1$edge[branch1,1]==tree1$edge[branch2,2]){
          branch2<-branch1
        }
        if(tree1$edge[branch1,2]==tree1$edge[branch2,1]){
          branch2<-branch1
        }
      }
      UncorVec[i]<-(BrlensMatrix[branch1,1]/BrlensMatrix[branch1,2])/(BrlensMatrix[branch2,1]/BrlensMatrix[branch2,2])
      if(!is.na(UncorVec[i])){
        if((log(UncorVec[i]))>0) {
          UncorVec[i]<-log(UncorVec[i])
        } else UncorVec[i]<-(-log(UncorVec[i]))
      }
    }
    histo <- hist(UncorVec,breaks=HistBreaks,plot=FALSE) 
    HistMatrix[a,] <-histo$counts   
}
FakeVec <- NA
SDVec <- NA
for (j in 1:(length(HistBreaks) -1)) {
  FakeVec <- append(FakeVec, round(mean(HistMatrix[,j])))
  SDVec <- append(SDVec, sd(HistMatrix[,j])) 
}
FakeVec <- na.omit(FakeVec)
SDVec <- na.omit(SDVec)

histAuto <- hist(AutocorVec, breaks=HistBreaks, plot=F)
histAuto <- histAuto$counts
histFake <- as.vector(FakeVec)
histData <- rep(NA, 2*length(histAuto))
histData[seq(from=1,to=length(histData),by=2)] <- histAuto
histData[seq(from=2,to=length(histData),by=2)] <- histFake
xLabels <- rep(NA, length(HistBreaks)*2-2)
xLabels[seq(from=1,to=length(xLabels),by=2)] <- HistBreaks[1:length(HistBreaks)-1]
barplot(histData, col=c("blue", "white"), names.arg=xLabels, axisnames=T, width=1.6, space=0.25, main="distribution of branch rates: random and parent-offspring pairs",xlab="log(( b[i]/t[i])/(b[j]/t[j]) )")

xers <- seq(from=3.2, to=80, by=4)
y1ers <- histFake - as.vector(SDVec)
for (l in 1:length(y1ers)) {
    if(y1ers[l] < 0) 
      y1ers[l] <- 0
}
y2ers <- histFake + as.vector(SDVec)
arrows(xers, y1ers, xers, y2ers, length=0.05, angle=90, code=3)

#-------------------------------------------------------------------------------
