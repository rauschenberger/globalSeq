## ----eval=FALSE----------------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  BiocInstaller::biocLite("globalSeq")

## ------------------------------------------------------------------------
library(globalSeq)

## ----eval=FALSE----------------------------------------------------------
#  attach(toydata)

## ----echo=FALSE----------------------------------------------------------
names <- names(toydata)
for(i in 1:length(names)){
    assign(names[i], toydata[[i]])
}
rm(names)

## ----eval=FALSE----------------------------------------------------------
#  utils::help(globalSeq)
#  utils::vignette("globalSeq")

## ------------------------------------------------------------------------
cbind(y,X)

## ------------------------------------------------------------------------
set.seed(1)
omnibus(y,X)

## ------------------------------------------------------------------------
rbind(y,offset)

## ------------------------------------------------------------------------
set.seed(1)
omnibus(y,X,offset=offset)

## ------------------------------------------------------------------------
rbind(y,group)

## ------------------------------------------------------------------------
set.seed(1)
omnibus(y,X,group=group)

## ------------------------------------------------------------------------
set.seed(1)
omnibus(y,X,phi=0)

## ------------------------------------------------------------------------
X1 <- X[,c(1:11,15)]
X2 <- X[,12:14]

## ------------------------------------------------------------------------
set.seed(1)
omnibus(y,list(X1,X2))

## ----results='hide'------------------------------------------------------
omnibus(y,X,kind=1) # crude permutation test
omnibus(y,X,kind=0.05) # interrupting permutation
omnibus(y,X,kind=0) # method of control variables

## ----results='hide',fig.width=5,fig.height=2.5,fig.show='hold'-----------
proprius(y,X,type="samples")
proprius(y,X,type="covariates")

## ----results='hide',fig.width=5,fig.height=3,fig.show='hold'-------------
proprius(y,X,type="covariates",alpha=0.05)

## ------------------------------------------------------------------------
cbind(Yloc,Ychr,Y)
cbind(Vloc,Vchr,V)
cbind(Wloc,Wchr,W)

## ----echo=FALSE----------------------------------------------------------
for(i in 1:length(Yloc)){
  cat(paste(names(Yloc[i]),": ",sep=""))
  cat(names(Vloc)[Yloc[i]-5 < Vloc & Vloc < Yloc[i]+5],"\n")
}

## ----results='hide',message=FALSE----------------------------------------
set.seed(1)
cursus(Y,Yloc,V,Vloc,window=5)

## ----results='hide',message=FALSE----------------------------------------
set.seed(1)
cursus(Y,Yloc,V,Vloc,window=5,Ychr,Vchr)

## ----results='hide',message=FALSE----------------------------------------
offset <- colSums(Y) # library sizes
set.seed(1)
cursus(Y,Yloc,V,Vloc,window=5,offset=offset)

## ----message=FALSE-------------------------------------------------------
set.seed(1)
cursus(Y,Yloc,list(V,W),list(Vloc,Wloc),list(5,50))

## ----eval=FALSE----------------------------------------------------------
#  list <- edgeR::DGEList(counts=Y)
#  list <- edgeR::calcNormFactors(list)
#  list <- edgeR::estimateDisp(list)
#  offset <- list$samples$norm.factors
#  phi <- list$tagwise.dispersion
#  
#  cursus(Y,Yloc,V,Vloc,window=5,offset=offset,phi=phi)

