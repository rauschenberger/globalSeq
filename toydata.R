
# This script simulates data for the vignette.

ls <- ls()
set.seed(pi)

# first set
y <- c(7,1,1,15,6,16,5,2,2,5)
group <- c(1,1,1,2,1,2,1,1,1,1)
offset <- c(15,1,2,28,10,28,9,2,5,7)
X <- matrix(stats::rbinom(n = 150,size = 1,prob = 0.5),nrow = 10,ncol = 15)
X[c(4,6),] <- stats::rbinom(30,size = 30,prob = 0.5)
X[,c(12,13,14)] <- stats::rbinom(n = 30,size = 1,prob = 0.5)
X[2,4:9] <- stats::rbinom(6,size = 30,prob = 0.5)
names(y) <- paste("ind",1:length(y),sep = "")
rownames(X) <- names(y)
colnames(X) <- paste("X",1:ncol(X),sep = "")

# second set
Y <- matrix(stats::rnbinom(5 * 7,mu = 20,size = 1 / 0.5),nrow = 7,ncol = 5)
Yloc <- sort(sample(1:100,size = nrow(Y)))
Ychr <- rep(1,nrow(Y))
V <- matrix(round(stats::rnorm(5 * 10),1),nrow = 10,ncol = 5)
Vloc <- sort(sample(1:100,size = nrow(V)))
Vchr <- rep(1,nrow(V))
W <- matrix(stats::rbinom(n = 5 * 8,size = 2,prob = 0.5),nrow = 8,ncol = 5)
Wloc <- sort(sample(1:100,size = nrow(W)))
Wchr <- rep(1,nrow(W))
colnames(Y) <- colnames(V) <- colnames(W) <- paste("ind",1:5,sep = "")
rownames(Y) <- names(Yloc) <- paste("gene",1:nrow(Y),sep = "")
rownames(V) <- names(Vloc) <- paste("V",1:nrow(V),sep = "")
rownames(W) <- names(Wloc) <- paste("W",1:nrow(W),sep = "")

toydata <- list(
    y = y,X = X,group = group,offset = offset,
    Y = Y,Yloc = Yloc,Ychr = Ychr,V = V,Vloc = Vloc,
    Vchr = Vchr,W = W,Wloc = Wloc,Wchr = Wchr
)

rm(list=setdiff(ls(),c(ls,"toydata")))
