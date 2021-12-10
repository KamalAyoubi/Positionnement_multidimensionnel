#Code R : MDS à partir d’une distance euclidienne

library(e1071)

Donnee1 <- read.csv2("FlyingMileage.csv", row.names=1 )

logtsDK_brut <- read.table("logtsDKedit.csv", header = T, sep = "\t",row.names=1)
logtsDK_qualit <- as.matrix(logtsDK_brut[,12:27])
logtsDK_q <- logtsDK_qualit
rownames(logtsDK_q)=logtsDK_brut[,28]

Donnee2 <- hamming.distance(logtsDK_q)

MDS <- function(Donnee){
  D <- as.matrix(Donnee)
  n <- ncol(D)
  
  
  # Matrice des distances  carrée
  D2 <- D*D
  
  # Matrice de centrage : $C = I_{n} - \frac{1}{n} 1_{n}$
  Mat_1 <- matrix(c(rep(1,n^2)),nrow=n)
  Mat_ind <- diag(n)
  
  C <- Mat_ind  - (1/(n)) * Mat_1
  
  # Matrice G des produits scalaires.
  G <- -0.5 * C %*% D2 %*% C
  
  
  # Valeurs propres et vecteurs propres
  Diago <- eigen(G)
  
  
  # Matrice Lambda^1/2 pour 2 dimensions
  Lambda_0.5 <- diag( c(sqrt(Diago$values[1]),sqrt(Diago$values[2])),2)
  
  
  # Vecteurs propres correspondants
  vect_propres <- Diago$vectors[,1:2]
  
  
  # Solution
  Y <- vect_propres %*% Lambda_0.5
  row.names(Y) <- colnames(Donnee)
  
  return(Y)
  
  Y <- MDS(Donnee1)
  YH <- MDS(Donnee1)
  Y
  axe2 <- - Y[,2]
  axe1 <- - Y[,1]
  
  
  plot(axe1,axe2,type="p", xlim=c(-1500,1500) , ylim=c(-700,700),)
  text(axe1,axe2,colnames(Donnee),pos=1,cex=0.8)
}