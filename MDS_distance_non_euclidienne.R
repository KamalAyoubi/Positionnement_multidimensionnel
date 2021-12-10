#Code R : MDS à partir d’une distance non euclidienne


diagonalisation  <- function(Donnee){
  D <- as.matrix(Donnee )
  n <- ncol(D)
  
  # Matrice des distances au carré
  D2 <- D*D
  
  # Matrice de centrage
  Mat_1 <- matrix(c(rep(1,n^2)),nrow=n)
  Mat_ind <- diag(n)
  
  C <- Mat_ind  - (1/(n)) * Mat_1
  
  # Matrice G des produits scalaires.
  G <- -0.5 * C %*% D2 %*% C
  
  # Valeurs propres et vecteurs propres
  Diago <- eigen(G)
  
  return(Diago)
}

approximation <- function(Donnee)
{ 
  # Valeurs propres
  dg<- diagonalisation(Donnee)
  val_p <-  dg$values
  
  # Le nombre des valeurs propres
  n <- length(val_p)
  
  # Le remplacement des valeurs propres négatives par des zéros
  for (i in 1:n) 
  {
    if(val_p[i]<0){ val_p[i] = 0}
    else{val_p[i] = val_p[i]}
    
  }
  return(val_p) 
}

quantite_dinfo <- function(k){
  
  # Les nouvelles valeurs propres
  new_valp <- approximation(Donnee)
  
  # La quantité d'information  non-euclidienne
  R2 <- sum(new_valp[1:k])/sum(new_valp)
  
  return(info) 
}

qualite_repr <- function(Donnee)
{
  
  dg<- diagonalisation(Donnee)
  val_p <-  dg$values
  new_valp <- approximation(Donnee)
  
  # calcul de la qualité d'approximation  euclidienne 
  quali <- sum(new_valp^2)/sum(val_p^2)
  
  return(quali) 
}

coordonnees <- function(Donnee){
  
  new_valp <- approximation(Donnee)
  dg <- diagonalisation(Donnee)
  
  # Matrice Lambda^1/2 pour 2 dimensions
  Lambda_0.5 <- diag( c(sqrt(new_valp[1]),sqrt(new_valp[2])),2)
  
  # Vecteurs propres correspondants
  vect_propres <- dg$vectors[,1:2]
  
  # Solution
  Y <- vect_propres %*% Lambda_0.5
  
  row.names(Y) <- colnames(Donnee)
  
  return(Y)
}

Y <- coordonnees(Donnee)

axe2 <-  Y[,2]
axe1 <-  - Y[,1]


plot(axe1,axe2,type="p", xlim=c(-1500,1500) , ylim=c(-700,700), main = "Représentation graphique des 10 villes américaines")

plot(axe1,axe2,type="p")
text(axe1,axe2,colnames(Donnee),pos=1,cex=0.8, col="blue")