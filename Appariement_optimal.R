#Code R: Application sur distances entre séquences obtenues par appariement optimal

# Charger TraMineR 
library(TraMineR)

# Charger le jeu de donnée mvad
data(mvad)
mvad.labels <- c("employment", "further education", "higher education",
                 "joblessness", "school", "training")
mvad.scodes <- c("EM","FE","HE","JL","SC","TR")

# Céer un objet séquences d'états
mvad.seq <- seqdef(mvad, 17:86, left="DEL")
seqistatd(mvad.seq[1:10, ])

# visualiser les 10 premières séquences
seqiplot(mvad.seq, with.legend= "right")

# La matrice des coûts de substitution 
couts2 <- seqsubm(mvad.seq,method="CONSTANT", cval=2)

# Le calcul de la matrice de distance entre les séquences en utilisant l'optimal matching
mvad.om2 <- seqdist(mvad.seq, method = "OM", indel = 1, sm = couts2)

# Effectuer MDS   
Y <- coordonnees(mvad.om2)
axe2 <- - Y[,2]
axe1 <- - Y[,1]
plot(axe1,axe2,type="p", xlim=c(-50,90) , ylim=c(-50,75))
text(axe1,axe2,colnames(mvad.om2),pos=1,cex=0.8)

# calculerl'indicateur de qualité C2
qualite_repr(mvad.om2)

# Effectuer une classification avec l'indice de Ward
seq.dist <- hclust(as.dist(mvad.om2), method = "ward.D2")

# Tracer le dendrogramme
plot(as.dendrogram(seq.dist), leaflab = "none", 
     main = "Dendrogramme de la hiérarchie basé sur l'indice de Ward")

# déterminer les sauts d'inertie.
R.inertie <- sort(seq.dist$height, decreasing = TRUE)
plot(R.inertie[1:20], type = "s", xlab = "Nombre de classes", ylab = "inertie")
points(c(2, 4), R.inertie[c(2, 4)], col = c( "red3", "blue3"), cex = 2, lwd = 3)

# une partition en 4 classes
nbcl = 4
seq.part <- cutree(seq.dist, nbcl)
seq.part <- factor(seq.part, labels = paste("classe", 1:nbcl, sep = "."))

# Couper le Dendrogramme de la hiérarchie en 4 classes 
plot(as.dendrogram(seq.dist), leaflab = "none",
     main = "Dendrogramme de la hiérarchie coupé en 4 classes ")
rect.hclust(seq.dist, 4, border = "green3")

# Appliquer le k-means
Centre <- apply ( mvad.om2 ,2,function ( x ) tapply (x , seq.part  , mean ))
seq.km <- kmeans ( mvad.om2 ,Centre )
seq.cluster <- factor(seq.km$cluster, labels = paste("classe", 1:nbcl, sep = "."))
table(seq.cluster)

# Tracer les clusters
library(factoextra)
fviz_cluster(seq.km, mvad.om2,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#999999"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)
# Déterminer le R2
(seq.km$betweenss/seq.km$totss)*100

# Tracer les chronogrammes
seqdplot(mvad.seq , group = seq.part, border = NA, main = "Chronogrammes")