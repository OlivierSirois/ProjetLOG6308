m = as.matrix(read.table("citeseer.rtable", check.names=F))
library("expm")
library("dplyr")
library("signal")
d <- 0.85
pr <- rep(1, dim(m)[1])

## Cosinus entre un vecteur v et chaque colonne dela matrice m
cosinus.vm <- function(v,m) {
  # On on met tous nos valeurs de NA a 0, sinon on va avoir des problemes de calculs avec des matrices sparse
  m[is.na(m)] <- 0
  v[is.na(v)] <- 0 ;
  # On calcule le cosinus entre le vecteur V et les colonnes de la matrice m en utilisant la formule vu en classe
  (v %*% m)/(sqrt(colSums(m^2)) * sqrt(sum(v^2)))
}

#Correlation entre la rangee v (v = index) et chaque colonne de la matrice m
corr.vm <- function(v,m) {
  # on centre les valeurs de la matrice m en fonction de la moyenne, on la renomme m.centre
  v.i <- rowMeans(m[,], na.rm=T)
  # on enleve les NA
  m[is.na(m)] <- 0
  m.centre <- m - v.i
  # on centre le vecteur v en fonction de sa moyenne
  v.index <- v
  v.index[is.na(v)] <- 0
  v.index <- v.index - mean(v, na.rm=T)
  # on retourne ensuite un vecteur correspondant entre le vecteur v et sa correlation avec chaque rangers de m
  return( (v.index%*%t(m.centre))/(sqrt(sum(v.index^2) * rowSums(m.centre^2))))
}

#cette fonction fais la somme de toutes les puissance de la matrice m allant de 0 jusqu'a n
sum.powers.matrix <- function(m, n) {
  powers <- c(1:n)
  res <- Reduce('+', lapply(powers, function(x) m %^% x))
  #res[res > 1] <- 1

  return(res)
}

#cette fonction fais la somme de toutes les puissance de la matrice m allant de 0 jusqu'a n
sum.powers.matrix.sinc <- function(m, n, freq) {
  if(n%%2 == 1){
    n <- n+1
  }
  powers <- c(1:(n/2))
  print(powers)
  x <- seq(-3*pi, 3*pi, length=n)
  
  x <- lapply(x, function(x) freq * sin(x*freq)/(x*pi) )
  #x <- replace(x, is.na(x), 1)
  plot(c(1:length(x)), x)
  sinc.half <- x[(n/2+1):n]
  
  #print(sinc.half[[1]])
  res.a <- Reduce('+', lapply(powers, function(x) {(m %^% x) * sinc.half[[x]]} ))
  res.b <- Reduce('+', lapply(powers, function(x) {(t(m) %^% x) * sinc.half[[x]]} ))
  res <- res.a + res.b
  #res <- Reduce('+', lapply(powers, function(x) {(res %^% x) * sinc.half[[x]]} ))
  print(res)
  
  return(res)
}

sinc <- function (x) sin(x) / x

#on enleve toutes les autoreferences dans la matrice referentielle. Ce probleme survient lorsqu'on fait une addition des différentes puissances de matrices.
#ce comportement est a évité car sa donne de l'importance artificielle à notre cible lorsqu'elle est référencé elle même par sa référence.
remove.autoreferences <- function(m) {
  res <- m
  diag(res) <- 0
  return (res)
}

# Cette fonction exécute l'algorithme PageRank jusqu'à ce que c'est valeurs soit stabilisés. On définie une stabilit. lorsque l'erreur moyenne absolue est moins de .0001
page.rank.until.stab <- function(m, d, pr){
  pr.next <- page.rank(m, d, pr)
  while(abs(mean(pr.next-pr)) > 0.0001){
    pr <- pr.next
    pr.next <- page.rank(m,d,pr)
  }
  return (pr)
}

#Cette fonction fait une itérations de l'algorithme page rank. Nous avons ajouter une petite modifications dans le cas ou la somme d'un colonne est de 0 (c'est à dire, une
# article non référencé) on remplace ensuite la valeur qui sera égale a Inf, par un 0
page.rank <- function(m, d, pr){
  denum <- (pr/colSums(m))
  denum[denum == Inf] <- 0
  (pr <- (1-d)/3 + (d * (m %*%denum)))
  #print(pr[450])
  return(pr)
}

pagerank.iteration <- function(refs, n, d, pr) {
  # Number of articles
  n.articles <- dim(refs)[1]

  # n-level references
  m <- sum.powers.matrix(refs, n)

  # Compute PageRank
  pr.res <- (1-d)/n.articles + (d * (m %*% (pr/colSums(m))))

  return(pr.res)
}

mat <- sum.powers.matrix.sinc(m, 15, 0.30)
