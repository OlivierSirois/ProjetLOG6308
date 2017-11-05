m = as.matrix(read.table("citeseer.rtable", check.names=F))
library("expm")
library("dplyr")
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
  res[res > 1] <- 1

  return(res)
}

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

#diag(m) <- 0


# On calcule notre domaine S comme étant tout les article sont référencés par notre origine. Pour faire cela, on regarde ceux qui on une valeur positive dans notre matrice référentielle
S <- which(m["422908",]==1)
#print(S)
#print(m[S,S])


pr.S <- rep(1, dim(m[S,S])[1])
#on calcule le page rank local du domaine S
pr.S <- page.rank.until.stab(m[S,S], d, pr.S)
#print(pr.S)
# Pour le domaine S prime, nous voulons aussi rajouter les références des références. Pour faire cela, nous faisont que prendre la somme des deux première puissance de la
# matrice référentielle. C'est a dire, la matrice référentielle elle-même et la deuxième puissance.
m.prime <- sum.powers.matrix(m,2)

#on enleve les auto-references
diag(m.prime) <- 0

S.prime <- which(m.prime["422908",]==1)
#on calcule le page rank du domaine S prime
pr.S.prime <- rep(1, dim(m[S.prime,S.prime])[1])

pr.S.prime <- page.rank.until.stab(m[S.prime, S.prime], d, pr.S.prime)

#S.prime <- remove.autoreferences(S.prime)


#On calcule le pagerank de tout nos articles (pas très long)
pr <- page.rank.until.stab(m, d, pr)
#on place les rankings PageRank de notre domaine dans un vecteur
S.rankings <- pr[S]
#pour le domaine S prime
S.prime.rankings <- pr[S.prime]


#On crée un dataframe avec nos données
S.loc.dat <- data.frame(pr.S)
S.loc.dat$article <- rownames(pr.S)
#meme chose pour S prime (local)
S.prime.loc.dat <- data.frame(pr.S.prime)
S.prime.loc.dat$article <- rownames(pr.S.prime)
#on trie pour S et S prime
S.loc.best <- S.loc.dat %>% select(pr.S, article) %>% arrange(desc(pr.S, arr.ind=T))
S.prime.loc.best <- S.prime.loc.dat %>% select(pr.S.prime, article) %>% arrange(desc(pr.S.prime, arr.ind=T))

# on effectue les meme calcules pour PageRank de maniere globale
S.dat <- data.frame(S.rankings, S)
S.dat$article = rownames(S.dat)
#même chose pour le domaine S prime
S.prime.dat <- data.frame(S.prime.rankings, S.prime)
S.prime.dat$article = rownames(S.prime.dat)
#on effectue un trie sur nos valeurs, on sort celles qui sont les plus hautes en premier
S.best <- S.dat %>% select(S.rankings, S, article) %>% arrange(desc(S.rankings, arr.ind=T))
#même chose pour le domaine S prime
S.prime.best <- S.prime.dat %>% select(S.prime.rankings, S.prime, article) %>% arrange(desc(S.prime.rankings, arr.ind=T))

## on calcule nos coefficients de correlation et du cosinus

corr.ratings <- corr.vm(m["422908", ], m[S,])
cos.ratings <- cosinus.vm(m["422908",], t(m[S,]))

# on remplace les NaN par 0
corr.ratings[is.nan(corr.ratings)] <- 0
cos.ratings[is.nan(cos.ratings)] <- 0

# on prends nos etiquettes
labels.corr <- colnames(corr.ratings)
labels.cos <- colnames(cos.ratings)

# on cree nos dataframes
df.corr <- data.frame(corr = as.vector(corr.ratings), article = labels.corr)
df.cos <- data.frame(cos = as.vector(cos.ratings), article = labels.cos)

#on trie
df.best.cos <- df.cos %>% select(cos, article) %>% arrange(desc(cos, arr.ind=T))
df.best.corr <- df.corr %>% select(corr, article) %>% arrange(desc(corr, arr.ind=T))

#print.data.frame(df.best.cos)
#print.data.frame(df.best.corr)

#print.data.frame(S.prime.best)
#print.data.frame(S.best)
#print(cl)

# ---------------------------- Validation croiée -------------------------------

# Calcul du cosinus entre un vecteur et les colonnes d'une matrice.
# Méthode adaptée pour résoudre les problèmes liés aux colonnes remplies de 0.
cosinus.vm <- function(v,m) {
  # On on met tous nos valeurs de NA a 0, sinon on va avoir des problemes de calculs avec des matrices sparse
  m[is.na(m)] <- 0
  v[is.na(v)] <- 0 ;
  # On calcule le cosinus entre le vecteur V et les colonnes de la matrice m en utilisant la formule vu en classe
  denom <- sqrt(colSums(m^2)) * sqrt(sum(v^2))
  denom[denom == 0] <- 1
  res <- (v %*% m) / denom

  return(res)
}

# Retourne les indices des n-premières valeurs du vecteur (ordre décroissant).
max.nindex <- function(m, n=5) {
    i <- order(m, decreasing=TRUE)
    return(i[1:n])
}

# Prédiction d'une valeur dans la matrice de références (approche item-item).
predict.value <- function(refs, user, item)
{
    # Moyenne de chaque item
    items.avg <- colMeans(refs, na.rm=T)

    # Similarité entre l'item étudié et les autres items
    items.similarity <- cosinus.vm(refs[,item], refs)

    # Sélection des 20 premiers voisins en fonction de leur similarité
    # (la distance n'était pas pertinente ici, vu le nombre de 0 dans la matrice
    # de références)
    n.voisins <- 20 + 1
    most.similar <- max.nindex(items.similarity, n.voisins)[-item]

    # Calcul du facteur de correction (inverse de la somme des similarités).
    # Prise en compte du cas où cette somme est nulle.
    correction.factor <- 1 / sum(items.similarity[most.similar])
    if (correction.factor == Inf)
    {
        correction.factor <- 1
    }

    # Prédiction de la valeur
    value <- items.avg[item] + correction.factor *
             sum(sapply(most.similar, function(x)
                items.similarity[x] * (refs[user, x] - items.avg[x])), na.rm=T)

    return(value)
}

# Calcul de la RMSE entre une matrice de résultats et une matrice cible
rmse <- function(results, target)
{
    res <- sqrt(sum(abs(results^2 - target^2)) / dim(target)[1]^2)

    return(res)
}

# Proportion de références de test (ici : 5%, les calculs sont vraiment longs
# pour 10%)
cross.validation.factor <- 0.05

# Nombre d'articles dans la base
nb.articles <- dim(m)[1]

# Sélection aléatoire des indices des références de la base de test
test.refs.row.indices <- sample(1:nb.articles, round(cross.validation.factor * nb.articles))
test.refs.col.indices <- sample(1:nb.articles, round(cross.validation.factor * nb.articles))

# Séparation de la base d'entraînement et de la base de test
m.test <- m[test.refs.row.indices, test.refs.col.indices]
m.training <- m
m.training[test.refs.row.indices, test.refs.col.indices] <- NA

results <-sapply(test.refs.row.indices, function(x)
            sapply(test.refs.col.indices, function(y)
                predict.value(as.matrix(m.training), x, y)))

print(results)

cross.validation.rmse <- rmse(as.matrix(results), as.matrix(m.test))

print(cross.validation.rmse)
