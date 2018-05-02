########################################
# Trabalho Final - INF-0611          
# Nome(s):  Felipe Wolff Ramos
#           Lucas Aoki Heredia
########################################

library(dtw)

read_csv <- function(filename) {
  db <- read.csv(filename, header = FALSE)
  series_size <- length(db[1, -1])
  colnames(db) <- c("especie", (1:series_size))
  return(db)
}

# funcao para normalizar
normalize <- function(serie) {
  serie <- as.numeric(serie)
  return((serie - mean(serie) ) / sd(serie))
}

dtw_dists<- function(ss1, ss2, dist.method = "Euclidean") {
  dist_matrix <- matrix(nrow = ncol(ss1), ncol = ncol(ss2))
  
  for (i in seq(ncol(ss1))) {
    for (j in seq(ncol(ss2))) {
      dist_matrix[i, j] = dtw(ss1[, i], ss2[, j], dist.method = dist.method, distance.only = TRUE)$distance
    }
  }
  return(dist_matrix)
}

closest_species <- function(db_species, dist_matrix) {
  species_order <- matrix(nrow = nrow(dist_matrix), ncol = ncol(dist_matrix))
  
  for (i in seq(nrow(dist_matrix))) {
    species_order[i,] <- db_species[order(dist_matrix[i,])]  
  }
  return(species_order)
}

# Calculo da Precisao
precision <- function(retrieved, relevant) {
  rate <- length(retrieved[retrieved == unique(relevant)]) / length(retrieved)
  return (rate)
}

# Calculo da Revocacao
recall <- function(retrieved, relevant) {
  rate <- length(retrieved[retrieved == unique(relevant)]) / length(relevant)
  return(rate)
}

# Funcao que aplica recall ou precision definido por fn
check_relevant <- function(ids, id_rel, fn, k) {
  if (is.null(k)) {
    k <- length(ids)
  }
  v <- c()
  for (i in k) {
    v <- c(v, fn(ids[1:i], id_rel))
  }
  return(v)
}

mean_precision_recall <- function(retrieved_ids, query_ids, k) {
  prec <- matrix(nrow = length(k), ncol = length(query_ids))
  rec  <- matrix(nrow = length(k), ncol = length(query_ids))
  for (i in seq(length(query_ids))) {
    prec[,i] <- check_relevant(retrieved_ids[i,], retrieved_ids[i, retrieved_ids[i,] == query_ids[i]], precision, k)
    rec[,i] <- check_relevant(retrieved_ids[i,], retrieved_ids[i, retrieved_ids[i,] == query_ids[i]], recall, k)
  }  
  return(data.frame(elements=k, mean_precision=rowMeans(prec), mean_recall=rowMeans(rec)))
}

db    <- read_csv("SwedishLeaf_TRAIN.csv")
query <- read_csv("SwedishLeaf_TEST.csv")
# normaliza cada série, mas retorna cada série por coluna
db_norm <- apply(db[-1], 1, normalize)
db_species <- db$especie
query_norm <- apply(query[-1], 1, normalize)
query_species <- query$especie
k <- seq(5, 100, 5)

# Returns distances matrix 625x500
# dtw_euclidean <- dtw_dists(query_norm, db_norm)
# dtw_closest_euclidean <- closest_species(db_species, dtw_euclidean)
# dtw_L1 <- dtw_dists(query_norm, db_norm, "L1")
# dtw_closest_L1 <- closest_species(db_species, dtw_L1)
dtw_L2 <- dtw_dists(query_norm, db_norm, "L2")
dtw_closest_L2 <- closest_species(db_species, dtw_L2)


##### função plotMeanPrecisionRecall roubada do seu código... tem que estar sourced de lá
# plotMeanPrecisionRecall(mean_precision_recall(dtw_closest_euclidean, query_species, k))
# plotMeanPrecisionRecall(mean_precision_recall(dtw_closest_L1, query_species, k))
plotMeanPrecisionRecall(mean_precision_recall(dtw_closest_L2, query_species, k))
