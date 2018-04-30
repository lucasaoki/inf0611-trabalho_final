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

db    <- read_csv("SwedishLeaf_TRAIN.csv")
query <- read_csv("SwedishLeaf_TEST.csv")
# normaliza cada série, mas retorna cada série por coluna
db_norm <- apply(db[-1], 1, normalize)

