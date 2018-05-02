library(ggplot2)
#setwd("/Users/seelengut/courses/machlearn/inf-0611/tfinal/inf0611-trabalho_final")
trainSet <- read.csv("SwedishLeaf_TRAIN.csv", header = FALSE, sep = ",")
testSet <- read.csv("SwedishLeaf_TEST.csv", header = FALSE, sep = ",")
#unname(unlist(trainSet[trainSet[[1]] == 4,][1,]))

recurrencePlot <- function(data, epsilon=0.1) {
  m <- matrix(ncol = length(data), nrow = length(data))
  for (i in c(1:length(data))) {
    m[i,] <- abs(data[i] - data)
  }
  return(m < epsilon)
}

calcLbpPos <- function(data, i, j) {
  dimenLen <- length(data[1,])
  multVect <- c(1,2,4,8,16,32,64,128)
  v <- rep(T,8)
  if (i == 1) {
    v <- v & c(F,F,T,T,T,T,T,F)
  }
  if (i == dimenLen) {
    v <- v & c(T,T,T,F,F,F,T,T)
  }
  if (j == 1) {
    v <- v & c(T,T,T,T,T,F,F,F)
  }
  if (j == dimenLen) {
    v <- v & c(T,F,F,F,T,T,T,T)
  }
  v[1] <- v[1] && data[i-1, j]
  v[2] <- v[2] && data[i-1, j+1]
  v[3] <- v[3] && data[i, j+1]
  v[4] <- v[4] && data[i+1, j+1]
  v[5] <- v[5] && data[i+1, j]
  v[6] <- v[6] && data[i+1, j-1]
  v[7] <- v[7] && data[i, j-1]
  v[8] <- v[8] && data[i-1,j-1]
  return(sum(v*multVect))
}

lbpMatrix <- function(data) {
  dataLen <- length(data[1,])
  multVect <- c(1,2,4,8,16,32,64,128)
  m <- matrix(ncol = dataLen, nrow = dataLen)
  for (i in c(1:dataLen)) {
    for (j in c(1:dataLen)) {
      if (!data[i,j]) {
        m[i,j] <- 255
      } else {
        m[i,j] <- calcLbpPos(data,i,j)
      }
    }
  }
  return(m)
}

lbp <- function(data) {
  m <- lbpMatrix(data)
  v <- rep(0,256)
  for (i in c(1:length(m[1,]))) {
    for (j in c(1:256)) {
      v[j] <- v[j] + sum(m[i,] == (j-1))
    }
  }
  return(v)
}

l2Dist <- function(v1, v2) {
  return(sqrt(sum((v1-v2)^2)))
}

l1Dist <- function(v1, v2) {
  return(sum(abs(v1-v2)))
}

applyStrategy1 <- function(dataSet) {
  result <- matrix(ncol = 257, nrow = nrow(dataSet))
  for (i in c(1:nrow(dataSet))) {
    v <- unname(unlist(dataSet[i,]))
    tag <- v[1]
    v <- v[2:length(v)]
    v <- lbp(recurrencePlot(v,0.2))
    v <- c(tag,v)
    result[i,] <- v
  }
  return(result)
}

getClosest <- function(train, test, distFunc) {
  result <- matrix(nrow=nrow(test), ncol=nrow(train)+1)
  for (i in c(1:nrow(test))) {
    v <- c(test[i,1])
    t <- test[i,][2:ncol(test)]
    for (j in c(1:nrow(train))) {
      v <- c(v,distFunc(train[j,][2:ncol(train)],t))
    }
    result[i,] <- v
  }
  return(result)
}

# Função que calcula a precisão com base nas distancias calculadas
# retornando o resultado ordenado pelas distâncias
precisionForKFirst <- function(distances, ids, relevant_id, k) {
  ord_ids <- ids[order(distances)][1:k]
  #print(ord_ids)
  max_occurrences <- min(k,sum(ids == relevant_id))
  total <- 0
  relevant <- 0
  for (id in ord_ids) {
    total <- total + 1
    if (relevant_id == id) {
      relevant <- relevant + 1
    }
    #if (relevant == max_occurrences) break
  }
  return(relevant/total)
}

averagePrecisionForKFirst <- function(distances, ids, relevant_id, k) {
  ord_ids <- ids[order(distances)][1:k]
  #print(ord_ids)
  max_occurrences <- min(k, sum(ids == relevant_id))
  total <- 0
  relevant <- 0
  meanSum <- 0
  for (id in ord_ids) {
    total <- total + 1
    if (relevant_id == id) {
      relevant <- relevant + 1
      meanSum <- meanSum + relevant/total
    }
    #if (relevant == max_occurrences) break
  }
  return(meanSum/total)
}

# Função que calcula a revocação com base nas distancias calculadas
# retornando o resultado ordenado pelas distâncias
recallForKFirst <- function(distances, ids, relevant_id, k) {
  ord_ids <- ids[order(distances)][1:k]
  #total <- min(k, sum(ids == relevant_id))
  total <- sum(ids == relevant_id)
  relevant <- 0
  for (i in c(1:length(ord_ids))) {
    if (relevant_id == ord_ids[i]) {
      relevant <- relevant + 1
    }
    if (relevant/total == 1) break
  }
  return(relevant/total)
}

meanPrecisionRecall <- function(distances, ids, k) {
  precision <- c()
  recall <- c()
  for (i in k) {
    p <- c()
    r <- c()
    for (row in c(1:nrow(distances))) {
      #p <- c(p, precisionForKFirst(distances[row,][2:ncol(distances)], ids, distances[row,][1], i))
      p <- c(p, averagePrecisionForKFirst(distances[row,][2:ncol(distances)], ids, distances[row,][1], i))
      r <- c(r, recallForKFirst(distances[row,][2:ncol(distances)], ids, distances[row,][1], i))
    }
    #print(p)
    precision <- c(precision, mean(p))
    recall <- c(recall, mean(r))
  }
  return(data.frame(elements=k, mean_precision=precision, mean_recall=recall))
}

# Função que gera o gráfico de precisao e revocacao medias x numero de elementos considerados.
plotMeanPrecisionRecall <- function(dataFrame) {
  plotMeansByK(dataFrame, "Precisão e Revocação Médias x Elementos Considerados", c("Precisão","Revocação"))
}

plotMeansByK <- function(dataFrame, title, keys) {
  keys <- rep(keys, times=1, each=nrow(dataFrame))
  data <- c()
  for (i in c(2:ncol(dataFrame))) {
    data <- c(data, dataFrame[[i]])
  }
  
  plot <- data.frame(elements=dataFrame[[1]], keys = keys, means = data)
  g <- ggplot(data = plot, aes(x=elements, y=data, colour=keys)) + ylim(0,1)
  g <- g + geom_point() + geom_line()
  
  g <- g + labs(title=title, y="Percentual", x="Elementos Considerados", colour="") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12))
  return(g)
}


# Estratégia 1 - Recurrence Plot com 1 bit + LBP
# Comparação feita através dos valores acumulados de cada uma das 256 possibilidades
st1Train <- applyStrategy1(trainSet)
st1Test <- applyStrategy1(testSet)
st1ClosestL1 <- getClosest(st1Train, st1Test, l1Dist)
st1ClosestL2 <- getClosest(st1Train, st1Test, l2Dist)
st1MeanPrecRecL1 <- meanPrecisionRecall(st1ClosestL1, trainSet[,1], seq(5,100,5))
st1MeanPrecRecL2 <- meanPrecisionRecall(st1ClosestL2, trainSet[,1], seq(5,100,5))
plotMeanPrecisionRecall(st1MeanPrecRecL1)
plotMeanPrecisionRecall(st1MeanPrecRecL2)

# Gráfico comparativo das revocações médias entre ambas estratégias e distâncias
meanRecallDf <- data.frame(elements=st1MeanPrecRecL1$elements, st1_l1 = st1MeanPrecRecL1$mean_recall, st1_l2 = st1MeanPrecRecL2$mean_recall)
plotMeansByK(meanRecallDf, "Comparação das Revocações Médias", c("Estratégia 1 - L1", "Estratégia 1 - L2"))

# Gráfico comparativo das precisões médias entre ambas estratégias e distâncias
meanPrecisionDf <- data.frame(elements=st1MeanPrecRecL1$elements, st1_l1 = st1MeanPrecRecL1$mean_precision, st1_l2 = st1MeanPrecRecL2$mean_precision)
plotMeansByK(meanPrecisionDf, "Comparação das Precisões Médias", c("Estratégia 1 - L1", "Estratégia 1 - L2"))
