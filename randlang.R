library("Matrix")
library("FactoMineR")
library("igraph")

load("bigwals-sample.RData")

data <- read.delim("bigwals-data.txt", row.names = 1, check.names = F)
booldata <- 1 - is.na(data)
tokeep <- data[rowSums(booldata) > ncol(data)/2,]
langdata <- read.delim("bigwals-langs.txt", row.names = 4)
langdata <- langdata[dimnames(tokeep)[[1]],]
cltreedata <- as.matrix(read.csv("bigwals-cltree.csv",row.names = 1))
cltree <- graph_from_edgelist(cltreedata[,1:2])
E(cltree)$weight <- as.numeric(cltreedata[,3])
weights <- c(1,unlist(sapply(colnames(tokeep),function (x) incident(cltree,x,mode = "in")$weight)))

ftokeep <- as.data.frame(lapply(tokeep,as.character))
levels <- lapply(ftokeep, levels)

dtokeep <- tab.disjonctif.prop(tokeep)
dtokeep[dtokeep > 0 & dtokeep < 1] <- NA
correl <- 2*sin((pi/6)*cor(dtokeep, use = "pairwise.complete.obs"))
correl[is.na(correl)] <- 0
correl <- nearPD(correl,T)$mat
transcorrel <- chol(correl)

cooklang <- function (){
  standard <- numeric()
  langs <- c()
  bgfac <- readline("Balanced, global, family, area, or country? [bgfac] ")
  
  if (bgfac == "b"){
    langs <- result
  } else if (bgfac == "g"){
    langs <- rownames(tokeep)
  } else if (bgfac == "f"){
    f <- readline(cat("Which family? (Type the name in full.) E.g.:\n",head(names(sort(table(langdata[,1]), decreasing = T))),sep = "\n"))
    langs <- langdata[,1] == f
  } else if (bgfac == "a"){
    a <- readline(cat("Which area? (Type the name in full.) E.g.:\n",names(sort(table(langdata[,2]), decreasing = T)),sep = "\n"))
    langs <- langdata[,2] == a
  } else if (bgfac == "c"){
    c <- readline(cat("Which country? (Type the name in full.) E.g.:\n",head(names(sort(table(langdata[,3]), decreasing = T))),sep = "\n"))
    langs <- langdata[,3] == c
  }
  
  standard <- colMeans(dtokeep[langs,],na.rm = T)
  randfeats <- apply(rnorm(ncol(dtokeep)) %*% transcorrel,1,pnorm) <= standard
  
  for (i in sample(which(randfeats == 1 & colnames(dtokeep) > 0))){
    # message(sum(colMeans((dtokeep[,i] - dtokeep) > 0, na.rm = T) == 0 & colnames(dtokeep) > 0)-1)
    randfeats[which(colMeans((dtokeep[,i] - dtokeep) > 0, na.rm = T) == 0 & colnames(dtokeep) > 0)] <- 2
    randfeats[i] <- 1
  }
  
  randlang <- numeric(ncol(tokeep))
  names(randlang) <- dimnames(tokeep)[[2]]
  thiscol <- 1
  for (l in seq_along(names(levels))){
    randlang[l] <- as.numeric(levels[l][[1]][which.max(randfeats[thiscol:(thiscol+length(levels[l][[1]])-1)])])
    thiscol <- thiscol+length(levels[l][[1]])
  }
  
  randlang <- as.matrix(randlang)
  
  comp <- as.matrix(sort(rowSums(abs(t(colMeans(tokeep[langs,],na.rm = T))/apply(tokeep,2,max,na.rm=T) - (tokeep[langs,]/apply(tokeep,2,max,na.rm=T)))*matrix(weights, nrow = 1), na.rm = T), decreasing = T))
  comp <- names(comp[1:(length(comp)/2),])
  
  langdiff <- rowMeans(abs(t(randlang)/apply(tokeep,2,max,na.rm=T) - (tokeep[comp,]/apply(tokeep,2,max,na.rm=T)))*weights, na.rm = T)
  message(paste("Your language resembles ",names(which.min(langdiff))," (distance = ",min(langdiff),").",sep = ""))
  
  randlang
}

save(result,tokeep,langdata,weights,levels,dtokeep,transcorrel,file = "randlang-sample.RData")