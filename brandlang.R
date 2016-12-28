library(pacman)
p_load(Matrix,FactoMineR,igraph,Hmisc,clusterGeneration)

load("bigwals-sample.RData")

bwdata <- read.delim("bigwals-data.txt", row.names = 1, check.names = F)
booldata <- 1 - is.na(bwdata)
bwtokeep <- bwdata[rowSums(booldata) > ncol(bwdata)/2,]
langdata <- read.delim("bigwals-langs.txt", row.names = 12)
langdata <- langdata[rownames(bwtokeep),]
featdata <- read.delim("bigwals-feature-labels.txt", header = F)
bwcltreedata <- as.matrix(read.csv("bigwals-cltree.csv",row.names = 1))
bwcltree <- graph_from_edgelist(bwcltreedata[,1:2])
E(bwcltree)$weight <- as.numeric(bwcltreedata[,3])
bwweights <- c(1,unlist(sapply(colnames(bwtokeep),function (x) incident(bwcltree,x,mode = "in")$weight)))

tokeep <- bwtokeep
weights <- bwweights

ftokeep <- as.data.frame(lapply(tokeep,as.character))
levels <- lapply(ftokeep, levels)

dtokeep <- tab.disjonctif.prop(tokeep)
dtokeep[dtokeep > 0 & dtokeep < 1] <- NA

condposs <- matrix(1, nrow = ncol(dtokeep), ncol = ncol(dtokeep))
condnegs <- matrix(1, nrow = ncol(dtokeep), ncol = ncol(dtokeep))
# condprobs <- outer(1:ncol(dtokeep), 1:ncol(dtokeep), FUN = Vectorize(function(x,y){sum(dtokeep[dtokeep[,x] == 1, y],na.rm = T)/sum(dtokeep[,x] == 1 & !is.na(dtokeep[,y]),na.rm = T)}))
for (i in 1:(nrow(condposs)-1)){
  for (j in (i+1):ncol(condposs)){
    condposs[i,j] <- sum(dtokeep[dtokeep[,j] == 1, i],na.rm = T)/sum(dtokeep[,j] == 1 & !is.na(dtokeep[,i]),na.rm = T)
    condnegs[i,j] <- sum(dtokeep[dtokeep[,j] == 0, i],na.rm = T)/sum(dtokeep[,j] == 0 & !is.na(dtokeep[,i]),na.rm = T)
  }
  if (i %% 10 == 0){
    message(paste("Column",i,"of",nrow(condposs)-1))
  }
}
ratios <- condposs/(condposs + condnegs)
ratios[is.nan(ratios)] <- 0.5

prior <- colMeans(dtokeep,na.rm = T)
thiscol <- 1
for (l in 1:ncol(tokeep)){
  if (sum(prior[thiscol:(thiscol + length(levels[[l]]) - 1)]) == 0 | is.nan(sum(prior[thiscol:(thiscol + length(levels[[l]]) - 1)]))){
    prior[thiscol:(thiscol + length(levels[[l]]) - 1)] <- 1/length(thiscol:(thiscol + length(levels[[l]]) - 1))
  }
  pick <- sample(1:length(levels[[l]]), 1, prob = prior[thiscol:(thiscol + length(levels[[l]]) - 1)])
  prior[thiscol:(thiscol + length(levels[[l]]) - 1)] <- 0
  prior[thiscol + pick - 1] <- 1
  if (thiscol+length(levels[[l]]) < length(prior)){
    remainder <- (thiscol+length(levels[[l]])):length(prior)
    ratio <- ratios[thiscol+pick-1,remainder]
    prior[remainder] <- (ratio * prior[remainder]) / ((ratio * prior[remainder]) + ((1 - ratio) * (1 - prior[remainder])))
  }
  thiscol <- thiscol + length(levels[[l]])
  # message(colnames(tokeep)[l])
}

randlang <- numeric(ncol(tokeep))
randlangfeats <- character(ncol(tokeep))
names(randlang) <- colnames(tokeep)
names(randlangfeats) <- colnames(tokeep)
thiscol <- 1
for (l in 1:ncol(tokeep)){
  randlangfeats[l] <- as.character(featdata[featdata[,1] == colnames(tokeep)[l],2][which(prior[thiscol:(thiscol+length(levels[[l]])-1)] == 1)])
  randlang[l] <- which(prior[thiscol:(thiscol+length(levels[[l]])-1)] == 1) - 1
  thiscol <- thiscol + length(levels[[l]])
}

langdiff <- rowMeans(abs(t(as.matrix(randlang))/apply(tokeep,2,max,na.rm=T) - (tokeep/apply(tokeep,2,max,na.rm=T)))*weights, na.rm = T)
message(paste("Your language resembles ",names(which.min(langdiff))," (distance = ",min(langdiff),").",sep = ""))
as.matrix(randlangfeats)

save(result,tokeep,langdata,featdata,weights,levels,dtokeep,likelihoods,file = "randlang-sample.RData")