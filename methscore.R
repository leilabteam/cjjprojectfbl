#!/usr/bin/env Rscript

setwd('E:/Rpractice/陈佳佳核糖体测序/fbl14')
W2xGFP_15end = read.table('5W2xGFP_1.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)

side_weight <- seq(0.5,1,by=0.1)

out = as.data.frame(matrix(nrow=0,ncol=4))
colnames(out) = c('end','gene','site','score')
nn=0
gg = unique(data[,1])
vec1 <- rep(0, times=1786) #gg[1]
vec2 = rep(0, times=1482)  ##gg[2]
vec3 =  rep(0, times=1306) ##gg[3]
tmp1 = data[data[,1]==gg[1],]
tmp2 = data[data[,1]==gg[2],]
tmp3 = data[data[,1]==gg[3],]
for(i in 1:dim(tmp1)[1]){
  s1 = tmp1[i,2]
  s2 = tmp1[i,3]-1
  for(j in s1:s2){
    vec1[j+1] = tmp1[i,4]
  }
}

for(i in 1:dim(tmp2)[1]){
  s1 = tmp2[i,2]
  s2 = tmp2[i,3]-1
  for(j in s1:s2){
    vec2[j+1] = tmp2[i,4]
  }
}

for(i in 1:dim(tmp3)[1]){
  s1 = tmp3[i,2]
  s2 = tmp3[i,3]-1
  for(j in s1:s2){
    vec3[j+1] = tmp3[i,4]
  }
}
# scoreC
for(j in 1:length(vec1)){
  nn = nn +1
  counts = vec1
  if(6< j & j < length(vec1)-6){
    b <- c(counts[j-6],counts[j-5],counts[j-4],counts[j-3],counts[j-2],counts[j-1])
    a <- c(counts[j+6],counts[j+5],counts[j+4],counts[j+3],counts[j+2],counts[j+1])
    b <- sum(b*side_weight)/sum(side_weight)
    a <- sum(a*side_weight)/sum(side_weight)
    if(a>0 & b>0){
      methscore1 <- counts[j]*2/(a+b)
      # methscore <- c(methscore,methscore1)
    } 
    else methscore1 <- 1
  }
  else methscore1 <- 1
  out[nn,1] = '3end'
  out[nn,2] = gg[1]
  out[nn,4] = 1-methscore1
  out[nn,3] = j
}

for(j in 1:length(vec2)){
  nn = nn +1
  counts = vec2
  if(6< j & j < length(vec2)-6){
    b <- c(counts[j-6],counts[j-5],counts[j-4],counts[j-3],counts[j-2],counts[j-1])
    a <- c(counts[j+6],counts[j+5],counts[j+4],counts[j+3],counts[j+2],counts[j+1])
    b <- sum(b*side_weight)/sum(side_weight)
    a <- sum(a*side_weight)/sum(side_weight)
    if(a>0 & b>0){
      methscore1 <- counts[j]*2/(a+b)
      # methscore <- c(methscore,methscore1)
    }
    else methscore1 <- 1
  }
  else methscore1 <- 1
  out[nn,1] = '3end'
  out[nn,2] = gg[2]
  out[nn,4] = 1-methscore1
  out[nn,3] = j
}

for(j in 1:length(vec3)){
  nn = nn +1
  counts = vec3
  if(6< j & j < length(vec3)-6){
    b <- c(counts[j-6],counts[j-5],counts[j-4],counts[j-3],counts[j-2],counts[j-1])
    a <- c(counts[j+6],counts[j+5],counts[j+4],counts[j+3],counts[j+2],counts[j+1])
    b <- sum(b*side_weight)/sum(side_weight)
    a <- sum(a*side_weight)/sum(side_weight)
    if(a>0 & b>0){
      methscore1 <- counts[j]*2/(a+b)
      # methscore <- c(methscore,methscore1)
    }
    else methscore1 <- 1
  }
  else methscore1 <- 1
  out[nn,1] = '3end'
  out[nn,2] = gg[3]
  out[nn,4] = 1-methscore1
  out[nn,3] = j
}

W2xGFP_15end = read.table('5W2xGFP_1.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)
data = read.table('5WT18send.bedGraph',sep="\t",stringsAsFactors=F)



W2xGFP_15end

gg = unique(W2xGFP_15end[,1])
vec1 <- rep(0, times=1763) #gg[1]
vec2 = rep(0, times=1053)  ##gg[2]
vec3 =  rep(0, times=1580) ##gg[3]
tmp1 = W2xGFP_15end[W2xGFP_15end[,1]==gg[1],]
tmp2 = W2xGFP_15end[W2xGFP_15end[,1]==gg[2],]
tmp3 = W2xGFP_15end[W2xGFP_15end[,1]==gg[3],]
for(i in 1:dim(tmp1)[1]){
  s1 = tmp1[i,2]
  s2 = tmp1[i,3]-1
  for(j in s1:s2){
    vec1[j+1] = tmp1[i,4]
  }
}

for(i in 1:dim(tmp2)[1]){
  s1 = tmp2[i,2]
  s2 = tmp2[i,3]-1
  for(j in s1:s2){
    vec2[j+1] = tmp2[i,4]
  }
}

for(i in 1:dim(tmp3)[1]){
  s1 = tmp3[i,2]
  s2 = tmp3[i,3]-1
  for(j in s1:s2){
    vec3[j+1] = tmp3[i,4]
  }
}
# scoreC
for(j in 1:length(vec1)){
  nn = nn +1
  counts = vec1
  if(6< j & j < length(vec1)-6){
    b <- c(counts[j-6],counts[j-5],counts[j-4],counts[j-3],counts[j-2],counts[j-1])
    a <- c(counts[j+6],counts[j+5],counts[j+4],counts[j+3],counts[j+2],counts[j+1])
    b <- sum(b*side_weight)/sum(side_weight)
    a <- sum(a*side_weight)/sum(side_weight)
    if(a>0 & b>0){
      methscore1 <- counts[j]*2/(a+b)
      # methscore <- c(methscore,methscore1)
    } 
    else methscore1 <- 1
  }
  else methscore1 <- 1
  out[nn,1] = '5end'
  out[nn,2] = gg[1]
  out[nn,4] = 1-methscore1
  out[nn,3] = j
}

for(j in 1:length(vec2)){
  nn = nn +1
  counts = vec2
  if(6< j & j < length(vec2)-6){
    b <- c(counts[j-6],counts[j-5],counts[j-4],counts[j-3],counts[j-2],counts[j-1])
    a <- c(counts[j+6],counts[j+5],counts[j+4],counts[j+3],counts[j+2],counts[j+1])
    b <- sum(b*side_weight)/sum(side_weight)
    a <- sum(a*side_weight)/sum(side_weight)
    if(a>0 & b>0){
      methscore1 <- counts[j]*2/(a+b)
      # methscore <- c(methscore,methscore1)
    }
    else methscore1 <- 1
  }
  else methscore1 <- 1
  out[nn,1] = '5end'
  out[nn,2] = gg[2]
  out[nn,4] = 1-methscore1
  out[nn,3] = j
}

for(j in 1:length(vec3)){
  nn = nn +1
  counts = vec3
  if(6< j & j < length(vec3)-6){
    b <- c(counts[j-6],counts[j-5],counts[j-4],counts[j-3],counts[j-2],counts[j-1])
    a <- c(counts[j+6],counts[j+5],counts[j+4],counts[j+3],counts[j+2],counts[j+1])
    b <- sum(b*side_weight)/sum(side_weight)
    a <- sum(a*side_weight)/sum(side_weight)
    if(a>0 & b>0){
      methscore1 <- counts[j]*2/(a+b)
      # methscore <- c(methscore,methscore1)
    }
    else methscore1 <- 1
  }
  else methscore1 <- 1
  out[nn,1] = '5end'
  out[nn,2] = gg[3]
  out[nn,4] = 1-methscore1
  out[nn,3] = j
}

write.table(out,file='score.xls',col.names=T,row.names=F,sep="\t",quote=F)

print(a)
print(b)

