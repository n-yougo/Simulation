x <- read.table("https://github.com/n-yougo/dataset/raw/master/singhTest_x.txt")
y <- read.table("https://github.com/n-yougo/dataset/raw/master/singhTest_y.txt")
x <- as.matrix(x)
y <- as.matrix(y)
p <- dim(x)[2]      #ｘの次元

#前立腺ガンの人数を数える
z <- matrix(nrow = length(y),ncol = 1)
for (i in 1:length(y)) {
  if (y[i,] == -1) {
    z[i,] <- 1
  }
  else {
    z[i,] <- 0
  }
}

n_0 <- length(y)
n_1 <- sum(z)        #がん患者の人数
n_2 <- n_0 - n_1     #健康な人の人数


#行列 x をがん患者とそうでない人に分割する
x_1 <- matrix(nrow = n_1,ncol = p)     #がん患者
x_2 <- matrix(nrow = n_2,ncol = p)     #健康な人
for (i in 1:n_1) {
  x_1[i,] <- x[which(y == -1)[i],]
}
for (i in 1:n_2) {
  x_2[i,] <- x[which(y == 1)[i],]
}

#行列 x をがん患者と健康な人の両方のデータが入った２つの行列に分ける
n_11 <- 13     #がん患者25人を分割
n_12 <- 12
n_21 <- 4      #健康な人9人を分割
n_22 <- 5

n_01 <- n_11 + n_21    #がん患者13人＋健康な人4人
n_02 <- n_12 + n_22    #がん患者12人＋健康な人5人

x_11 <- head(x_1,n = n_11)   #がん患者の行列上半分
x_12 <- tail(x_1,n = n_12)   #がん患者の行列下半分
x_21 <- head(x_2,n = n_21)   #健康な人の行列上半分
x_22 <- tail(x_2,n = n_22)   #健康な人の行列下半分

s_1 <- rbind(x_11,x_21)      #両方のデータが入った行列
s_2 <- rbind(x_12,x_22)      #両方のデータが入った行列

#双対標本分散共分散行列を作る
A_1 <- matrix(nrow = n_01,ncol = p)     #s_1の平均値ベクトル
A_2 <- matrix(nrow = n_02,ncol = p)     #s_2の平均値ベクトル
for (j in 1:n_01) {
  for (i in 1:p) {
    A_1[j,i] <- mean(s_1[,i])
  }
}
for (j in 1:n_02) {
  for (i in 1:p) {
    A_2[j,i] <- mean(s_2[,i])
  }
}

SoD_1 <- sqrt(n_01*n_02)^(-1)*(s_1 - A_1) %*% t(s_2 - A_2)    #分割して作った双対標本分散共分散行列


#クロスデータ行列法？かな？
S2_1 <- SoD_1 %*% t(SoD_1)
eigen(S2_1)$value
S2_2 <- t(SoD_1) %*% SoD_1
eigen(S2_2)$value

U <- array(dim = c(17,17,2))      #S2_1,S2_2の固有ベクトルをarray配列に
U[,,1] <- eigen(S2_1)$vector
U[,,2] <- eigen(S2_2)$vector

u <- array(dim = c(17,3,2))     #S2_1,S2_2の固有ベクトルのうち前から3つを取り出す
for (i in 1:2) {
  for (j in 1:3) {
    u[,j,i] <- U[,j,i]
  }
}

######符号をmathematicaと合わせた。
#for (j in 1:3){
#  for (i in 1:2){
#    print(u[,j,i])
#  }
#}

#u[,1,1] <- (-1) * u[,1,1] 
#u[,1,2] <- (-1) * u[,1,2]
#u[,2,2] <- (-1) * u[,2,2]
#u[,3,1] <- (-1) * u[,3,1]
######

s <- c(1:3)                 #符号の調整?
for (j in 1:3) {
  s[j] <- sign(u[,j,1] %*% (s_1 - A_1) %*% t(s_2 - A_2) %*% u[,j,2])
}
for (i in 1:3){
  u[,i,2] <- s[i] * u[,i,2]
}


t <- matrix(nrow = 3,ncol = n_0)     #????
for (j in 1:3){
  for (k in 1:n_0) {
    if (k <= n_01) {
      t[j,k] <- u[,j,1][k] * sqrt(eigen(S2_1)$value[j] * n_01)
    }
    else {
      t[j,k] <- u[,j,2][k-n_01] * sqrt(eigen(S2_1)$value[j] * n_02)
    }
  }
}

###第１主成分軸、第２主成分軸にデータをプロット
R <- matrix(nrow = 3,ncol = 2)     #plotの範囲を確定
for (i in 1:3){
  R[i,1] <- max(t[i,])
  R[i,2] <- min(t[i,])
}

#すべてのデータの散布図
plot(t[1,],t[2,],xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[2,2],R[2,1]))

#各n_ijのデータの散布図
plot(t[1,][1:n_11],t[2,][1:n_11],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[2,2],R[2,1]),
     col=2,pch=1,cex=2,ann=F)

plot(t[1,][(n_11+1):n_01],t[2,][(n_11+1):n_01],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[2,2],R1[2,1]),
     col=4,pch=4,cex=2,ann=F)

plot(t[1,][(n_01+1):(n_01+n_12)],t[2,][(n_01+1):(n_01+n_12)],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[2,2],R[2,1]),
     col=2,pch=1,cex=2,ann=F)

plot(t[1,][(n_01+n_12+1):n_0],t[2,][(n_01+n_12+1):n_0],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[2,2],R[2,1]),
     col=4,pch=4,cex=2,ann=F)

#上の４つのグラフを重ねる
plot(t[1,][1:n_11],t[2,][1:n_11],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[2,2],R[2,1]),
     col=2,pch=1,cex=2,ann=F)
par(new=T)
plot(t[1,][(n_11+1):n_01],t[2,][(n_11+1):n_01],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[2,2],R[2,1]),
     col=4,pch=4,cex=2,ann=F)
par(new=T)
plot(t[1,][(n_01+1):(n_01+n_12)],t[2,][(n_01+1):(n_01+n_12)],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[2,2],R[2,1]),
     col=2,pch=1,cex=2,ann=F)
par(new=T)
plot(t[1,][(n_01+n_12+1):n_0],t[2,][(n_01+n_12+1):n_0],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[2,2],R[2,1]),
     col=4,pch=4,cex=2,ann=F)


###第１主成分軸と第３主成分軸にデータをプロット

#すべてのデータの散布図
plot(t[1,],t[3,],xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[3,2],R[3,1]))

#各n_ijのデータの散布図
plot(t[1,][1:n_11],t[3,][1:n_11],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[3,2],R[3,1]),
     col=2,pch=1,cex=2,ann=F)

plot(t[1,][(n_11+1):n_01],t[3,][(n_11+1):n_01],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[3,2],R[3,1]),
     col=4,pch=4,cex=2,ann=F)

plot(t[1,][(n_01+1):(n_01+n_12)],t[3,][(n_01+1):(n_01+n_12)],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[3,2],R[3,1]),
     col=2,pch=1,cex=2,ann=F)
s
plot(t[1,][(n_01+n_12+1):n_0],t[3,][(n_01+n_12+1):n_0],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[3,2],R[3,1]),
     col=4,pch=4,cex=2,ann=F)

#上の４つのグラフを重ねる
plot(t[1,][1:n_11],t[3,][1:n_11],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[3,2],R[3,1]),
     col=2,pch=1,cex=2,ann=F)
par(new=T)
plot(t[1,][(n_11+1):n_01],t[3,][(n_11+1):n_01],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[3,2],R[3,1]),
     col=4,pch=4,cex=2,ann=F)
par(new=T)
plot(t[1,][(n_01+1):(n_01+n_12)],t[3,][(n_01+1):(n_01+n_12)],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[3,2],R[3,1]),
     col=2,pch=1,cex=2,ann=F)
par(new=T)
plot(t[1,][(n_01+n_12+1):n_0],t[3,][(n_01+n_12+1):n_0],
     xlim = 1.5*c(R[1,2],R[1,1]),ylim = 2*c(R[3,2],R[3,1]),
     col=4,pch=4,cex=2,ann=F)

