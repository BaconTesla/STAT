library(ggplot2)
# make a fake data matrix
data.matrix <- matrix(nrow = 100, ncol = 10)
colnames(data.matrix) <- c(
  paste('wt', 1:5, sep = ''),
  paste('ko', 1:5, sep = '')
)
rownames(data.matrix) <- paste('gene', 1:100, sep = '')
for (i in 1:100) {
  wt.values <- rpois(5, lambda = sample(x=10:100, size=1))
  ko.values <- rpois(5, lambda = sample(x=10:100, size=1))
  
  data.matrix[i,] <- c(wt.values, ko.values)
}
head(data.matrix)
dim(data.matrix)

pca <- prcomp(t(data.matrix), scale = TRUE)

plot(pca$x[,1], pca$x[,2],main = 'pca')


pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main = 'Scree Plot', xlab = 'Principal Component', ylab = 'Percent Variation')

# 使用ggplot2绘图
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data


ggplot(data = pca.data, aes(x = X, y = Y, label = Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep='')) +
  ylab(paste("PC2- ", pca.var.per[2], "%", sep='')) +
  theme_bw() +
  ggtitle("My PCA Graph")


loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores)
gene_score_ranked <- sort(gene_scores, decreasing = TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes


# 使用svd也可以完成该任务
# svd 奇异值分解

# v就是prcomp函数返回的rotation，一个特征向量矩阵
# u类似与prcomp函数返回的x，不过被压缩为单位向量，每一列为一个PC，每一行是一个样本。可以通过乘d来将其展开
# d类似于prcomp函数返回的sdev，不过不是按照无偏估计量缩放的，sdev=sqrt(ss(fit)/(n-1)) d=sqrt(ss(fit))

svd.stuff <- svd(scale(t(data.matrix), center = TRUE))

svd.data <- data.frame(Sample = colnames(data.matrix),
                       X=(svd.stuff$u[,1] * svd.stuff$d[1]),
                       Y=(svd.stuff$u[,2] * svd.stuff$d[2])
)
svd.data

svd.df <- ncol(data.matrix) - 1
svd.var <- svd.stuff$d^2 / svd.df
svd.var.per <- round(svd.var/sum(svd.var)*100, 1)

ggplot(svd.data, aes(x = X, y = Y, label = Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", svd.var.per[1], '%', sep='')) +
  ylab(paste("PC2 - ", svd.var.per[2], '%', sep='')) +
  theme_bw() +
  ggtitle("svd(scale(t(data.matrix), center = TRUE))")



# use eigen do the same thing again!
cov.mat <- cov(scale(t(data.matrix), center = TRUE))
dim(cov.mat)

eigen.stuff <- eigen(cov.mat, symmetric = TRUE)
dim(eigen.stuff$vectors)
head(eigen.stuff$vectors[,1:2])

eigen.pcs <- t(t(eigen.stuff$vectors) %*% t(scale(t(data.matrix), center = TRUE)))
eigen.pcs[,1:2]

eigen.data <- data.frame(Sample = rownames(eigen.pcs),
                         X = (-1 * eigen.pcs[,1]),
                         Y = (eigen.pcs[,2]))
eigen.pcs

eigen.var.per <- round(eigen.stuff$values/sum(eigen.stuff$values)*100,1)

ggplot(data = eigen.data, aes(x = X, y = Y, label = Sample)) +
  geom_text() + 
  xlab(paste("PC1 - ", eigen.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", eigen.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("eigen on cov(t(data.matrix))")












