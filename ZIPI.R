# 载入必要的库
library(Hmisc)
library(reshape2)
library(igraph)
library(ggplot2)

# 读取四组微生物数据
bac <- read.csv('bacteria.csv', row.name = 1, check.names = FALSE)
arc <- read.csv('archaea.csv', row.name = 1, check.names = FALSE)
euk <- read.csv('eukaryote.csv', row.name = 1, check.names = FALSE)
vir <- read.csv('virus.csv', row.name = 1, check.names = FALSE)

# 假设这些数据是丰度数据，进行标准化
standardize_data <- function(df) {
  return(scale(df, center = TRUE, scale = TRUE))
}

# 对四组微生物数据进行标准化
bacteria <- standardize_data(bac)
archaea <- standardize_data(arc)
eukaryote <- standardize_data(euk)
virus <- standardize_data(vir)

# 转置数据
bacteria <- data.frame(t(bacteria), check.names = FALSE)
archaea <- data.frame(t(archaea), check.names = FALSE)
eukaryote <- data.frame(t(eukaryote), check.names = FALSE)
virus <- data.frame(t(virus), check.names = FALSE)

# 合并所有微生物数据为一个数据框
combined_data <- cbind(bacteria, archaea, eukaryote, virus)

# 计算相关性
rcorr_result <- rcorr(as.matrix(combined_data), type = "spearman")
r <- rcorr_result$r
p <- rcorr_result$P

# 设置阈值
threshold_r <- 0.7
threshold_p <- 0.05

# 筛选相关性矩阵
r[abs(r) < threshold_r] <- 0
p <- p.adjust(p, method = 'BH')
p[p >= threshold_p] <- -1
p[p < threshold_p & p >= 0] <- 1
p[p == -1] <- 0
z <- r * p

# 将邻接矩阵转换为 igraph 网络对象
g <- graph_from_adjacency_matrix(as.matrix(z), weighted = TRUE, mode = 'undirected')

# 自相关处理
g <- simplify(g)

# 删除孤立节点
g <- delete_vertices(g, names(degree(g)[degree(g) == 0]))

# 设置边权重
E(g)$weight <- abs(E(g)$weight)

# 添加微生物分类属性（假设你有一个合适的分类数据框）
tax <- read.csv("net_tax.csv", row.name = 1)
tax <- tax[as.character(V(g)$name), ]
unmatched_names <- setdiff(as.character(V(g)$name), rownames(tax))
if (length(unmatched_names) > 0) {
  cat("以下节点在 tax 中找不到对应的分类信息:\n")
  print(unmatched_names)
} else {
  tax <- tax[match(as.character(V(g)$name), rownames(tax)), ]
  V(g)$Domain <- tax$Domain
  V(g)$Kingdom <- tax$Kingdom
  V(g)$Phylum <- tax$Phylum
  V(g)$Class <- tax$Class
  V(g)$Order <- tax$Order
  V(g)$Family <- tax$Family
  V(g)$Genus <- tax$Genus
  V(g)$Species <- tax$Species
}

# 绘制网络图
plot(g, vertex.label = V(g)$name, edge.width = E(g)$weight, vertex.size = 5)

# 设置随机种子以确保结果的一致性
set.seed(123)

# g是网络图对象
communities <- cluster_fast_greedy(g)

comm_membership <- membership(communities)

# 计算每个节点的度
node_degree <- degree(g)

# 计算每个模块的平均度和标准差
module_mean <- tapply(node_degree, comm_membership, mean)
module_sd <- tapply(node_degree, comm_membership, sd)

# 计算Zi，并处理标准差为0的情况
Zi <- ifelse(module_sd[comm_membership] != 0,
             (node_degree - module_mean[comm_membership]) / module_sd[comm_membership],
             0)

# 创建一个矩阵来存储每个节点与每个模块的连接数
module_connections <- matrix(0, nrow=length(V(g)), ncol=max(comm_membership))

# 计算每个节点与每个模块的连接数
for (v in V(g)) {
  neighbors_of_v <- neighbors(g, v)
  modules_of_neighbors <- comm_membership[neighbors_of_v]
  module_connections[v, ] <- table(factor(modules_of_neighbors, levels=1:max(comm_membership)))
}

# 计算Pi
Pi <- 1 - rowSums((module_connections / node_degree)^2)

# 将结果整合到一个数据框
zi_pi_metrics <- data.frame(
  name = V(g)$name,
  Zi = Zi,
  Pi = Pi
)

# 节点分类
zi_pi_metrics$type <- ifelse(zi_pi_metrics$Zi > 2.5 & zi_pi_metrics$Pi < 0.62, "Module hubs",
                             ifelse(zi_pi_metrics$Zi < 2.5 & zi_pi_metrics$Pi > 0.62, "Connectors",
                                    ifelse(zi_pi_metrics$Zi > 2.5 & zi_pi_metrics$Pi > 0.62, "Network hubs",
                                           "Peripherals")))
# 添加 OTU 相关的物种分类信息
zi_pi_metrics <- merge(zi_pi_metrics, tax, by.x = "name", by.y = "row.names", all.x = TRUE)

# 导出CSV
write.csv(zi_pi_metrics, "zi_pi_metrics.csv", row.names = FALSE)

p <- ggplot(zi_pi_metrics, aes(x = Pi, y = Zi, color = type)) +
  geom_point(alpha = 0.7, size = 3) +
  theme_minimal() +
  labs(
    x = "Among-module conectivities (Pi)",
    y = "Within-module conectivities (Zi)",
    color = "Node Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  geom_vline(xintercept = 0.62, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 2.5, linetype = "dashed", alpha = 0.5)

p
