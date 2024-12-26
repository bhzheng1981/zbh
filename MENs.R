#这个代码要求四组数据的数据量一样
library(Hmisc)
library(reshape2)
library(igraph)
# 读取四组微生物数据
bac <- read.csv('bacteria.csv', row.name = 1, check.names = FALSE)  # 细菌丰度表
arc <- read.csv('archaea.csv', row.name = 1, check.names = FALSE)    # 古菌丰度表
euk <- read.csv('eukaryote.csv', row.name = 1, check.names = FALSE) # 真核生物丰度表
vir <- read.csv('virus.csv', row.name = 1, check.names = FALSE)        # 病毒丰度表
# 假设 bacteria 是原始数据框，包含绝对丰度
standardize_data <- function(df) {
  # Z-score 标准化
  return(scale(df, center = TRUE, scale = TRUE))
}

# 对四组微生物数据进行标准化
bacteria<- standardize_data(bac)
archaea <- standardize_data(arc)
eukaryote <- standardize_data(euk)
virus <- standardize_data(vir)

# 转置数据
bacteria <- data.frame(t(bacteria), check.names = FALSE)
archaea <- data.frame(t(archaea), check.names = FALSE)
eukaryote <- data.frame(t(eukaryote), check.names = FALSE)
virus <- data.frame(t(virus), check.names = FALSE)

# 计算各组之间的相关性（例如：Spearman相关系数）
library(Hmisc)  # 确保加载Hmisc包以使用rcorr函数
bacteria_archaea_corr <- rcorr(as.matrix(bacteria), as.matrix(archaea), type = 'spearman')
bacteria_eukaryote_corr <- rcorr(as.matrix(bacteria), as.matrix(eukaryote), type = 'spearman')
bacteria_virus_corr <- rcorr(as.matrix(bacteria), as.matrix(virus), type = 'spearman')
archaea_eukaryote_corr <- rcorr(as.matrix(archaea), as.matrix(eukaryote), type = 'spearman')
archaea_virus_corr <- rcorr(as.matrix(archaea), as.matrix(virus), type = 'spearman')
eukaryote_virus_corr <- rcorr(as.matrix(eukaryote), as.matrix(virus), type = 'spearman')

# 细菌-古菌
# 提取相关系数矩阵，保留细菌和古菌的相关性r 和显著性 p 值矩阵
r_bacteria_archaea <- bacteria_archaea_corr$r[colnames(bacteria), colnames(archaea)]
p_bacteria_archaea <- bacteria_archaea_corr$P[colnames(bacteria), colnames(archaea)]


# 细菌-真核生物
r_bacteria_eukaryote <- bacteria_eukaryote_corr$r[colnames(bacteria), colnames(eukaryote)]
p_bacteria_eukaryote <- bacteria_eukaryote_corr$P[colnames(bacteria), colnames(eukaryote)]

# 细菌-病毒
r_bacteria_virus <- bacteria_virus_corr$r[colnames(bacteria), colnames(virus)]
p_bacteria_virus <- bacteria_virus_corr$P[colnames(bacteria), colnames(virus)]

# 古菌-真核生物
r_archaea_eukaryote <- archaea_eukaryote_corr$r[colnames(archaea), colnames(eukaryote)]
p_archaea_eukaryote <- archaea_eukaryote_corr$P[colnames(archaea), colnames(eukaryote)]

# 古菌-病毒
r_archaea_virus <- archaea_virus_corr$r[colnames(archaea), colnames(virus)]
p_archaea_virus <- archaea_virus_corr$P[colnames(archaea), colnames(virus)]

# 真核生物-病毒
r_eukaryote_virus <- eukaryote_virus_corr$r[colnames(eukaryote), colnames(virus)]
p_eukaryote_virus <- eukaryote_virus_corr$P[colnames(eukaryote), colnames(virus)]

# 阈值筛选：例如 r>=0.7 且 p<0.05
threshold_r <- 0.7
threshold_p <- 0.05

# 筛选细菌-古菌
r_bacteria_archaea[abs(r_bacteria_archaea) < threshold_r] <- 0
p_bacteria_archaea <- p.adjust(p_bacteria_archaea, method = 'BH')    # 使用 BH 法校正 p 值
p_bacteria_archaea[p_bacteria_archaea >= threshold_p] <- -1
p_bacteria_archaea[p_bacteria_archaea < threshold_p & p_bacteria_archaea >= 0] <- 1
p_bacteria_archaea[p_bacteria_archaea == -1] <- 0
z_bacteria_archaea <- r_bacteria_archaea * p_bacteria_archaea

# 筛选细菌-真核生物
r_bacteria_eukaryote[abs(r_bacteria_eukaryote) < threshold_r] <- 0
p_bacteria_eukaryote <- p.adjust(p_bacteria_eukaryote, method = 'BH') 
p_bacteria_eukaryote[p_bacteria_eukaryote >= threshold_p] <- -1
p_bacteria_eukaryote[p_bacteria_eukaryote < threshold_p & p_bacteria_eukaryote >= 0] <- 1
p_bacteria_eukaryote[p_bacteria_eukaryote == -1] <- 0
z_bacteria_eukaryote <- r_bacteria_eukaryote * p_bacteria_eukaryote

# 筛选细菌-病毒
r_bacteria_virus[abs(r_bacteria_virus) < threshold_r] <- 0
p_bacteria_virus <- p.adjust(p_bacteria_virus, method = 'BH') 
p_bacteria_virus[p_bacteria_virus >= threshold_p] <- -1
p_bacteria_virus[p_bacteria_virus < threshold_p & p_bacteria_virus >= 0] <- 1
p_bacteria_virus[p_bacteria_virus == -1] <- 0
z_bacteria_virus <- r_bacteria_virus * p_bacteria_virus

# 筛选古菌-真核生物
r_archaea_eukaryote[abs(r_archaea_eukaryote) < threshold_r] <- 0
p_archaea_eukaryote <- p.adjust(p_archaea_eukaryote, method = 'BH') 
p_archaea_eukaryote[p_archaea_eukaryote >= threshold_p] <- -1
p_archaea_eukaryote[p_archaea_eukaryote < threshold_p & p_archaea_eukaryote >= 0] <- 1
p_archaea_eukaryote[p_archaea_eukaryote == -1] <- 0
z_archaea_eukaryote <- r_archaea_eukaryote * p_archaea_eukaryote

# 筛选古菌-病毒
r_archaea_virus[abs(r_archaea_virus) < threshold_r] <- 0
p_archaea_virus <- p.adjust(p_archaea_virus, method = 'BH') 
p_archaea_virus[p_archaea_virus >= threshold_p] <- -1
p_archaea_virus[p_archaea_virus < threshold_p & p_archaea_virus >= 0] <- 1
p_archaea_virus[p_archaea_virus == -1] <- 0
z_archaea_virus <- r_archaea_virus * p_archaea_virus

# 筛选真核生物-病毒
r_eukaryote_virus[abs(r_eukaryote_virus) < threshold_r] <- 0
p_eukaryote_virus <- p.adjust(p_eukaryote_virus, method = 'BH') 
p_eukaryote_virus[p_eukaryote_virus >= threshold_p] <- -1
p_eukaryote_virus[p_eukaryote_virus < threshold_p & p_eukaryote_virus >= 0] <- 1
p_eukaryote_virus[p_eukaryote_virus == -1] <- 0
z_eukaryote_virus <- r_eukaryote_virus * p_eukaryote_virus

# 合并四组微生物之间的相关系数矩阵为一个大矩阵
len <- ncol(bacteria) + ncol(archaea) + ncol(eukaryote) + ncol(virus)
z <- data.frame(matrix(data=NA, nrow = len, ncol = len, byrow = FALSE, dimnames = NULL))

rownames(z) <- c(colnames(bacteria), colnames(archaea), colnames(eukaryote), colnames(virus))
colnames(z) <- c(colnames(bacteria), colnames(archaea), colnames(eukaryote), colnames(virus))

# 对称地填充相关系数矩阵
z[colnames(bacteria), colnames(archaea)] <- z_bacteria_archaea
z[colnames(archaea), colnames(bacteria)] <- t(z_bacteria_archaea)

z[colnames(bacteria), colnames(eukaryote)] <- z_bacteria_eukaryote
z[colnames(eukaryote), colnames(bacteria)] <- t(z_bacteria_eukaryote)

z[colnames(bacteria), colnames(virus)] <- z_bacteria_virus
z[colnames(virus), colnames(bacteria)] <- t(z_bacteria_virus)

z[colnames(archaea), colnames(eukaryote)] <- z_archaea_eukaryote
z[colnames(eukaryote), colnames(archaea)] <- t(z_archaea_eukaryote)

z[colnames(archaea), colnames(virus)] <- z_archaea_virus
z[colnames(virus), colnames(archaea)] <- t(z_archaea_virus)

z[colnames(eukaryote), colnames(virus)] <- z_eukaryote_virus
z[colnames(virus), colnames(eukaryote)] <- t(z_eukaryote_virus)

z[colnames(bacteria), colnames(bacteria)] <- 0
z[colnames(archaea), colnames(archaea)] <- 0
z[colnames(eukaryote), colnames(eukaryote)] <- 0
z[colnames(virus), colnames(virus)] <- 0

# 将邻接矩阵转换为 igraph 网络对象
g <- graph_from_adjacency_matrix(as.matrix(z), weighted = TRUE, mode = 'undirected')

# 自相关处理
g <- simplify(g)

# 孤立节点的删除
g <- delete_vertices(g, names(degree(g)[degree(g) == 0]))

# 设置边权重为绝对值，并保留原相关系数
E(g)$correlation <- E(g)$weight

# 添加微生物分类属性
tax <- read.csv("net_tax.csv", row.name = 1)
tax <- tax[as.character(V(g)$name), ]
# 检查 tax 中的行名称是否包含所有图中的节点名称
unmatched_names <- setdiff(as.character(V(g)$name), rownames(tax))
if (length(unmatched_names) > 0) {
  cat("以下节点在 tax 中找不到对应的分类信息:\n")
  print(unmatched_names)
} else {
  # 匹配并排序 tax 的行顺序
  tax <- tax[match(as.character(V(g)$name), rownames(tax)), ]
  
  # 添加分类信息
  V(g)$Domain <- tax$Domain
  V(g)$Kingdom <- tax$Kingdom  # 确保列名一致
  V(g)$Phylum <- tax$Phylum
  V(g)$Class <- tax$Class
  V(g)$Order <- tax$Order
  V(g)$Family <- tax$Family
  V(g)$Genus <- tax$Genus
  V(g)$Species <- tax$Species
}
# 检查节点的分类属性是否已成功添加
print(head(V(g)$Domain))
print(head(V(g)$Kingdom))

# ##############统计下在Domain分类学水平下不同类群之间边的正相关比例

# 1. 获取边的数据
edges <- as.data.frame(as_edgelist(g), stringsAsFactors = FALSE)  # 获取边的列表
edges$weight <- E(g)$weight  # 添加边的权重

# 2. 获取节点的Domain信息
node_domains <- data.frame(Name = V(g)$name, Domain = V(g)$Domain, stringsAsFactors = FALSE)

# 3. 合并边和节点Domain信息
edges <- merge(edges, node_domains, by.x = "V1", by.y = "Name")
edges <- merge(edges, node_domains, by.x = "V2", by.y = "Name", suffixes = c("_from", "_to"))

# 4. 统计正相关边
positive_edges <- edges[edges$weight > 0, ]  # 直接根据权重筛选正相关边

# 5. 计算每个Domain之间的正相关边数量
positive_count <- as.data.frame(table(positive_edges$Domain_from, positive_edges$Domain_to))
colnames(positive_count) <- c("Domain_from", "Domain_to", "Positive_Count")

# 6. 统计每个Domain的边的总数
total_edges_count <- as.data.frame(table(edges$Domain_from, edges$Domain_to))
colnames(total_edges_count) <- c("Domain_from", "Domain_to", "Total_Count")

# 7. 合并正相关边数和总边数
combined_table <- merge(positive_count, total_edges_count, by = c("Domain_from", "Domain_to"), all = TRUE)

# 8. 计算正相关比例（如果需要）
combined_table$Positive_Ratio <- combined_table$Positive_Count / combined_table$Total_Count

# 输出结果
print(combined_table)
write.csv(combined_table, "domain_positive_edges.csv", row.names = FALSE)
##################################################
# 1. 获取边的数据
edges <- as.data.frame(as_edgelist(g), stringsAsFactors = FALSE)  # 获取边的列表
edges$weight <- E(g)$weight  # 添加边的权重

# 2. 获取节点的Genus信息
node_Genus <- data.frame(Name = V(g)$name, Genus = V(g)$Genus, stringsAsFactors = FALSE)

# 3. 合并边和节点Genus信息
edges <- merge(edges, node_Genus, by.x = "V1", by.y = "Name")
edges <- merge(edges, node_Genus, by.x = "V2", by.y = "Name", suffixes = c("_from", "_to"))

# 2. 重新筛选与 Genus 为 Cylindrospermopsis 相关的边
cylindrospermopsis_edges <- edges[edges$Genus_from == "Cylindrospermopsis" | edges$Genus_to == "Cylindrospermopsis", ]
print("cylindrospermopsis_edges 数据:")
print(head(cylindrospermopsis_edges))
# 3. 计算与 Cylindrospermopsis 相关的边的总数
total_edges_cylindrospermopsis <- nrow(cylindrospermopsis_edges)
cat("与 Genus 为 Cylindrospermopsis 相关的边的总数:", total_edges_cylindrospermopsis, "\n")

# 4. 统计每个属（Genus）与 Cylindrospermopsis 的正相关边数量
positive_edges_with_genus <- cylindrospermopsis_edges[cylindrospermopsis_edges$weight > 0, ]
positive_edges_with_genus$Other_Genus <- ifelse(positive_edges_with_genus$Genus_from == "Cylindrospermopsis", 
                                                positive_edges_with_genus$Genus_to, 
                                                positive_edges_with_genus$Genus_from)

# 5. 统计与其他属之间的正相关边数量
genus_positive_count <- as.data.frame(table(positive_edges_with_genus$Other_Genus))
colnames(genus_positive_count) <- c("Genus", "Positive_Count")

# 6. 打印正相关边的结果
print("与 Cylindrospermopsis 存在正相关边的属及其正相关边数量:")
print(genus_positive_count)
write.csv(genus_positive_count, " Cylindrospermopsis 存在正相关边的属.csv", row.names = FALSE)
# 7. 统计每个属（Genus）与 Cylindrospermopsis 的负相关边数量
negative_edges_with_genus <- cylindrospermopsis_edges[cylindrospermopsis_edges$weight < 0, ]
negative_edges_with_genus$Other_Genus <- ifelse(negative_edges_with_genus$Genus_from == "Cylindrospermopsis", 
                                                negative_edges_with_genus$Genus_to, 
                                                negative_edges_with_genus$Genus_from)

# 8. 统计与其他属之间的负相关边数量
genus_negative_count <- as.data.frame(table(negative_edges_with_genus$Other_Genus))
colnames(genus_negative_count) <- c("Genus", "Negative_Count")

# 9. 打印负相关边的结果
print("与 Cylindrospermopsis 存在负相关边的属及其负相关边数量:")
print(genus_negative_count)
write.csv(genus_negative_count, " Cylindrospermopsis 存在负相关边的属.csv", row.names = FALSE)
###################################################获取 Cylindrospermopsis种水平的相关关系
# 1. 获取边的数据
edges <- as.data.frame(as_edgelist(g), stringsAsFactors = FALSE)  # 获取边的列表
edges$weight <- E(g)$weight  # 添加边的权重

# 2. 获取节点的 Species 信息
node_Species <- data.frame(Name = V(g)$name, Species = V(g)$Species, stringsAsFactors = FALSE)

# 3. 合并边和节点 Species 信息
edges <- merge(edges, node_Species, by.x = "V1", by.y = "Name")
edges <- merge(edges, node_Species, by.x = "V2", by.y = "Name", suffixes = c("_from", "_to"))

# 4. 定义目标物种列表（全字符匹配）
target_species <- c("Cylindrospermopsis_sp._CR12", "Cylindrospermopsis_curvispora", 
                    "unclassified_g__Cylindrospermopsis", "Cylindrospermopsis_raciborskii")

# 5. 筛选与目标物种相关的边
cylindrospermopsis_edges <- edges[edges$Species_from %in% target_species | edges$Species_to %in% target_species, ]

# 6. 计算与目标物种相关的边的总数
total_edges_cylindrospermopsis <- nrow(cylindrospermopsis_edges)
cat("与目标物种相关的边的总数:", total_edges_cylindrospermopsis, "\n")

# 7. 统计正相关边的数量
positive_edges_with_species <- cylindrospermopsis_edges[cylindrospermopsis_edges$weight > 0, ]
positive_edges_with_species$Other_Species <- ifelse(positive_edges_with_species$Species_from %in% target_species, 
                                                    positive_edges_with_species$Species_to, 
                                                    positive_edges_with_species$Species_from)

# 8. 统计与其他物种的正相关边数量
species_positive_count <- as.data.frame(table(positive_edges_with_species$Other_Species))
colnames(species_positive_count) <- c("Species", "Positive_Count")

# 9. 统计负相关边的数量
negative_edges_with_species <- cylindrospermopsis_edges[cylindrospermopsis_edges$weight < 0, ]
negative_edges_with_species$Other_Species <- ifelse(negative_edges_with_species$Species_from %in% target_species, 
                                                    negative_edges_with_species$Species_to, 
                                                    negative_edges_with_species$Species_from)

# 10. 统计与其他物种的负相关边数量
species_negative_count <- as.data.frame(table(negative_edges_with_species$Other_Species))
colnames(species_negative_count) <- c("Species", "Negative_Count")

# 11. 输出结果到CSV文件
write.csv(species_positive_count, "Species_正相关边的物种.csv", row.names = FALSE)
write.csv(species_negative_count, "Species_负相关边的物种.csv", row.names = FALSE)



##########################################################
# 将网络导出为 GraphML 格式，方便在 Gephi 中可视化
# 设置边权重为绝对值，并保留原相关系数，进行网络属性的计算
E(g)$weight <- abs(E(g)$weight)
write_graph(g, '四组微生物关系.graphml', format = 'graphml')
# 查看网络图
plot(g)

###################################################
# 计算网络属性
# 计算网络属性
nodes_num <- gorder(g)  # 节点数量
edges_num <- gsize(g)    # 边数量
average_degree <- mean(degree(g))  # 平均度
nodes_connectivity <- vertex_connectivity(g)  # 节点连通度
edges_connectivity <- edge_connectivity(g)  # 边连通度

# 使用新函数计算平均路径长度和图密度
average_path_length <- mean_distance(g)  # 平均路径长度
graph_density <- edge_density(g)  # 图密度

# 其他属性计算
graph_diameter <- diameter(g)  # 网络直径
clustering_coefficient <- transitivity(g)  # 聚类系数

# 使用新函数计算介数中心性和度中心性
betweenness_centralization <- centr_betw(g)$centralization  # 介数中心性
degree_centralization <- centr_degree(g)$centralization  # 度中心性

# 模块性
modularity <- modularity(cluster_fast_greedy(g))  # 模块性

# 计算新的网络属性
# eigenvector_centrality特征向量中心性评估节点的影响力，不仅考虑它直接连接的节点数量，还考虑其连接节点的影响力。
eigenvector_centrality <- eigen_centrality(g)$vector
#Closeness Centrality（接近中心性）评估节点的影响力，不仅考虑它直接连接的节点数量，还考虑其连接节点的影响力。
closeness_centrality <- closeness(g)
#Assortativity（同配性）评估网络中节点是否倾向于连接相似类型的节点（例如高度节点是否倾向于连接其他高度节点）。
assortativity <- assortativity_degree(g)
#Global Efficiency（全局效率）衡量网络整体的效率，表示节点间信息传递的有效性。
global_efficiency <- mean(1 / distances(g))
#Edge Betweenness Centrality（边介数中心性）衡量一条边在多少对节点之间的最短路径中起关键作用，常用于识别桥接不同模块的边。
edge_betweenness_centrality <- edge_betweenness(g)
#Degree Distribution（度分布）网络中节点度的分布，可以通过查看网络的度分布图来判断网络是否呈现无尺度特性。
degree_distribution <- degree_distribution(g)

# 小世界性 Small-worldness，需要结合平均路径长度和聚类系数计算。衡量网络是否具有高聚类系数和较短的平均路径长度，类似于小世界网络的特性。
clustering_coeff <- transitivity(g, type = "global")
average_path_len <- mean_distance(g)
random_graph <- sample_gnm(vcount(g), ecount(g), directed = FALSE)
small_worldness <- (clustering_coeff / transitivity(random_graph)) / (average_path_len / mean_distance(random_graph))

# 局部效率Local Efficiency衡量局部子图中的节点间传递信息的有效性。
local_efficiency <- sapply(V(g), function(v) {
  subgraph <- induced_subgraph(g, ego(g, 1, v)[[1]])
  return(mean(1 / distances(subgraph)))
})

# 平均聚类系数Average Clustering Coefficient衡量网络中节点形成团体或三角形结构的倾向性。
avg_clustering_coefficient <- transitivity(g, type = "average")

# 社区检测Community Detection用于识别网络中的社区或模块结构，
communities <- cluster_walktrap(g)
modularity_score <- modularity(communities)

#打印结果
print(paste("Nodes:", nodes_num))
print(paste("Edges:", edges_num))
print(paste("Average Degree:", average_degree))
print(paste("Node Connectivity:", nodes_connectivity))
print(paste("Edge Connectivity:", edges_connectivity))
print(paste("Average Path Length:", average_path_length))
print(paste("Graph Diameter:", graph_diameter))
print(paste("Graph Density:", graph_density))
print(paste("Clustering Coefficient:", clustering_coefficient))
print(paste("Betweenness Centralization:", betweenness_centralization))
print(paste("Degree Centralization:", degree_centralization))
print(paste("Modularity:", modularity_score))
print(paste("Assortativity:", assortativity))
print(paste("Global Efficiency:", global_efficiency))
print(paste("Average Clustering Coefficient:", avg_clustering_coefficient))
print(paste("Small Worldness:", small_worldness))

print(paste("Closeness Centrality:", closeness_centrality))
print(paste("Eigenvector Centrality:", eigenvector_centrality))
print(paste("Edge Betweenness Centrality:", edge_betweenness_centrality))

