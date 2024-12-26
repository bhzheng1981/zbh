library(BiocManager)
library(plspm)
#install.packages("plsdepot")
library(plsdepot)
data = read.csv("sem.csv", header = TRUE, row.names = 1)
data_scaled <- scale(data)
# 对 Alga_density 变量进行对数变换

########################################################
# 构建路径模型矩阵
WT = c(0, 0, 0, 0, 0, 0, 0)
PO4_P= c(0, 0, 0, 0, 0, 0, 0)
Archaea = c(1, 1, 0, 0, 0, 0, 0)
Eukaryota = c(1, 1, 1, 0, 0, 0, 0)
Viruses = c(1, 1, 1, 1, 0, 0, 0)
Bacteria = c(1, 1, 1, 1, 1, 0, 0)
Alga_density = c(1, 1, 1, 1, 1, 1,  0)
# matrix (by row binding)
edu_path = rbind(WT, PO4_P, Archaea, Eukaryota, Viruses, Bacteria, Alga_density)
colnames(edu_path) = rownames(edu_path)
# plot the inner matrix
innerplot(edu_path, box.size = 0.1)
#设置隐变量对应的显变量数据和模型类型
# outer model
edu_blocks = list(1, 2, 3, 4, 5, 6, 7)
# modes (reflective blocks)
edu_modes = rep("A", 7)

# 确认路径矩阵的行名和列名
colnames(edu_path) = rownames(edu_path)
#运行模型
# apply plspm
edu_pls1 = plspm(data_scaled, edu_path, edu_blocks, modes = edu_modes)
# print edu_pls1
edu_pls1
#查看模型的全部结果
#summary()函数湖展示PLS_PM全部结果
summary(edu_pls1)
plot(edu_pls1,box.size = 0.1, box.col = "gray95", lcol = "gray")
#显变量之间应该一致的
edu_pls1$unidim
#下面出图看看
plot(edu_pls1, what = "loadings")
#check unidimensionality：小于0.7的就代表这些变量中存在问题,用下图进行检查
library(ggplot2)
# barchart of loadings
ggplot(data = edu_pls1$outer_model, aes(x = name, y = loading, fill = block)) +  
  geom_bar(stat = "identity" , position = "dodge") +  
  # threshold line (to peek acceptable loadings above 0.7)
  geom_hline(yintercept = 0.7, color = "gray50" ) +
  # add title
  ggtitle("Barchart of Loadings") +
  # rotate x-axis names
  theme(axis.text.x = element_text(angle = 90))
################检测显变量是否合适，观察对角线值是否都大于同行的其他值edu_pls2$crossloadingspath=edu_pls2$path_coefs   path1=abs(path)
#路径效应指数
edu_pls1$path_coefs
#显变量对隐变量的解释，loading
edu_pls1$outer_model
#这里是R方
aa= summary(edu_pls1)
aa$inner_summary
#显著性
edu_pls1$inner_model
#提取模型拟合度
edu_pls1$gof
#绘制路径图
innerplot(edu_pls1)
#提取影响
edu_pls1$effects
#查看模型全部结果
summary(edu_pls1)
###############展示不同隐变量之间的贡献关系
#出图，添加路径尺度大小
plot(edu_pls1, arr.pos = 0.35,box.size = 0.1, box.col = "gray95", lcol = "gray")
Paths =edu_pls1$path_coefs
arrow_lwd =10* round(Paths, 2)
plot(edu_pls1, arr.pos = 0.35, arr.lwd = arrow_lwd)
plot(edu_pls1,box.size = 0.1, box.col = "gray95", lcol = "gray",arr.pos = 0.35,arr.lwd = arrow_lwd)
#效应可视化
good_rows = c(2:21)#print(edu_pls1$effects)根据这个代码修改前面的行号
#
path_effs = as.matrix(edu_pls1$effects[good_rows, 2:3])
rownames(path_effs) = edu_pls1$effects[good_rows, 1]
# setting margin size
op = par(mar = c(8, 3, 1, 0.5))
# barplots of total effects (direct + indirect)
barplot(t(path_effs), border = NA, col = c("#9E9AC8", "#DADAEB"),
        las = 2, cex.names = 0.8, cex.axis = 0.8,
        legend = c("Direct", "Indirect"),
        args.legend = list(x = "top", ncol = 2, border = NA,
                           bty = "n", title = "Effects"))
# resetting default margins
par(op)
#############提取模型的R和Q，计算模型拟合度指数GoF
# 观测值
observed <- data$Alga_density

# 预测值
predicted <- edu_pls1$predictions[, "Alga_density"]

# 打印 inner_summary
cat("Inner Summary:\n")
print(edu_pls1$inner_summary)

# 提取 R²
r_squared <- edu_pls1$inner_summary[edu_pls1$inner_summary$Type == "Endogenous" & 
                                      rownames(edu_pls1$inner_summary) == "Alga_density", "R2"]

# 打印 R²
cat("提取的 R²:", r_squared, "\n")

# 手动计算 R²
ss_total <- sum((observed - mean(observed))^2)
ss_residual <- sum((observed - predicted)^2)
r_squared_manual <- 1 - (ss_residual / ss_total)

# 打印手动计算的 R²
cat("手动计算的 R²:", r_squared_manual, "\n")

# 计算 Q²
q_squared <- 1 - (sum((observed - predicted)^2) / sum((observed - mean(observed))^2))

# 打印 Q²
cat("Q²:", q_squared, "\n")

# 计算 GoF
gof <- sqrt(r_squared * q_squared)

# 输出 GoF
if (!is.na(gof)) {
  print(gof)
} else {
  cat("GoF 计算出错，返回值为 NA。\n")
}



##########boot检验后，我们对于各种指标就会得到误差无置信区间。
#包括显变量对隐变量的影响指标weigt和loading命令调取
#：$boot$weigts和$boot$loadings，path路径$boot$paths，R2$boot$rsq,
#以及隐变量之间的影响$boot$total.efs，每个值都有标准误和在95%区间的最低和最高值。
#文章上使用的人还不多。
# boot检验
#foot_val =  plspm(data_scaled, edu_path, edu_blocks, modes = edu_modes,
#                  boot.val = TRUE, br = 200)
#foot_val$boot



