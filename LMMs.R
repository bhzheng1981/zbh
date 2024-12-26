# 加载必要的包
library(lme4)
library(lmerTest)
library(ggplot2)
library(car)    # 用于Levene's检验
library(gridExtra)

# 读取CSV文件
dataset <- read.csv("env.csv")

# 查看数据结构
str(dataset)# 对 Algal_density 进行 log10 转换
dataset$Algal_density <- log10(dataset$Algal_density)

# 数据预处理（标准化所有数值变量，不包括 Algal_density）
vars_to_scale <- c("WT", "DO", "pH", "Cond", "CODMn", "TN", "TDN", "NO3_N", "NH4_N", "TP", "TDP", "PO4_P", "PP")
dataset[vars_to_scale] <- scale(dataset[vars_to_scale])

# 检查数据结构
str(dataset)

# 正态性检验（Shapiro-Wilk检验）
shapiro_wt <- shapiro.test(dataset$WT)
shapiro_po4 <- shapiro.test(dataset$PO4_P)
shapiro_algal <- shapiro.test(dataset$Algal_density)

cat("Shapiro-Wilk Test p-values:\n")
cat("WT: ", shapiro_wt$p.value, "\n")
cat("PO4_P: ", shapiro_po4$p.value, "\n")
cat("Algal_density: ", shapiro_algal$p.value, "\n")

# 绘制Q-Q图检查正态性
par(mfrow = c(1, 3))  # 设置三幅图同时展示
qqnorm(dataset$Algal_density); qqline(dataset$Algal_density, col = "blue", main = "Q-Q Plot of Algal Density")
qqnorm(dataset$WT); qqline(dataset$WT, col = "blue", main = "Q-Q Plot of WT")
qqnorm(dataset$PO4_P); qqline(dataset$PO4_P, col = "blue", main = "Q-Q Plot of PO4_P")
# 确保 Sample 被转换为因子类型
dataset$Sample <- as.factor(dataset$Sample)
# 方差齐性检验（Levene检验）
leveneTest_result <- leveneTest(Algal_density ~ Sample, data = dataset)
cat("Levene's Test p-value: ", leveneTest_result$p.value, "\n")
leveneTest_result
# 构建混合效应模型
model1 <- lmer(Algal_density ~ WT + (1 | Sample), data = dataset)
model2 <- lmer(Algal_density ~ PO4_P + (1 | Sample), data = dataset)
model3 <- lmer(Algal_density ~ WT + PO4_P + (1 | Sample), data = dataset)
model4 <- lmer(Algal_density ~ WT * PO4_P + (1 | Sample), data = dataset)

# 模型比较
anova(model1, model2, model3, model4)

# 使用模式4作为最终模型
final_model <- model4

# 检查模型结果
summary(final_model)

# 分析固定效应（Type III ANOVA）
anova(final_model)

# 检查模型残差的正态性
par(mfrow = c(1, 1))  # 恢复单图模式
qqnorm(resid(final_model))
qqline(resid(final_model), col = "red", main = "Q-Q Plot of Model Residuals")
shapiro.test(resid(final_model))
# 创建一个数据框用于绘图（生成一系列 WT 和 PO4_P 的组合）
new_data <- expand.grid(
  WT = seq(min(dataset$WT), max(dataset$WT), length.out = 100),
  PO4_P = c(min(dataset$PO4_P), median(dataset$PO4_P), max(dataset$PO4_P))  # 取 PO4_P 的最小值、中位数和最大值
)

# 基于模型生成预测值
new_data$Algal_density_pred <- predict(final_model, newdata = new_data, re.form = NA)

# 绘制交互效应图
ggplot(new_data, aes(x = WT, y = Algal_density_pred, color = factor(PO4_P))) +
  geom_line(linewidth = 1) +  # 将 size 改为 linewidth
  labs(
    title = "Interaction Effect of WT and PO4_P on Algal Density",
    x = "WT (Water Temperature)",
    y = "Predicted Algal Density (log10)",
    color = "PO4_P Level"
  ) +
  theme_minimal()







