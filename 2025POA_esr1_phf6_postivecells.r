hypomap = readRDS('f:/hypoMap.rds')
head(hypomap@meta.data)

colnames(hypomap@meta.data)
table(hypomap@meta.data$Dataset)


Moffit10x = subset(hypomap,subset = Dataset == 'Moffit10x')
table(Moffit10x$Sample_ID)
table(Moffit10x$Sex)


POA_2postive <- subset(Moffit10x, subset = Esr1 > 0)
POA_2postive <- subset(POA_2postive, subset = Phf6 > 0)


setwd('e:/Personal_data/WYX/huangju/')
merged_seurat = readRDS('merged_seurat.rds')
phf6 = subset(merged_seurat,subset = Phf6 > 0)
phf6 <- subset(phf6, subset = Esr1 > 0)

POA_2postive = merge(phf6,POA_2postive)

POA_2postive = NormalizeData(POA_2postive)
POA_2postive = FindVariableFeatures(POA_2postive,selection.method = "vst", nfeatures = 2000)
POA_2postive = ScaleData(POA_2postive) 
POA_2postive = RunPCA(POA_2postive, npcs = 30, verbose = FALSE)

POA_2postive <- IntegrateLayers(
  object = POA_2postive, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

POA_2postive <- RunUMAP(POA_2postive, dims = 1:10, reduction = 'integrated.cca')
POA_2postive <- FindNeighbors(POA_2postive, dims = 1:10, reduction = 'integrated.cca')
POA_2postive <- FindClusters(POA_2postive, resolution = 0.08, reduction = 'integrated.cca')

POA_2postive$Sex[is.na(POA_2postive$Sex)] = 'M'
POA_2postive@meta.data[POA_2postive$Condition %in% c('FR','FNR'),]$Sex = 'F'
DimPlot(POA_2postive,group.by = 'seurat_clusters',split.by = 'Sex')

FeaturePlot(POA_2postive,
            features = c("Syt1"),
            max.cutoff = 2.5,split.by = 'Sex',reduction = 'umap')

VlnPlot(POA_2postive,features = c('Esr1',"Phf6",'Syt1','Slc17a6','Slc32a1','Gad1'),group.by = 'RNA_snn_res.0.08',stack = T,pt.size = 0)

POA_2postive$neutype = 'Glu'
POA_2postive@meta.data[POA_2postive@meta.data$RNA_snn_res.0.08 %in%c(0,1),]$neutype = 'Gaba'


pdf('./Phf6_Slc17a6_Slc32a1_Gad1_featureplot.pdf',height = 10,width = 6)
FeaturePlot(POA_2postive,
            features = c("Phf6",'Slc17a6','Slc32a1','Gad1'),order = T,
            split.by = 'Sex',reduction = 'umap')
dev.off()


DimPlot(POA_2postive,group.by = 'neutype',split.by = 'Sex')




length(POA_2postive@meta.data[POA_2postive$neutype == 'Glu'& POA_2postive$Sex == 'M',]$neutype)
length(POA_2postive@meta.data[POA_2postive$neutype == 'Gaba'& POA_2postive$Sex == 'M',]$neutype)

length(POA_2postive@meta.data[POA_2postive$neutype == 'Glu'& POA_2postive$Sex == 'F',]$neutype)
length(POA_2postive@meta.data[POA_2postive$neutype == 'Gaba'& POA_2postive$Sex == 'F',]$neutype)

counts <- c(Glu = 964, Gaba = 2269)
pct   <- round(100 * counts / sum(counts), 1)
labs  <- paste0(names(counts), ": ", counts, " (", pct, "%)")
data <- data.frame(
  category = names(counts),
  count = counts,
  pct = pct,
  label = labs
)
p1 =ggplot(data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +  # 去除背景和坐标轴
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  ggtitle("Female")


counts <- c(Glu = 858, Gaba = 2551)
pct   <- round(100 * counts / sum(counts), 1)
labs  <- paste0(names(counts), ": ", counts, " (", pct, "%)")
data <- data.frame(
  category = names(counts),
  count = counts,
  pct = pct,
  label = labs
)
p2 = ggplot(data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +  # 去除背景和坐标轴
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  ggtitle("Male")

p1+p2



table(POA_2postive@meta.data$Sample_ID)

POA_2postive@meta.data$Sample_ID[is.na(POA_2postive@meta.data$Sample_ID)] = 'GSE183093'

POA_2postive$sexsample = paste(POA_2postive$Sex,POA_2postive$Sample_ID,sep = '_')

POA_2postive@meta.data$Condition[is.na(POA_2postive@meta.data$Condition)] = '-'

POA_2postive@meta.data[POA_2postive@meta.data$Condition=='FNR',]$sexsample = 'FNR_GSE183093'
table(POA_2postive$sexsample)
table(POA_2postive$neutype)



# 计算 sexsample 和 neutype 的交叉表
cross_tab <- table(POA_2postive$sexsample, POA_2postive$neutype)
# 查看交叉表
print(cross_tab)
# 计算每个 sexsample 中 neutype 的比例
prop_cross_tab <- prop.table(cross_tab, 1) * 100  # 按行（sexsample）计算比例
print(prop_cross_tab)



# 删除极值 (Glu 为 0 或 100 的行)
filtered_data <- as.data.frame(prop_cross_tab)
filtered_data <- filtered_data[!filtered_data$Freq %in% c(0,100), ]

# 将数据转换为长格式，以便 ggplot 使用
library(tidyr)

# 只保留 Glu 的数据
glu_data <- filtered_data[filtered_data$Var2 == "Glu", ]

# 进行 t-test 检验，假设你想要比较 F 和 M 两组
# 提取 F 和 M 组数据
f_data <- glu_data[grep("^F", glu_data$Var1), ]
m_data <- glu_data[grep("^M", glu_data$Var1), ]

# 进行 t-test 检验
t_test_result <- t.test(f_data$Freq, m_data$Freq)
t_test_result

# 查看 p-value 和显著性
p_value <- t_test_result$p.value
significance <- ifelse(p_value < 0.05, "*", "ns")  # 如果 p < 0.05，标注为 *，否则为 ns

glu_data$Sex = c(rep('F',4),rep('M',3))
# 绘制柱状图
library(ggplot2)
# 假设 glu_data 中有 Sex (性别组) 和 Freq (Glu 比例) 数据
ggplot(glu_data, aes(x = Sex, y = Freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Freq, 1)), position = position_dodge(width = 0.8), vjust = -0.5) +  # 添加百分比标签
  geom_point(aes(x = Sex, y = Freq), color = "black", size = 3, position = position_dodge(width = 0.8)) +  # 添加散点
  labs(x = "Sex Sample", y = "Glu Proportion (%)", title = "Glu Proportions by Sex Sample") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) +  # 设置 y 轴比例从 0 到 50%
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # x 轴标签旋转 45 度
  # 根据 t-test 显著性添加星号
  annotate("text", x = 2, y = max(glu_data$Freq) + 5, label = significance, size = 6, color = "red")


