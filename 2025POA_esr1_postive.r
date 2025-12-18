# 导入必要的库
library(dplyr)
library(pheatmap)# 加载 Nebulosa 包
library(Nebulosa)
library(Seurat)
library(readr)
library(data.table)
# 读取计数矩阵
rawcount <- fread("./GSE183093/GSE183093_POAcountmatrix.tsv", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
metadata <- fread("./GSE183093/GSE183093_POA_Cell_metadata.tsv.gz", header = TRUE)
rawcount[1:5,1:5]

POA_count <- as.data.frame(rawcount)

POA_count[1:5,1:5]
tail(colnames(POA_count))


rownames(POA_count) = POA_count$barcode
POA_count <- POA_count[,-1]


POA_count[1:5,1:5]

metadata[1:5,]

# # 转置矩阵
POA_count_t <- t(POA_count)

# 创建 Seurat 对象
# meta_list <- pull(metadata, "Barcode")
# metadata <- metadata[, -1]
# rownames(metadata) <- meta_list
POA_seurat <- CreateSeuratObject(counts = POA_count_t, meta.data = as.data.frame(metadata), min.cells = 3, min.features = 200)

# 对 Seurat 对象进行数据归一化
POA_seurat <- POA_seurat %>% NormalizeData()
POA_seurat@meta.data[1:5,]

# 提取 Esr1 大于 0 的子集
POA_subset_Esr1 <- subset(POA_seurat, subset = Esr1 > 0)

# 读取另一个数据集
POA_1 <- Read10X(data.dir = "./GSE113576/")
POA_1 <- CreateSeuratObject(counts = POA_1, project = "POA", min.cells = 3, min.features = 200)

# 对第二个数据集进行通用分析管道
POA_1 <- POA_1 %>% NormalizeData()
POA_subset_1 <- subset(POA_1, subset = Esr1 > 0)

# 提取共同基因
genes1 <- rownames(POA_subset_Esr1)
genes2 <- rownames(POA_subset_1)
common_genes <- intersect(genes1, genes2)
POA_subset_Esr1_c <- subset(POA_subset_Esr1, features = common_genes)
POA_subset_1_c <- subset(POA_subset_1, features = common_genes)

# 合并变换后的数据
merged_seurat <- merge(x = POA_subset_Esr1_c, y = POA_subset_1_c)

# 对合并后的数据进行标准分析
merged_seurat <- merged_seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)


library(future)
# check the current active plan
plan()
# change the current plan to access parallelization
plan("multisession", workers = 32)
plan()

saveRDS(merged_seurat,'merged_seurat.rds')
# qsave(merged_seurat,'merged_seurat.rds')

setwd('e:/Personal_data/WYX/huangju/')
merged_seurat = readRDS('merged_seurat.rds')
# 使用 CCAIntegration 方法进行整合
merged_seurat <- IntegrateLayers(
  object = merged_seurat, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

merged_seurat <- JoinLayers(merged_seurat)

ElbowPlot(merged_seurat,reduction = "pca")

# 运行 UMAP、找邻居、聚类
merged_seurat <- RunUMAP(merged_seurat, dims = 1:20, reduction = "integrated.cca")
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:20, reduction = "integrated.cca")
options(future.globals.maxSize = 1e11)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.13, reduction = "integrated.cca")


merged_seurat@meta.data$seurat_clusters = merged_seurat@meta.data$RNA_snn_res.0.13
merged_seurat@meta.data$seurat_clusters = paste0('c',merged_seurat@meta.data$seurat_clusters)

# merged_seurat@meta.data[merged_seurat@meta.data$seurat_clusters == 'c2',]$seurat_clusters = 'c10'
# merged_seurat@meta.data[merged_seurat@meta.data$seurat_clusters == 'c1',]$seurat_clusters = 'c2'
# merged_seurat@meta.data[merged_seurat@meta.data$seurat_clusters == 'c10',]$seurat_clusters = 'c1'



library(ggplot2)
DefaultAssay(merged_seurat) = 'INTE'
DotPlot(merged_seurat, features = "Phf6", dot.scale = 12,group.by = 'seurat_clusters') +  scale_color_gradientn(colours  = colorRampPalette(c("white",'#F8E4D9','#E9AD93', "#D22D25")) (100))

pdf('D:/PersonalData/WYX/huangju/202503.phf6.dotplot.pdf',width = 3.5,height = 6)
DotPlot(merged_seurat, features = "Phf6", dot.scale = 12,group.by = 'seurat_clusters') +  scale_color_gradientn(colours  = colorRampPalette(c("white",'#F8E4D9','#E9AD93', "#D22D25")) (100))
dev.off()




DotPlot(merged_seurat,group.by = 'seurat_clusters', features = c('Esr1','Drp2','Pde1c','Serpine2','Slc17a6','St18','Bcl11b','Emx2os','Opn5','Zeb2','Otp',"Phf6",'Slc32a1'), dot.scale = 12) +  scale_color_gradientn(colours  = colorRampPalette(c("blue", "red")) (100))+RotatedAxis()
DotPlot(merged_seurat,group.by = 'seurat_clusters', features = c('Esr1','Slc17a6',"Phf6",'Slc32a1'), dot.scale = 12) +  scale_color_gradientn(colours  = colorRampPalette(c("white",'#F8E4D9','#E9AD93', "#D22D25")) (100))+RotatedAxis()




phf6 = subset(merged_seurat,subset = Phf6 > 0)

# 运行 UMAP、找邻居、聚类
phf6 <- RunUMAP(phf6, dims = 1:20, reduction = "integrated.cca")
phf6 <- FindNeighbors(phf6, dims = 1:20, reduction = "integrated.cca")
# options(future.globals.maxSize = 1e11)
phf6 <- FindClusters(phf6, resolution = 0.1, reduction = "integrated.cca")
phf6 <- FindClusters(phf6, resolution = 0.2, reduction = "integrated.cca")

FeaturePlot(
  phf6,
  features = c("Esr1","Phf6",'Slc17a6','Slc32a1'),max.cutoff = 2.5,split.by = 'Condition')

FeaturePlot(
  phf6,
  features = c("Phf6",'Slc17a6','Slc32a1'),max.cutoff = 2.5)

unique(phf6$Condition)
Idents(phf6) = 'Condition'
phf6=subset(phf6, idents =  c("Male", "FR", "FNR"))

phf6@meta.data$sex = 'Female'
phf6@meta.data[phf6@meta.data$Condition=='Male',]$sex = 'Male'

pdf('featureplot.pdf',width = 5,height = 10)
FeaturePlot(
  phf6,
  features = c("Phf6",'Slc17a6','Slc32a1','Gad1'),max.cutoff = 2.5,split.by = 'sex')
dev.off()

DimPlot(phf6,group.by = "RNA_snn_res.0.1",split.by = 'sex')
DimPlot(phf6,group.by = "RNA_snn_res.0.2")
DimPlot(phf6,group.by = "RNA_snn_res.0.13")

phf6$neutype = 'Glu'
phf6@meta.data[phf6@meta.data$RNA_snn_res.0.1 %in%c(0,1,4,5),]$neutype = 'Gaba'

VlnPlot(phf6,features = c("Phf6",'Slc17a6','Slc32a1','Gad1'),group.by = 'RNA_snn_res.0.1',stack = T,pt.size = 0)
DimPlot(phf6,group.by = "neutype")


table(phf6$orig.ident)

length(phf6@meta.data[phf6$neutype == 'Glu'& phf6$sex == 'Male',]$neutype)
length(phf6@meta.data[phf6$neutype == 'Gaba'& phf6$sex == 'Male',]$neutype)

length(phf6@meta.data[phf6$neutype == 'Glu'& phf6$sex == 'Female',]$neutype)
length(phf6@meta.data[phf6$neutype == 'Gaba'& phf6$sex == 'Female',]$neutype)



counts <- c(Glu = 859, Gaba = 2022)

pct   <- round(100 * counts / sum(counts), 1)
labs  <- paste0(names(counts), ": ", counts, " (", pct, "%)")

pie(counts,
    labels = labs,
    clockwise = TRUE,
    col = c("#6BAED6", "#FD8D3C"),
    main = 'Male')




counts <- c(Glu = 991, Gaba = 1976)

pct   <- round(100 * counts / sum(counts), 1)
labs  <- paste0(names(counts), ": ", counts, " (", pct, "%)")

pie(counts,
    labels = labs,
    clockwise = TRUE,
    col = c("#6BAED6", "#FD8D3C"),
    main = 'Female')



idents(phf6) = 'sex'
phf5.f = subset(phf6, idents = 'Female')
phf5.m = subset(phf6, idents = 'Male')

FeaturePlot(
  phf5.f,
  features = c('Esr1',"Phf6",'Slc17a6','Gad1'),max.cutoff = 2.5)

FeaturePlot(
  phf5.m,
  features = c('Esr1',"Phf6",'Slc17a6','Gad1'),max.cutoff = 2.5)

pdf('D:/PersonalData/WYX/huangju/Esr1_Phf6_Slc17a6_Gad1_featureplot_female_exp2.pdf',height = 7.4,width = 9)
FeaturePlot(
  phf5.f,
  features = c('Esr1',"Phf6",'Slc17a6','Gad1'),max.cutoff = 2)
dev.off()

pdf('D:/PersonalData/WYX/huangju/Esr1_Phf6_Slc17a6_Gad1_featureplot_male_exp2.pdf',height = 7.4,width = 9)
FeaturePlot(
  phf5.m,
  features = c('Esr1',"Phf6",'Slc17a6','Gad1'),max.cutoff = 2)
dev.off()


pdf('D:/PersonalData/WYX/huangju/Esr1_Phf6_Slc17a6_Slc32a1_featureplot_female_exp2.pdf',height = 7.4,width = 9)
FeaturePlot(
  phf5.f,
  features = c('Esr1',"Phf6",'Slc17a6','Slc32a1'),max.cutoff = 2)
dev.off()

pdf('D:/PersonalData/WYX/huangju/Esr1_Phf6_Slc17a6_Slc32a1_featureplot_male_exp2.pdf',height = 7.4,width = 9)
FeaturePlot(
  phf5.m,
  features = c('Esr1',"Phf6",'Slc17a6','Slc32a1'),max.cutoff = 2)
dev.off()



pdf('D:/PersonalData/WYX/huangju/Esr1_Phf6_Slc17a6_Gad1_featureplot_female_exp2.pdf',height = 7.4,width = 9)
FeaturePlot(
  phf5.f,
  features = c('Esr1',"Phf6",'Slc17a6','Gad1'),max.cutoff = 2)
dev.off()

pdf('D:/PersonalData/WYX/huangju/Esr1_Phf6_Slc17a6_Gad1_featureplot_male_exp2.pdf',height = 7.4,width = 9)
FeaturePlot(
  phf5.m,
  features = c('Esr1',"Phf6",'Slc17a6','Gad1'),max.cutoff = 2)
dev.off()


pdf('D:/PersonalData/WYX/huangju/Esr1_Phf6_Slc17a6_Slc32a1_featureplot_female_exp2.pdf',height = 7.4,width = 9)
FeaturePlot(
  phf5.f,
  features = c('Esr1',"Phf6",'Slc17a6','Slc32a1'),max.cutoff = 2)
dev.off()

pdf('D:/PersonalData/WYX/huangju/Esr1_Phf6_Slc17a6_Slc32a1_featureplot_male_exp2.pdf',height = 7.4,width = 9)
FeaturePlot(
  phf5.m,
  features = c('Esr1',"Phf6",'Slc17a6','Slc32a1'),max.cutoff = 2)
dev.off()
FeaturePlot(
  phf6,
  features = c('Syt1'),max.cutoff = 2)





setwd('D:/PersonalData/WYX/huangju/')


for (i in c(0.05,0.1,0.2)) {
  phf6 = FindClusters(phf6,resolution = i)
  print(i)
  
  pdf(paste0('phf6_res',i,'.pdf'),width = 4.5,height = 4)
  print(DimPlot(phf6,group.by = paste0('RNA_snn_res.',i)))
  dev.off()
  
  marker= FindAllMarkers(phf6, group.by = paste0('RNA_snn_res.',i),only.pos = T,logfc.threshold = 1)
  marker= marker.res0.1[marker.res0.1$p_val_adj<0.05,]
  writexl::write_xlsx(marker,paste0('phf6_res.',i,'.xlsx'))
}


saveRDS(phf6,file = '202506phf6.rds')
saveRDS(phf6,file = '202506phf6_2.rds')

phf6 <- phf6 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 10, verbose = FALSE)


phf6 <- RunUMAP(phf6, dims = 1:10, reduction = "integrated.cca")
phf6 <- FindNeighbors(phf6, dims = 1:10, reduction = "integrated.cca")


Idents(phf6) = 'sex'
phf5.f = subset(phf6, idents = 'Female')
phf5.m = subset(phf6, idents = 'Male')

pdf('D:/PersonalData/WYX/huangju/Esr1_Phf6_Slc17a6_Gad1_featureplot_female_exp2.pdf',height = 7.4,width = 9)
FeaturePlot(
  phf5.f,
  features = c('Esr1',"Phf6",'Slc17a6','Gad1'),max.cutoff = 2)
dev.off()

pdf('D:/PersonalData/WYX/huangju/Esr1_Phf6_Slc17a6_Gad1_featureplot_male_exp2.pdf',height = 7.4,width = 9)
FeaturePlot(
  phf5.m,
  features = c('Esr1',"Phf6",'Slc17a6','Gad1'),max.cutoff = 2)
dev.off()


pdf('D:/PersonalData/WYX/huangju/Esr1_Phf6_Slc17a6_Slc32a1_featureplot_female_exp2.pdf',height = 7.4,width = 9)
FeaturePlot(
  phf5.f,
  features = c('Esr1',"Phf6",'Slc17a6','Slc32a1'),max.cutoff = 2)
dev.off()

pdf('D:/PersonalData/WYX/huangju/Esr1_Phf6_Slc17a6_Slc32a1_featureplot_male_exp2.pdf',height = 7.4,width = 9)
FeaturePlot(
  phf5.m,
  features = c('Esr1',"Phf6",'Slc17a6','Slc32a1'),max.cutoff = 2)
dev.off()

for (i in c(0.05,0.1,0.2)) {
  phf6 = FindClusters(phf6,resolution = i)
  print(i)
  
  pdf(paste0('phf6_res',i,'.pdf'),width = 4.5,height = 4)
  print(DimPlot(phf6,group.by = paste0('RNA_snn_res.',i)))
  dev.off()
  
  marker= FindAllMarkers(phf6, group.by = paste0('RNA_snn_res.',i),only.pos = T,logfc.threshold = 1)
  marker= marker[marker$p_val_adj<0.05,]
  writexl::write_xlsx(marker,paste0('phf6_res.',i,'.xlsx'))
}


marker= FindAllMarkers(phf6, group.by = paste0('RNA_snn_res.0.1'),only.pos = T,logfc.threshold = 1)
marker= marker[marker$p_val_adj<0.05,]

marker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 2) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3


DotPlot(phf6,c('Esr1','Phf6',top3$gene),group.by = 'RNA_snn_res.0.1')+RotatedAxis()
VlnPlot(phf6,c('Esr1','Phf6',top5$gene),group.by = 'RNA_snn_res.0.1',stack=T)+RotatedAxis()

selectgene = c('Esr1','Phf6', "Ntng1","Nts","Samd3","Adamts2","Emx2os","Bcl11b" ,"Lhx8","Pantr1",
               "Slc17a6","Hs3st2","Ptprb", "Nfix","Moxd1","St18")

VlnPlot(phf6,selectgene,group.by = 'RNA_snn_res.0.1',stack=T)+RotatedAxis()





# 寻找差异基因
POA_markers <- FindAllMarkers(merged_seurat, only.pos = TRUE, min.pct = 0.4, logfc.threshold = 1)

# 对每个细胞簇进行排序
POA_markers_sorted <- POA_markers %>%
  arrange(cluster, p_val_adj)

# 选择每个细胞簇中 avg_log2FC > 1 且 p_val_adj 最小的唯一基因
selected_markers <- POA_markers_sorted %>%
  group_by(cluster) %>%
  filter(avg_log2FC> 1.2) %>%
  slice_min(order_by = p_val_adj) %>%
  distinct(cluster, .keep_all = TRUE)

# 选择每个细胞簇中的 marker 基因
selected_genes <- c(selected_markers$gene)

# 使用 VlnPlot 绘制堆叠小提琴图
# VlnPlot(
#   merged_seurat,
#   features = selected_genes,
#   stack = TRUE,
#   group.by = "seurat_clusters",
#   pt.size = 0,
#   combine = TRUE
# ) +
#   ggtitle(NULL) +  # Remove the title
#   theme(
#     axis.title.x = element_blank(),  # Remove x-axis label
#     axis.text.x = element_blank(),    # Remove x-axis ticks
#     axis.ticks.x = element_blank(),   # Remove x-axis tick marks
#     strip.text.x = element_text(size=8,face = "plain",angle = +45, hjust = +0.5)  # Rotate individual violin plot titles counterclockwise by 20 degrees and adjust horizontal justification
#   )

DefaultAssay(merged_seurat) = "RNA"
DotPlot(merged_seurat,assay = "RNA" ,features = c("Esr1","Phf6",selected_genes), cols = c("white", "red"), dot.scale = 12)+
  RotatedAxis()

dotPlot <- DotPlot(merged_seurat, features = "Phf6", cols = c("white", "red"), dot.scale = 12)
dotPlot

plot_data <- data.frame(dotPlot$data,stringsAsFactors = FALSE) %>%
  select(pct.exp, id)

# 将百分比转换为小数
pct <- plot_data$pct.exp

plot_data$id <- as.numeric(plot_data$id)

#piechart

# 设置图形布局，2行5列
par(mfrow=c(2, 5))

# 指定颜色
custom_colors <- c()

# 为每个百分比绘制饼图
for (i in 1:length(pct)) {
  pie(c(pct[i], 100 - pct[i]), labels = "", col = c(custom_colors[i], "white"),lwd = 1,border = "#bcbdc0")
  title(paste("Pie Chart ", i))
}

# 恢复默认图形布局
par(mfrow=c(1, 1))

# seurat_clusters <- as.character(plot_data$id)  # Convert to character if not already

# Create a pie chart using ggplot2

p<-FeaturePlot(
  merged_seurat,
  features = "Phf6",
  label = F
)+
  plot_annotation(theme = theme(plot.title = element_text(size = 12,hjust = 0.5, face = "bold")))+
  scale_color_gradientn(colours  = c("#BDBDBD", "lightgrey", "#d1c4e9", "#7e57c2", "#5e35b1","#4527a0","blue")
)

# a<-p+
#   labs(title = "")
# a
DimPlot(merged_seurat, label = F) 

# 使用 plot_density 函数绘制密度图
plot_density(merged_seurat, features = "Phf6")

# Assuming that 'Phf6' is a column in your Seurat object representing the expression of the Phf6 gene

# Create a binary indicator for Phf6 expression greater than 0
merged_seurat$Phf6_positive <- as.integer(merged_seurat$Phf6 > 0)

# 
