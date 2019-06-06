---
title: Analysis of single cell RNA-seq data
tag: single cell
categories: single_cell
---


### 1 Tabula Muris
#### 1.1 介绍
他们将高通量但低覆盖率的10倍数据与低通量高覆盖率的数据结合起来，FACS-sortedcell+Smartseq2。
#### 1.22 数据下载  
下载FACS data
<!--more-->
```
wget https://ndownloader.figshare.com/files/10038307
unzip 10038307
wget https://ndownloader.figshare.com/files/10038310
mv 10038310 FACS_metadata.csv
wget https://ndownloader.figshare.com/files/10039267
mv 10039267 FACS_annotations.csv
```
下载 10X data

```
wget https://ndownloader.figshare.com/files/10038325
unzip 10038325
wget https://ndownloader.figshare.com/files/10038328
mv 10038328 droplet_metadata.csv
wget https://ndownloader.figshare.com/files/10039264
mv 10039264 droplet_annotation.csv
```
下载好的数据如下所示：
![DownLoadData](https://note.youdao.com/yws/public/resource/12beaf1b878bc1059015bf8c805c27d5/xmlnote/A5BCB784861A4B3CB16854DDB816A6AC/889)  
现在有两个文件夹：“FACS”和“droplet”，每一个文件夹包括一个annotation和metadata文件。要检查这些文件，您可以使用head查看文本文件的前几行（按“q”键退出）：

```
head -n 10 droplet_metadata.csv
```

```
## channel,mouse.id,tissue,subtissue,mouse.sex
## 10X_P4_0,3-M-8,Tongue,NA,M
## 10X_P4_1,3-M-9,Tongue,NA,M
## 10X_P4_2,3-M-8/9,Liver,hepatocytes,M
## 10X_P4_3,3-M-8,Bladder,NA,M
## 10X_P4_4,3-M-9,Bladder,NA,M
## 10X_P4_5,3-M-8,Kidney,NA,M
## 10X_P4_6,3-M-9,Kidney,NA,M
## 10X_P4_7,3-M-8,Spleen,NA,M
## 10X_P7_0,3-F-56,Liver,NA,F
```
#### 1.3 读取这些数据集（Smartseq2）
现在，我们可以从逗号分隔的文件中读取相关的计数矩阵。然后检查产生的dataframe：

```
dat = read.delim("FACS/Kidney-counts.csv", sep=",", header=TRUE)
dat[1:5,1:5]
```

```
##               X A14.MAA000545.3_8_M.1.1 E1.MAA000545.3_8_M.1.1
## 1 0610005C13Rik                       0                      0
## 2 0610007C21Rik                       1                      0
## 3 0610007L01Rik                       0                      0
## 4 0610007N19Rik                       0                      0
## 5 0610007P08Rik                       0                      0
##   M4.MAA000545.3_8_M.1.1 O21.MAA000545.3_8_M.1.1
## 1                      0                       0
## 2                      0                       0
## 3                      0                       0
## 4                      0                       0
## 5                      0                       0
```
我们可以看到dataframe的第一列是基因名，所以首先我们把它们移到行名中这样我们就有了一个数字矩阵：

```
dim(dat)
```

```
## [1] 23433   866
```

```
rownames(dat) <- dat[,1]
dat <- dat[,-1]
```
由于这是Smartseq2的数据集，所以他可能包含spike-ins：

```
rownames(dat)[grep("^ERCC-", rownames(dat))]
```

```
##  [1] "ERCC-00002" "ERCC-00003" "ERCC-00004" "ERCC-00009" "ERCC-00012"
##  [6] "ERCC-00013" "ERCC-00014" "ERCC-00016" "ERCC-00017" "ERCC-00019"
## [11] "ERCC-00022" "ERCC-00024" "ERCC-00025" "ERCC-00028" "ERCC-00031"
## [16] "ERCC-00033" "ERCC-00034" "ERCC-00035" "ERCC-00039" "ERCC-00040"
## [21] "ERCC-00041" "ERCC-00042" "ERCC-00043" "ERCC-00044" "ERCC-00046"
## [26] "ERCC-00048" "ERCC-00051" "ERCC-00053" "ERCC-00054" "ERCC-00057"
## [31] "ERCC-00058" "ERCC-00059" "ERCC-00060" "ERCC-00061" "ERCC-00062"
## [36] "ERCC-00067" "ERCC-00069" "ERCC-00071" "ERCC-00073" "ERCC-00074"
## [41] "ERCC-00075" "ERCC-00076" "ERCC-00077" "ERCC-00078" "ERCC-00079"
## [46] "ERCC-00081" "ERCC-00083" "ERCC-00084" "ERCC-00085" "ERCC-00086"
## [51] "ERCC-00092" "ERCC-00095" "ERCC-00096" "ERCC-00097" "ERCC-00098"
## [56] "ERCC-00099" "ERCC-00104" "ERCC-00108" "ERCC-00109" "ERCC-00111"
## [61] "ERCC-00112" "ERCC-00113" "ERCC-00116" "ERCC-00117" "ERCC-00120"
## [66] "ERCC-00123" "ERCC-00126" "ERCC-00130" "ERCC-00131" "ERCC-00134"
## [71] "ERCC-00136" "ERCC-00137" "ERCC-00138" "ERCC-00142" "ERCC-00143"
## [76] "ERCC-00144" "ERCC-00145" "ERCC-00147" "ERCC-00148" "ERCC-00150"
## [81] "ERCC-00154" "ERCC-00156" "ERCC-00157" "ERCC-00158" "ERCC-00160"
## [86] "ERCC-00162" "ERCC-00163" "ERCC-00164" "ERCC-00165" "ERCC-00168"
## [91] "ERCC-00170" "ERCC-00171"
```
现在，我们可以从列名中提取该数据的大部分元数据：

```
cellIDs <- colnames(dat)
cell_info <- strsplit(cellIDs, "\\.")
Well <- lapply(cell_info, function(x){x[1]})
Well <- unlist(Well)
Plate <- unlist(lapply(cell_info, function(x){x[2]}))
Mouse <- unlist(lapply(cell_info, function(x){x[3]}))
```
我们可以检查每个元数据分类的分布：

```
summary(factor(Mouse))
```

```
## 3_10_M 3_11_M 3_38_F 3_39_F  3_8_M  3_9_M 
##    104    196    237    169     82     77
```
我们还可以检查是否有任何技术因素被混淆：

```
table(Mouse, Plate)
```

```
##         Plate
## Mouse    B001717 B002775 MAA000545 MAA000752 MAA000801 MAA000922
##   3_10_M       0       0         0       104         0         0
##   3_11_M       0       0         0         0       196         0
##   3_38_F     237       0         0         0         0         0
##   3_39_F       0     169         0         0         0         0
##   3_8_M        0       0        82         0         0         0
##   3_9_M        0       0         0         0         0        77
```
最后，我们将推断的单元类型注释，并将它们与表达式矩阵中的单元相匹配：

```
ann <- read.table("FACS_annotations.csv", sep=",", header=TRUE)
ann <- ann[match(cellIDs, ann[,1]),]
celltype <- ann[,3]
```
#### 1.4 构建一个scater object
为了创建一个SingleCellExperiment对象，我们必须把所有的单元注释放在一个dataframe中，因为实验性的批处理（pcr板）与供体鼠标完全混淆了，我们只保留其中一个。

```
require("SingleCellExperiment")
```

```
cell_anns <- data.frame(mouse = Mouse, well=Well, type=celltype)
rownames(cell_anns) <- colnames(dat)
sceset <- SingleCellExperiment(assays = list(counts = as.matrix(dat)), colData=cell_anns)
```
最后，如果数据集包含了“spike”——我们在singlecellexperobject对象中有一个隐藏的变量来跟踪它们：

```
isSpike(sceset, "ERCC") <- grepl("ERCC-", rownames(sceset))
```
#### 1.5 读取数据集（10X）
由于其大小和稀疏的10X数据（高达90%的表达式矩阵可能是0），它通常存储为稀疏矩阵。CellRanger的默认输出格式是一个.mtx文件，它将这个稀疏矩阵存储为一列行坐标、一列列坐标，以及一列大于0的表达式值。注意，如果您查看.mtx文件，您将看到两个标题行，后面是一行，详细描述了完整矩阵的行数、列数和计数。因为只有坐标被存储在.mtx文件中，所以每一行和列的名称必须分别存储在“genes.tsv”和“barcode.tsv”。

```
require("Matrix")
```

```
cellbarcodes <- read.table("droplet/Kidney-10X_P4_5/barcodes.tsv")
genenames <- read.table("droplet/Kidney-10X_P4_5/genes.tsv")
molecules <- Matrix::readMM("droplet/Kidney-10X_P4_5/matrix.mtx")
```
现在我们将添加适当的行和列名。但是，如果您检查读取的cellbarcode，您将看到它们只是与每个单元相关联的条形码序列。这是一个问题，因为每批10个数据都使用相同的条形码池，因此如果我们需要将来自多个10个批次的数据组合起来，那么这些数据将不会是唯一的。因此，我们将把批号附加到每个细胞的条形码上：

```
head(cellbarcodes)
```

```
##                   V1
## 1 AAACCTGAGATGCCAG-1
## 2 AAACCTGAGTGTCCAT-1
## 3 AAACCTGCAAGGCTCC-1
## 4 AAACCTGTCCTTGCCA-1
## 5 AAACGGGAGCTGAACG-1
## 6 AAACGGGCAGGACCCT-1
```

```
rownames(molecules) <- genenames[,1]
colnames(molecules) <- paste("10X_P4_5", cellbarcodes[,1], sep="_")
```
现在，让我们获取这些数据的元数据和注释：

```
meta <- read.delim("droplet_metadata.csv", sep=",", header=TRUE)
head(meta)
```

```
##    channel mouse.id  tissue   subtissue mouse.sex
## 1 10X_P4_0    3-M-8  Tongue        <NA>         M
## 2 10X_P4_1    3-M-9  Tongue        <NA>         M
## 3 10X_P4_2  3-M-8/9   Liver hepatocytes         M
## 4 10X_P4_3    3-M-8 Bladder        <NA>         M
## 5 10X_P4_4    3-M-9 Bladder        <NA>         M
## 6 10X_P4_5    3-M-8  Kidney        <NA>         M
```
在这里我们可以看到,我们需要使用“10 x_p4_5”找到这批的元数据,还要注意,老鼠ID是不同的格式在此元数据表为连字符而不是下划线和性别在中间的ID。检查附带论文的方法部分我们知道相同的8个老鼠用于droplet和plate-based技术。因此，我们需要修复小鼠的id，使其与FACS实验中使用的id一致。
```
meta[meta$channel == "10X_P4_5",]
```

```
##    channel mouse.id tissue subtissue mouse.sex
## 6 10X_P4_5    3-M-8 Kidney      <NA>         M
```

```
mouseID <- "3_8_M"
```
注意：根据你所分配的组织，你可能有来自混合样本的10倍数据：例如，鼠标id=3-m-5/6。您仍然应该重新格式化这些内容以保持一致，但是它们不会匹配来自FACS数据的鼠标id，这可能会影响您的下游分析。如果这些小鼠不是来自于近亲繁殖的病毒株，那么就有可能将单个细胞分配给特定的老鼠，使用的是无证明的snps，但这超出了本课程的范围。

```
ann <- read.delim("droplet_annotation.csv", sep=",", header=TRUE)
head(ann)
```

```
##                        cell  tissue cell_ontology_class
## 1 10X_P4_3_AAAGTAGAGATGCCAG Bladder    mesenchymal cell
## 2 10X_P4_3_AACCGCGTCCAACCAA Bladder    mesenchymal cell
## 3 10X_P4_3_AACTCCCGTCGGGTCT Bladder    mesenchymal cell
## 4 10X_P4_3_AACTCTTAGTTGCAGG Bladder        bladder cell
## 5 10X_P4_3_AACTCTTTCATAACCG Bladder    mesenchymal cell
## 6 10X_P4_3_AAGACCTAGATCCGAG Bladder    endothelial cell
##                      cell_ontology_term_iri cell_ontology_id
## 1 http://purl.obolibrary.org/obo/CL_0008019       CL:0008019
## 2 http://purl.obolibrary.org/obo/CL_0008019       CL:0008019
## 3 http://purl.obolibrary.org/obo/CL_0008019       CL:0008019
## 4 http://purl.obolibrary.org/obo/CL_1001319       CL:1001319
## 5 http://purl.obolibrary.org/obo/CL_0008019       CL:0008019
## 6 http://purl.obolibrary.org/obo/CL_0000115       CL:0000115
```
再一次，你会发现在注释和cell条形码之间有一个细微的区别，在匹配之前，我们必须纠正它们。

```
ann[,1] <- paste(ann[,1], "-1", sep="")
ann_subset <- ann[match(colnames(molecules), ann[,1]),]
celltype <- ann_subset[,3]
```
现在让我们构建单元元数据dataframe：

```
cell_anns <- data.frame(mouse = rep(mouseID, times=ncol(molecules)), type=celltype)
rownames(cell_anns) <- colnames(molecules);
```
#### 1.6 构建一个scater object
既然我们已经读入10X的数据，我们需要将它们组合成一个SingCellExperiment对象。首先，我们将检查基因名称是否相同，并且在所有批次中都是相同的：

```
identical(rownames(molecules1), rownames(molecules2))
## [1] TRUE
identical(rownames(molecules1), rownames(molecules3))
## [1] TRUE
```
现在我们来检查一下是否有重复的CellIDs：

```
sum(colnames(molecules1) %in% colnames(molecules2))
## [1] 0
sum(colnames(molecules1) %in% colnames(molecules3))
## [1] 0
sum(colnames(molecules2) %in% colnames(molecules3))
## [1] 0
```
一切都OK，我们可以把它们组合起来：
```
all_molecules <- cbind(molecules1, molecules2, molecules3)
all_cell_anns <- as.data.frame(rbind(cell_anns1, cell_anns2, cell_anns3))
all_cell_anns$batch <- rep(c("10X_P4_5", "10X_P4_6","10X_P7_5"), times = c(nrow(cell_anns1), nrow(cell_anns2), nrow(cell_anns3)))
```
现在我们来构建一个SingleCellExperiments对象。SingleCellExperiment类的一个优点是,它能够将数据存储为正常数据矩阵或稀疏矩阵格式,以及HDF5格式其允许大型non-sparse矩阵存储以及以一种有效的方式访问磁盘,而不是整个加载到RAM中。

```
require("SingleCellExperiment")
require("scater")
all_molecules <- as.matrix(all_molecules)
sceset <- SingleCellExperiment(assays = list(counts = as.matrix(all_molecules)), colData=all_cell_anns)
```
因为这是10倍的数据，所以它不会包含spike-in，所以我们只保存数据：

```
saveRDS(sceset, "kidney_droplet.rds")
```
### 2 cleaning the Expression Matrix
#### 2.1 ExpressionQC（UMI）
##### 2.1.1 介绍
一旦基因表达被量化，它就被总结为一个表达矩阵，其中每一行对应一个基因（或转录），每一列对应一个单元格。应该检查这个矩阵，以去除那些没有被检测到的质量好的细胞，这些细胞在读取QC或映射质量控制步骤中没有被检测到。在这个阶段，如果不能去除低质量的细胞，可能会增加技术噪音，这种噪音有可能掩盖下游分析中感兴趣的生物信号。  
由于目前还没有标准的方法来执行scRNASeq，因此在这里提出的各种质量控制措施的预期值可能会因实验而异。因此，为了执行QC，我们将寻找与数据集其余部分相关的异常值的单元，而不是与独立的质量标准进行比较。因此，在比较使用不同协议收集的数据集的质量度量时，应该小心谨慎。
##### 2.1.2 Tung dataset
为了说明细胞的质量控制，[我们考虑了一个由三个不同的个体产生的诱导多功能干细胞的数据集](http://jdblischak.github.io/singleCellSeq/analysis/)（Tung等人，2017年），在芝加哥大学的Yoav Gilad实验室。这些实验是在FluidigmC1平台上进行的，并促进了对唯一的分子标识符（UMIs）和ERCC spike的量化

```
library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)
```
加载数据和注释文件：

```
molecules <- read.table("tung/molecules.txt", sep = "\t")
anno <- read.table("tung/annotation.txt", sep = "\t", header = TRUE)
```
检查表达式矩阵的一小部分

```
head(molecules[ , 1:3])
```

```
##                 NA19098.r1.A01 NA19098.r1.A02 NA19098.r1.A03
## ENSG00000237683              0              0              0
## ENSG00000187634              0              0              0
## ENSG00000188976              3              6              1
## ENSG00000187961              0              0              0
## ENSG00000187583              0              0              0
## ENSG00000187642              0              0              0
```

```
head(anno)
```

```
##   individual replicate well      batch      sample_id
## 1    NA19098        r1  A01 NA19098.r1 NA19098.r1.A01
## 2    NA19098        r1  A02 NA19098.r1 NA19098.r1.A02
## 3    NA19098        r1  A03 NA19098.r1 NA19098.r1.A03
## 4    NA19098        r1  A04 NA19098.r1 NA19098.r1.A04
## 5    NA19098        r1  A05 NA19098.r1 NA19098.r1.A05
## 6    NA19098        r1  A06 NA19098.r1 NA19098.r1.A06
```
我们通过使用SingleCellExperiments（SCE）和scater包来对分析进行标准化。首先，创建SCE对象：

```
umi <- SingleCellExperiment(
    assays = list(counts = as.matrix(molecules)), 
    colData = anno
)
```
移除在任何细胞中没有表达的基因：

```
keep_feature <- rowSums(counts(umi) > 0) > 0
umi <- umi[keep_feature, ]
```
定义控制features（genes）-ERCC spike和线粒体基因：

```
isSpike(umi, "ERCC") <- grepl("^ERCC-", rownames(umi))
isSpike(umi, "MT") <- rownames(umi) %in% 
    c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
    "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
    "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
    "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
    "ENSG00000198840")
```
计算质量指标:

```
umi <- calculateQCMetrics(
    umi,
    feature_controls = list(
        ERCC = isSpike(umi, "ERCC"), 
        MT = isSpike(umi, "MT")
    )
)
```
##### 2.1.3 Cell QC
###### 2.1.3.1 文库大小（Library size）
接下来，我们考虑每个样本检测到的RNA分子总数（如果我们使用的是读计数而不是UMI计数，这将是read的总数）。没有多少reads/molecules的cells可能会被破坏或无法捕获一个细胞，因此应该被移除。

```
hist(
    umi$total_counts,
    breaks = 100
)
abline(v = 25000, col = "red")

```
![histplot](https://note.youdao.com/favicon.ico)
###### 2.1.3.2 Detected genes
除了为每个样本确保足够的排序深度之外，我们还希望确保read分布在transcriptome中。因此，我们计算每个样本中检测到的唯一基因的总数。

```
hist(
    umi$total_features,
    breaks = 100
)
abline(v = 7000, col = "red")
```
从图中我们得出结论，大多数细胞有7000-10000个检测到的基因，这对于high-depth的scRNA-seq来说是正常的。然而，这取决于实验流程和测序深度。例如，基于droplet的方法或具有较低测序深度的样本通常检测每个细胞的基因更少。上面的图中最显著的特征是在分布的左边的“重尾（heavy tail）”。如果在整个细胞中检测率是相等的，那么分布应该是近似正常的。因此，我们将这些细胞移出分布的尾部（少于7000个检测到的基因）。

```
filter_by_total_counts <- (umi$total_counts > 25000)
```


```
table(filter_by_total_counts)
```


```
## filter_by_total_counts
## FALSE  TRUE 
##    46   818
```

```
filter_by_total_features <- (umi$total_features > 7000)
```

```
## filter_by_expr_features
## FALSE  TRUE 
##   116   748
```

###### 2.1.3.3 ERCCs和MTS
另一种衡量细胞质量的指标是ERCC spike-in rna和内源性rna之间的比率。这个比率可以用来估计被捕获的细胞中RNA的总量。具有高水平 spike-in RNA的细胞具有较低的RNA起始数量，这可能是由于细胞死亡或压力导致RNA被降解。

```
plotPhenoData(
    umi,
    aes_string(
        x = "total_features",
        y = "pct_counts_MT",
        colour = "batch"
    )
)
```

![image](https://hemberg-lab.github.io/scRNA.seq.course/14-exprs-qc_files/figure-html/mt-vs-counts-1.png)

```
plotPhenoData(
    umi,
    aes_string(
        x = "total_features",
        y = "pct_counts_ERCC",
        colour = "batch"
    )
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/14-exprs-qc_files/figure-html/ercc-vs-counts-1.png)


```
filter_by_MT <- (umi$total_counts_MT/umi$total_counts<0.1)
```

```
table(filter_by_MT)
##filter_by_MT
##FALSE  TRUE 
##   31   833 
```

```
filter_by_ERCC<- (umi$batch!="NA19098.r2")
```

```
## filter_by_ERCC
## FALSE  TRUE 
##    96   768
```
##### 2.1.4 Cell filtering
###### 2.1.4.1 Manual
现在我们可以根据我们之前的分析定义一个细胞过滤：

```
umi$use <- (
    # sufficient features (genes)
    filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # sufficient endogenous RNA
    filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    filter_by_MT
)
```

```
table(umi$use)
## 
## FALSE  TRUE 
##   207   657
```
###### 2.1.4.2 Automatic
scater的另一种选择是在一组QC指标上进行PCA，然后使用自动离群检测来识别潜在的问题细胞。
默认情况下，以下指标用于基于PCA的异常检测:
- pct_counts_top_100_features
- total_features
- pct_counts_feature_controls
- n_detected_feature_controls
- log10_counts_endogenous_features
- log10_counts_feature_controls
scater首先创建一个矩阵其中行表示cell，列表示不同的QC指标.在这里，PCA图提供了由其质量指标排序的细胞的2D表示。然后使用mvoutlier包中的方法检测出异常值

```
umi <- plotPCA(
    umi,
    size_by = "total_features", 
    shape_by = "use",
    pca_data_input = "pdata",
    detect_outliers = TRUE,
    return_SCE = TRUE
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/14-exprs-qc_files/figure-html/auto-cell-filt-1.png)
##### 2.1.5 Compare filterings
从limma包中使用venn计数和vennDiagram函数来制作维恩图。
##### 2.1.6 Gene analysis
###### 2.1.6.1 Gene expression
除了去除质量不佳的细胞外，在我们怀疑技术人工制品可能扭曲了结果的情况下，排除基因是一个好主意。此外，对基因表达谱的检查可以提供关于如何改进实验过程的见解。  
考虑前50个表达基因所消耗的reads通常是有益的。

```
plotQC(umi, type = "highest-expression")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/14-exprs-qc_files/figure-html/top50-gene-expr-1.png)
这些分布是相对平坦的表示（但不是保证！）对这些cell的完整的转录组的良好覆盖。然而，在前15个基因中，有几只“spike-ins”的基因表明，如果实验要重复的话，对spike-ins的稀释度可能会更大。
###### 2.1.6.2 Gene filtering
去除那些表达水平被认为是“不可检测”的基因是一个好主意。如果1个基因的转录本在至少两个细胞中被检测到，我们就认为该基因是可以被检测的。如果我们考虑的是reads，而不是UMI，那么一个合理的阈值就是至少需要在两个cell中至少有5个reads。然而，在这两种情况下，阈值都依赖于测序深度。重要的是要记住，基因必须在细胞过滤后进行过滤，因为有些基因只能在质量较差的细胞中被检测出来（注：colData（umi）$use过滤器应用于umi数据集）

```
filter_genes <- apply(
    counts(umi[ , colData(umi)$use]), 
    1, 
    function(x) length(x[x > 1]) >= 2
)
rowData(umi)$use <- filter_genes
```

```
table(filter_genes)
```

```
## filter_genes
## FALSE  TRUE 
##  4660 14066
```
##### 2.1.7 Save the Data
QCed dataset数据集的维度（不要忘记我们上面定义的基因过滤参数）：

```
dim(umi[rowData(umi)$use, colData(umi)$use])
```

```
## [1] 14066   657
```
让我们用log-treansformed count创建一个额外的slot（我们将在下一章中需要它），并从还原的slot中删除保存的PCA结果。

```
assay(umi, "logcounts_raw") <- log2(counts(umi) + 1)
reducedDim(umi) <- NULL
```
保存数据

```
saveRDS(umi, file = "tung/umi.rds")
```
#### 2.2 Expression QC（Reads）

```
library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)
```

同上面的一节
#### 2.3 Data visualization
##### 2.3.1 介绍
single-cell RNA-seq 的一个重要方面是控制批处理影响。批处理效果是在处理过程中添加到样品中的由人为产生的影响。例如，如果在不同的实验室或在同一实验室的不同日子里准备了两组样本，那么我们可能会观察到在一起处理的样本之间有更大的相似性。在最坏的情况下，批处理可能被误认为是真正的生物变异。Tung数据使我们能够以可控的方式探索这些问题，因为这些样本的处理方式的一些重要方面已经被记录下来。理想情况下，我们希望看到来自同一分组的批次，以及对应于每个个体的不同组。

```
library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)
umi <- readRDS("tung/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control
```
##### 2.3.2 PCA Plot
对数据进行概述的最简单方法是使用主成分分析转换它，然后将前两个主要组件可视化。
主成分分析（PCA）是一个统计过程，它使用一个转换将一组观察转换成一组称为主成分（pc）的线性无关变量的值。主成分的数量小于或等于原始变量的数量。
从数学上说，pc对应于协方差矩阵的特征向量。特征向量是按特征值排序的所以第一个主成分在数据中尽可能多地反映了数据的可变性，而每一个后续的组件在约束条件下都有最高的方差它与前面的组件是正交的
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/pca.png)
###### 2.3.2.1 Before QC
在做log转换之前：

```
plotPCA（
umi[endog_genens,],
exprs_values="counts",
colour_by="batch",
size_by="total_features",
shape_by="individiual"
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/16-exprs-overview_files/figure-html/expr-overview-pca-before-qc1-1.png)
做log转换之后：

```
plotPCA(
    umi[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/16-exprs-overview_files/figure-html/expr-overview-pca-before-qc2-1.png)
很明显，对数转换对我们的数据是有利的——它减少了第一个主成分的方差，并且已经分离了一些生物效应。此外，它使表达式的分布更加正常。在接下来的分析和章节中，我们将在默认情况下使用日志转换的原始计数。
但是，请注意，仅仅一个对数转换不足以解释单元之间的不同技术因素（例如，排序深度）。因此，请不要使用logcounts_raw来进行下游分析，而是将其作为最小的合适的数据使用，它不仅是对数转换的，而且还被库的大小（如CPM normalisation）规范化。在本课程中，我们仅为演示目的使用logcounts-raw格式！
###### 2.3.2.2 After QC
plotPCA(
    umi.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
![image](https://hemberg-lab.github.io/scRNA.seq.course/16-exprs-overview_files/figure-html/expr-overview-pca-after-qc-1.png)
比较上面两图，很明显，在质量控制之后，NA19098.r2中的细胞不再形成一组异常值。  
默认情况下，scater使用最多的500个最可变的基因来计算PCA。这可以通过改变ntop参数进行调整

```
plotPCA(umi.qc[endog_genes,],exprs_values="logcounts_raw",colour_by="batch",size_by="total_features",shape_by="indivaidual",ntop=50)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/16-exprs-overview_files/figure-html/expr-overview-pca-after-qc-exercise1-2-1.png)
##### 2.3.3 tSNE map

替代PCA的将scRNASeq数据可视化另一种方法是tSNE图。tSNE（t-分布式随机邻居嵌入）将维度减少（例如PCA）与在邻近的网络上随机漫步，将高维数据（即我们的14,214维表达式矩阵）映射到二维空间，同时保持细胞之间的距离。与PCA相比，tSNE是一种随机算法，这意味着在同一数据集上多次运行该方法将导致不同的绘图块与PCA相比，tSNE是一种随机算法，这意味着在同一数据集上多次运行该方法将导致不同的绘图块。由于该算法的非线性和随机性，tSNE更难直观地解释tSNE。为了确保重现性，我们在下面的代码中修复了随机数字生成器的“seed”，这样我们就可以得到相同的图。
###### 2.3.3.1 Before QC

```
plotTSNE(
    umi[endog_genes, ],
    exprs_values = "logcounts_raw",
    perplexity = 130,
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    rand_seed = 123456
)
```

![image](https://hemberg-lab.github.io/scRNA.seq.course/16-exprs-overview_files/figure-html/expr-overview-tsne-before-qc-1.png)
###### 2.3.3.2 After QC

```
plotTSNE(
    umi.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    perplexity = 130,
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    rand_seed = 123456
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/16-exprs-overview_files/figure-html/expr-overview-tsne-after-qc-1.png)

解释PCA和tSNE Plot通常具有挑战性，由于它们的随机和非线性特性，他们不太直观。然而，在这种情况下，很明显，它们提供了类似的数据图像。比较上面两图，同样清楚的是来自NA19098.r2的样本在QC过滤之后不再是异常值。
此外，tSNE要求你提供perplexity值，它反映了用来建立近邻网络的邻居的数量;高值创造了一个密集的网络，将细胞聚集在一起，而低的值使网络更加稀疏，使得细胞群相互分离。
#### 2.4 Data visualization（Reads）

```
library(scater)
options(stringsAsFactors = FALSE)
reads <- readRDS("tung/reads.rds")
reads.qc <- reads[rowData(reads)$use, colData(reads)$use]
endog_genes <- !rowData(reads.qc)$is_feature_control
```

```
plotPCA(
    reads[endog_genes, ],
    exprs_values = "counts",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/17-exprs-overview-reads_files/figure-html/expr-overview-pca-before-qc-reads1-1.png)

```
plotPCA(
    reads[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/17-exprs-overview-reads_files/figure-html/expr-overview-pca-before-qc-reads2-1.png)

```
plotPCA(
    reads.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/17-exprs-overview-reads_files/figure-html/expr-overview-pca-after-qc-reads-1.png)

```
plotTSNE(
    reads[endog_genes, ],
    exprs_values = "logcounts_raw",
    perplexity = 130,
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    rand_seed = 123456
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/17-exprs-overview-reads_files/figure-html/expr-overview-tsne-before-qc-reads-1.png)

```
plotTSNE(
    reads.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    perplexity = 130,
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    rand_seed = 123456
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/17-exprs-overview-reads_files/figure-html/expr-overview-tsne-after-qc-reads-1.png)
tSNE map of the tung data
![image](https://hemberg-lab.github.io/scRNA.seq.course/17-exprs-overview-reads_files/figure-html/expr-overview-tsne-after-qc-exercise2-1-1.png)
tSNE map of the tung data (perplexity = 10)
![image](https://hemberg-lab.github.io/scRNA.seq.course/17-exprs-overview-reads_files/figure-html/expr-overview-tsne-after-qc-exercise2-2-1.png)
tSNE map of the tung data (perplexity = 200)
#### 2.5 Identifying confounding factors
##### 2.5.1 Introduction
在前几章中，我们考虑了批处理效果，在这一章中，我们将继续探索如何识别和移除实验中的影响。我们将继续使用scater包，因为它提供了一组专门用于测试和解释变量的质量控制的方法。此外，我们将继续使用前一章中使用的Blischak数据。

```
library(scater, quietly = TRUE)
options(stringsAsFactors = FALSE)
umi <- readRDS("tung/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control
```
umi.qc数据集包含经过细胞和基因过滤的数据。
##### 2.5.2 Correlations with PCs
我们再回顾一下PCA 图以及QCed数据集

```
plotPCA(
    umi.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "batch",
    size_by = "total_features"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/18-confounders_files/figure-html/confound-pca-1.png)
scater允许一个识别与实验和QC变量相关的主要成分
让我们来测试一些变量是否与任何pc相关联。
###### 2.5.2.1 Detected genes

```
plotQC(
    umi.qc[endog_genes, ],
    type = "find-pcs",
    exprs_values = "logcounts_raw",
    variable = "total_features"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/18-confounders_files/figure-html/confound-find-pcs-total-features-1.png)
事实上，我们可以看到PC1几乎可以完全由检测到的基因的数量来解释。事实上，在上面的PCA图上也可以看出。
scater也可以对每一个变量计算mariginal R2，当拟合线性模型时，每个基因的表达值都与这个变量相对应，并显示出基因边缘的密度图。

```
plotQC(
    umi.qc[endog_genes, ],
    type = "expl",
    exprs_values = "logcounts_raw",
    variables = c(
        "total_features",
        "total_counts",
        "batch",
        "individual",
        "pct_counts_ERCC",
        "pct_counts_MT"
    )
)

```
![image](https://hemberg-lab.github.io/scRNA.seq.course/18-confounders_files/figure-html/confound-find-expl-vars-1.png)
这一分析表明，检测到的基因（再次）和测序深度（计数的数量）的数量对许多基因具有重要的解释力，因此，这些变量是在标准化过程或下游的统计模型中进行调整的候选目标哦。ERCCs的表达似乎也是一个重要的解释变量，上面的一个显著特征是，批处理比individual更能解释。
#### 2.6 Identifying confounding factors (Reads)
处理同上
#### 2.7 Normalization theory
##### 2.7.1 介绍
在前一章中，我们确定了重要的混淆因素和解释变量。scater允许在随后的统计模型中解释这些变量，或者在需要的情况下使用normaliseExprs（）对其进行条件去除。这可以通过为normaliseExprs（）提供一个设计矩阵来实现。
##### 2.7.2 文库大小
由于scRNA-seq数据通常是在高度多路复用的平台上进行测序的，所以库的大小各不相同，从每个单元中获得的总读数可能会有很大的不同。一些(如量化方法当确定基因表达估计时（如Cufflinks、RSEM）将库的大小合并起来，因此处理这些软件输出的数据时不需要这种标准化。然而，如果使用另一种量化方法，那么必须通过将表达式矩阵的每一列乘以或除以一个“标准化因子”来修正库的大小，这是相对于其他cell的库大小的估计。许多修正库大小的方法已经被开发出来用于bulkRNA-seq，并且可以同样应用于scRNA-seq（例如。UQ，SF，CPM，RPKM，FPKM，TPM）。
##### 2.7.3 Normalisations
###### 2.7.3.1 CPM
使这些数据规范化的最简单方法是将每一列除以它的总数，然后乘以100 000，将它转换为每百万（CPM）计数。请注意，为了纠正总细胞RNA的含量，应该将spike排除在总表达的计算中，因此我们只使用内源性基因。

```
calc_cpm <-
function (expr_mat, spikes = NULL) 
{
    norm_factor <- colSums(expr_mat[-spikes, ])
    return(t(t(expr_mat)/norm_factor)) * 10^6
}
```
CPM的一个潜在缺点是，如果你的样本包含了在细胞中高度表达和不同表达的基因。在这种情况下，细胞内的总分子可能依赖于这些基因在细胞中是否开启，而由总分子正常化可能会隐藏这些基因的差异表达，或者为剩余的基因错误地创造出不同的表达。注意，RPKM、FPKM和TPM是CPM的变体，它根据各自的基因/转录本的长度进一步调整计数。
###### 2.7.3.2 RLE（SF）
size factr（SF）是由DESeq提出并推广的。首先计算出所有细胞中每个基因的几何平均值。每个细胞的大小因子是表达与基因的几何平均比率的基因的中值。这种方法的一个缺点是，由于它使用的是几何平均值，所以在所有的细胞中，只有非零表达值的基因可以被用于计算，这使得对于大型低深度的scRNASeq实验来说是不可取的。edgeR&scater称该方法为“relative log expression”。

```
calc_sf <-
function (expr_mat, spikes = NULL) 
{
    geomeans <- exp(rowMeans(log(expr_mat[-spikes, ])))
    SF <- function(cnts) {
        median((cnts/geomeans)[(is.finite(geomeans) & geomeans > 
            0)])
    }
    norm_factor <- apply(expr_mat[-spikes, ], 2, SF)
    return(t(t(expr_mat)/norm_factor))
}
```
###### 2.7.3.3 UQ
上四分位数（upperquartile，UQ）是由（Bullard et al.2010）提出的。这里每一列都除以每个库的75%的计数。通常，计算出的分位数是由单元间的中值缩放的，以保持表达式的绝对水平相对一致。这种方法的一个缺点是，对于低深度的scRNASeq实验，大量未被发现的基因可能导致75%的量子数为零（或接近它）。这个限制可以通过概括这个想法和使用一个更高的分位数来克服（例如。99%的分位数是scater的默认值），或者在计算75%的分位数之前排除0。

```
calc_uq <-
function (expr_mat, spikes = NULL) 
{
    UQ <- function(x) {
        quantile(x[x > 0], 0.75)
    }
    uq <- unlist(apply(expr_mat[-spikes, ], 2, UQ))
    norm_factor <- uq/median(uq)
    return(t(t(expr_mat)/norm_factor))
}
```
###### 2.7.3.4 TMM
另一种叫做TMM的方法是由（罗宾逊和Oshlack 2010）提出的m值（对参考）的加权平均化平均值。问题的m值是基因的log2在单个细胞之间的变化。一个cell用作参考，然后计算每个单元格的m值与这个引用相比较。然后，通过删除顶部和底部30%的值来调整这些值，并且通过加权它们来计算剩余值的平均值，以考虑对数尺度对方差的影响。每个非引用cell乘以计算因子。这种方法的两个潜在问题是在修剪后留下的非零基因，以及大多数基因没有差异表达的假设。
###### 2.7.3.5 scran
scran包实现了针对single-cell数据（l.Lun、巴赫和Marioni 2016）的CPM的个性化。简单地说，这个方法解决了单个cell的大量零值的问题，通过将cell计算一个标准化因子（类似于CPM）来计算每个pool的总和。由于一个cell存在于许多不同的pool中，所以可以使用线性代数从特定的因子的集合中分离出特定于cell的因素。
###### 2.7.3.6 Downsampling
修正库大小的最后一种方法是对表达式矩阵进行Downsampling，这样每个cell的总molecules大约是相同的。这种方法的好处是，零值将被向下的采样引入，从而消除由于检测到的基因数量不同而产生的任何偏差。然而，主要的缺点是过程不是确定的，所以每次向下抽样运行时，结果表达式矩阵略有不同。因此，经常分析必须在多个向下取样上运行，以确保结果是健壮的。
```
Down_Sample_Matrix <-
function (expr_mat) 
{
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
        prob <- min_lib_size/sum(x)
        return(unlist(lapply(x, function(y) {
            rbinom(1, y, prob)
        })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
}
```
##### 2.7.4 Effectiveness
为了比较不同的标准化方法的效率，我们将使用PCA图的视觉检查，通过scater的plotRLE（）函数来计算细胞相关的相对对数表达式。也就是说，有许多（很少）的细胞比大多数基因的中位表达要高（更低），从而在细胞中产生一个正的（负的）RLE，而标准化细胞的RLE接近于零。

```
calc_cell_RLE <-
function (expr_mat, spikes = NULL) 
{
    RLE_gene <- function(x) {
        if (median(unlist(x)) > 0) {
            log((x + 1)/(median(unlist(x)) + 1))/log(2)
        }
        else {
            rep(NA, times = length(x))
        }
    }
    if (!is.null(spikes)) {
        RLE_matrix <- t(apply(expr_mat[-spikes, ], 1, RLE_gene))
    }
    else {
        RLE_matrix <- t(apply(expr_mat, 1, RLE_gene))
    }
    cell_RLE <- apply(RLE_matrix, 2, median, na.rm = T)
    return(cell_RLE)
}
```
请注意，RLE、TMM和UQ大小因子方法是为bulkRNA-seq数据开发的，根据实验环境的不同，可能不适合single-cell RNA-seq数据，因为它们的基本假设可能会受到问题的违反。
注意，scater是edgeR的calcNormFactors函数的包装器，它实现了几个库大小的标准化方法，使我们可以很容易地将这些方法应用到我们的数据中。
#### 2.8 Nomalization practice（UMI）

```
install.packages("devtools")
devtools::install_github("hemberg-lab/scRNA.seq.funcs")
library(scRNA.seq.funcs)
library(scater)
library(scran)
options(stringsAsFactors = FALSE)
set.seed(1234567)
umi <- readRDS("tung/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control
```
##### 2.8.1 Raw

```
plotPCA(
    umi.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/20-exprs-norm_files/figure-html/norm-pca-raw-1.png)
##### 2.8.2 CPM

```
logcounts(umi.qc) <- log2(calculateCPM(umi.qc, use.size.factors = FALSE) + 1)
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/20-exprs-norm_files/figure-html/norm-pca-cpm-1.png)

```
plotRLE(
    umi.qc[endog_genes, ], 
    exprs_mats = list(Raw = "logcounts_raw", CPM = "logcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "batch"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/20-exprs-norm_files/figure-html/norm-ours-rle-cpm-1.png)
##### 2.8.3 Size-factor（RLE）

```
umi.qc <- normaliseExprs(
    umi.qc,
    method = "RLE", 
    feature_set = endog_genes,
    return_log = TRUE,
    return_norm_as_exprs = TRUE
)
```

```
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/20-exprs-norm_files/figure-html/norm-pca-rle-1.png)

```
plotRLE(
    umi.qc[endog_genes, ], 
    exprs_mats = list(Raw = "logcounts_raw", RLE = "logcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "batch"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/20-exprs-norm_files/figure-html/norm-ours-rle-rle-1.png)
##### 2.8.4 Upperquantile

```
umi.qc <- normaliseExprs(
    umi.qc,
    method = "upperquartile", 
    feature_set = endog_genes,
    p = 0.99,
    return_log = TRUE,
    return_norm_as_exprs = TRUE
)
```

```
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/20-exprs-norm_files/figure-html/norm-pca-uq-1.png)

```
plotRLE(
    umi.qc[endog_genes, ], 
    exprs_mats = list(Raw = "logcounts_raw", UQ = "logcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "batch"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/20-exprs-norm_files/figure-html/norm-ours-rle-uq-1.png)
##### TMM

```
umi.qc <- normaliseExprs(
    umi.qc,
    method = "TMM",
    feature_set = endog_genes,
    return_log = TRUE,
    return_norm_as_exprs = TRUE
)
```

```
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/20-exprs-norm_files/figure-html/norm-pca-tmm-1.png)

```
plotRLE(
    umi.qc[endog_genes, ], 
    exprs_mats = list(Raw = "logcounts_raw", TMM = "logcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "batch"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/20-exprs-norm_files/figure-html/norm-ours-rle-tmm-1.png)
##### 2.8.6 scran

```
qclust <- quickCluster(umi.qc, min.size = 30)
umi.qc <- computeSumFactors(umi.qc, sizes = 15, clusters = qclust)
umi.qc <- normalize(umi.qc)
```

```
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/20-exprs-norm_files/figure-html/norm-pca-lsf-1.png)

```
plotRLE(
    umi.qc[endog_genes, ], 
    exprs_mats = list(Raw = "logcounts_raw", scran = "logcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "batch"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/20-exprs-norm_files/figure-html/norm-ours-rle-scran-1.png)
scran有时会计算出负或零大小的因素。这些将完全扭曲标准化的表达式矩阵。我们可以检查scran的大小因素，如：

```
summary(sizeFactors(umi.qc))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.4646  0.7768  0.9562  1.0000  1.1444  3.4348
```
对于这个数据集，所有大小的因素都是合理的，所以我们已经完成了。如果您发现scran已经计算出了负大小的因素，那么尝试增加cluster和pool的大小，直到它们都是正的。
##### 2.8.7 Downsampling

```
logcounts(umi.qc) <- log2(Down_Sample_Matrix(counts(umi.qc)) + 1)
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/20-exprs-norm_files/figure-html/norm-pca-downsample-1.png)

```
plotRLE(
    umi.qc[endog_genes, ], 
    exprs_mats = list(Raw = "logcounts_raw", DownSample = "logcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "batch"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/20-exprs-norm_files/figure-html/norm-ours-rle-downsample-1.png)
##### 2.8.8 Normalisation for gene/transcript length
有些方法结合了库的大小和片段/基因长度的标准化，例如：
- RPKM：每千碱基的read（单端测序）
- FPKM：每千碱基的Fragment（与RPKM相同但是是用于双端测序的数据，确保映射到同一片段的成对结束不会被计数两次）
- TPM-每千碱基的Transcript（与RPKM相同，但是，标准化的顺序是颠倒的——长度优先，顺序深度第二）
这些方法不适用于我们的数据集，因为包含了UMI的文字记录的末尾被优先测序。此外，一般情况下，这些数据只能使用来自对齐的BAM文件的适当的量化软件，而不是从读取计数中计算出来，因为通常只有一部分的整个遗传/转录本是被测序的，而不是整个长度。如果有疑问，检查基因/转录长度和表达水平之间的关系。
然而，这里我们展示了如何使用scater来计算这些normalisations。首先，我们需要在kilobase中找到有效的转录长度。然而，我们的数据集只包含了基因id，因此我们将使用基因长度而不是转录本。scater使用biomaRt，它允许通过其他属性对基因进行注释：

```
umi.qc <- getBMFeatureAnnos(
    umi.qc,
    filters = "ensembl_gene_id", 
    attributes = c(
        "ensembl_gene_id",
        "hgnc_symbol",
        "chromosome_name",
        "start_position",
        "end_position"
    ), 
    feature_symbol = "hgnc_symbol",
    feature_id = "ensembl_gene_id",
    biomart = "ENSEMBL_MART_ENSEMBL", 
    dataset = "hsapiens_gene_ensembl",
    host = "www.ensembl.org"
)

# If you have mouse data, change the arguments based on this example:
# getBMFeatureAnnos(
#     object,
#     filters = "ensembl_transcript_id",
#     attributes = c(
#         "ensembl_transcript_id",
#         "ensembl_gene_id", 
#         "mgi_symbol",
#         "chromosome_name",
#         "transcript_biotype",
#         "transcript_start",
#         "transcript_end",
#         "transcript_count"
#     ),
#     feature_symbol = "mgi_symbol",
#     feature_id = "ensembl_gene_id",
#     biomart = "ENSEMBL_MART_ENSEMBL",
#     dataset = "mmusculus_gene_ensembl",
#     host = "www.ensembl.org"
# )
```
有些基因没有被注释，因此我们把它们过滤掉：
```
umi.qc.ann <- umi.qc[!is.na(rowData(umi.qc)$ensembl_gene_id), ]
```
现在我们通过使用endposition和startposition来计算kilobase的总基因长度：
```
eff_length <- 
    abs(rowData(umi.qc.ann)$end_position - rowData(umi.qc.ann)$start_position) / 1000
```

```
plot(eff_length, rowMeans(counts(umi.qc.ann)))
```
基因长度和平均表达之间没有关系，所以FPKM s和TPM对于这个数据集是不合适的。但无论如何我们都会证明他们的。
这里要注意的是计算总基因长度而不是exon总长度。许多基因会包含大量的内含子，所以它们的长度将与我们计算的结果大不相同。请把我们的计算看作是近似值。如果你想使用总exon长度，请参考[这一页](https://www.biostars.org/p/83901/)。
现在，我们已经准备好了进行规范：

```
tpm(umi.qc.ann) <- log2(calculateTPM(umi.qc.ann, eff_length) + 1)
```

```
plotPCA(
    umi.qc.ann,
    exprs_values = "tpm",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```

```
tpm(umi.qc.ann) <- log2(calculateFPKM(umi.qc.ann, eff_length) + 1)
```

```
plotPCA(
    umi.qc.ann,
    exprs_values = "tpm",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```
PCA寻找细胞之间的差异。每个基因的细胞的基因长度是相同的，因此FPKM几乎与CPM图相同（它只是旋转），因为它首先执行CPM然后标准化基因长度。然而，TPM是不同的，因为它在执行CPM之前根据它们的长度对基因进行加权。
#### 2.10 Dealing with confounders
##### 2.10.1 Introduction
在上一章中，我们对库大小进行了规范化，有效地将其作为混淆器删除。现在我们将考虑从我们的数据中删除其他不太明确的混淆器。技术混杂因素（又称批次效应）可能源于试剂，分离方法，进行实验的实验室/实验者，甚至实验的实施日期/时间的差异。特别是技术混杂因素和批次效应的会计是一个涉及实验设计原则的大课题。在这里，我们讨论了在实验设计合适时可以考虑混杂因素的方法。

从根本上说，技术混杂因素的解释涉及识别并且理想地去除表达数据中与感兴趣的生物信号无关（即混淆）的变异源。存在各种方法，其中一些使用加入或管家基因，其中一些使用内源基因。
###### 2.10.1.1 Advantages and disadvantages of using spike-ins to remove confounders
使用加标作为对照基因是有吸引力的，因为在我们的实验中向每个细胞添加了相同量的ERCC（或其他）spike-in。原则上，我们观察到的这些基因的所有变异都是由技术噪声引起的;而内源基因受技术噪声和生物变异性的影响。通过将模型拟合到spike-in内并从内源基因“减去”这种噪声，可以消除技术噪声。基于这个前提有几种方法可用（例如BASiCS，scLVM，RUVg）;每个使用不同的噪声模型和不同的拟合程序。或者，可以鉴定表现出超出技术噪声的显着变异的基因（例如，与中值的距离，高度可变的基因）。然而，使用加标用于标准化存在问题（特别是ERCC，源自细菌序列），包括由于各种原因，它们的变异性实际上可能高于内源基因的变异性。
鉴于使用spike-ins的问题，通常可以通过使用内源基因获得更好的结果。在我们有大量内源基因的情况下，平均而言，细胞之间并没有系统地变化，我们期望技术效应影响大量基因（这是一个非常常见和合理的假设），那么这些方法（例如， RUVs方法）可以很好地执行。

```
library(scRNA.seq.funcs)
library(RUVSeq)
library(scater)
library(SingleCellExperiment)
library(scran)
library(kBET)
library(sva) # Combat
library(edgeR)
set.seed(1234567)
options(stringsAsFactors = FALSE)
umi <- readRDS("tung/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control
erccs <- rowData(umi.qc)$is_feature_control

qclust <- quickCluster(umi.qc, min.size = 30)
umi.qc <- computeSumFactors(umi.qc, sizes = 15, clusters = qclust)
umi.qc <- normalize(umi.qc)
```
##### 2.10.2 Remove Unwanted Variation
导致技术噪声的因素经常表现为“批次效应”，其中在不同日期或不同技术人员处理的细胞系统地彼此不同。通常可以使用相同的工具或轻微的变体来消除技术噪音并校正批量效应。我们将考虑Remove Unwanted Variation（RUVSeq）。简而言之，RUVSeq的工作原理如下。对于ñ样品和Ĵ基因，考虑以下广义线性模型（GLM），其中RNA-Seq读数计数在已知的感兴趣的协变量和不想要的变异的未知因素上回归：logE[Y|W,X,O]=Wα+Xβ+O

这里，ÿ是ñ×Ĵ观察到的基因水平读数计数矩阵，w 是一个ñ×ķ对应于“不想要的变化”因素的矩阵0是ñ×Ĵ偏移矩阵可以设置为零或使用其他一些归一化过程（例如上四分位标准化）进行估计。同时评估W，α，β,k是不可行的。对于一个给定的
k，下面的三种方法来估计不需要的变化的因素W
:
- RUVg使用阴性对照基因（例如ERCC），假定在样品中具有恒定表达;
- RUVs使用中心（技术）重复/阴性对照样品，其中感兴趣的协变量是恒定的;
- RUVr使用残差，例如来自对感兴趣的协变量的计数的第一遍GLM回归。
我们将集中讨论前两种方法。
###### 2.10.1.1 RUVg

```
ruvg <- RUVg(counts(umi.qc), erccs, k = 1)
assay(umi.qc, "ruvg1") <- log2(
    t(t(ruvg$normalizedCounts) / colSums(ruvg$normalizedCounts) * 1e6) + 1
)
ruvg <- RUVg(counts(umi.qc), erccs, k = 10)
assay(umi.qc, "ruvg10") <- log2(
    t(t(ruvg$normalizedCounts) / colSums(ruvg$normalizedCounts) * 1e6) + 1
)
```
###### 2.10.2.2

```
scIdx <- matrix(-1, ncol = max(table(umi.qc$individual)), nrow = 3)
tmp <- which(umi.qc$individual == "NA19098")
scIdx[1, 1:length(tmp)] <- tmp
tmp <- which(umi.qc$individual == "NA19101")
scIdx[2, 1:length(tmp)] <- tmp
tmp <- which(umi.qc$individual == "NA19239")
scIdx[3, 1:length(tmp)] <- tmp
cIdx <- rownames(umi.qc)
ruvs <- RUVs(counts(umi.qc), cIdx, k = 1, scIdx = scIdx, isLog = FALSE)
assay(umi.qc, "ruvs1") <- log2(
    t(t(ruvs$normalizedCounts) / colSums(ruvs$normalizedCounts) * 1e6) + 1
)
ruvs <- RUVs(counts(umi.qc), cIdx, k = 10, scIdx = scIdx, isLog = FALSE)
assay(umi.qc, "ruvs10") <- log2(
    t(t(ruvs$normalizedCounts) / colSums(ruvs$normalizedCounts) * 1e6) + 1
)
```
##### 2.10.3 Combat
如果您有平衡设计的实验，Combat可用于消除批量效应，同时通过使用mod参数指定生物效应来保留生物效应。然而，Tung数据包含多个实验重复而不是平衡设计，因此使用mod1来保持生物变异性将导致错误。

```
combat_data <- logcounts(umi.qc)
mod_data <- as.data.frame(t(combat_data))
# Basic batch removal
mod0 = model.matrix(~ 1, data = mod_data) 
# Preserve biological variability
mod1 = model.matrix(~ umi.qc$individual, data = mod_data) 
# adjust for total genes detected
mod2 = model.matrix(~ umi.qc$total_features, data = mod_data)
assay(umi.qc, "combat") <- ComBat(
    dat = t(mod_data), 
    batch = factor(umi.qc$batch), 
    mod = mod0,
    par.prior = TRUE,
    prior.plots = FALSE
)
```
##### 2.10.4 mnnCorrect
mnnCorrect（Haghverdi等人，2017）假设每批次彼此共享至少一种生物学条件。因此，它适用于各种平衡的实验设计。然而，Tung数据包含每个个体而非平衡批次的多个重复，因此我们将分别标准化每个个体。请注意，由于混淆的实验设计，这将消除同一个体中批次之间的批次效应，但不会消除不同个体中批次之间的批次效应。

因此，我们将合并来自每个个体的复制品以形成三个批次。

```
do_mnn <- function(data.qc) {
    batch1 <- logcounts(data.qc[, data.qc$replicate == "r1"])
    batch2 <- logcounts(data.qc[, data.qc$replicate == "r2"])
    batch3 <- logcounts(data.qc[, data.qc$replicate == "r3"])
    
    if (ncol(batch2) > 0) {
        x = mnnCorrect(
          batch1, batch2, batch3,  
          k = 20,
          sigma = 0.1,
          cos.norm.in = TRUE,
          svd.dim = 2
        )
        res1 <- x$corrected[[1]]
        res2 <- x$corrected[[2]]
        res3 <- x$corrected[[3]]
        dimnames(res1) <- dimnames(batch1)
        dimnames(res2) <- dimnames(batch2)
        dimnames(res3) <- dimnames(batch3)
        return(cbind(res1, res2, res3))
    } else {
        x = mnnCorrect(
          batch1, batch3,  
          k = 20,
          sigma = 0.1,
          cos.norm.in = TRUE,
          svd.dim = 2
        )
        res1 <- x$corrected[[1]]
        res3 <- x$corrected[[2]]
        dimnames(res1) <- dimnames(batch1)
        dimnames(res3) <- dimnames(batch3)
        return(cbind(res1, res3))
    }
}

indi1 <- do_mnn(umi.qc[, umi.qc$individual == "NA19098"])
indi2 <- do_mnn(umi.qc[, umi.qc$individual == "NA19101"])
indi3 <- do_mnn(umi.qc[, umi.qc$individual == "NA19239"])

assay(umi.qc, "mnn") <- cbind(indi1, indi2, indi3)

# For a balanced design: 
#assay(umi.qc, "mnn") <- mnnCorrect(
#    list(B1 = logcounts(batch1), B2 = logcounts(batch2), B3 = logcounts(batch3)),  
#    k = 20,
#    sigma = 0.1,
#    cos.norm = TRUE,
#    svd.dim = 2
#)
```
##### 2.10.5 GLM

一般线性模型是Combat的更简单版本。如果你有一个平衡的设计，它可以纠正批量，同时保留生物效应。在混淆/复制设计中，生物效应将不适合/保留。与mnnCorrect类似，我们可以分别从每个个体中删除批量效应，以保持个体之间的生物（和技术）差异。出于演示目的，我们将天真地纠正所有共同创建的批处理效果：

```
glm_fun <- function(g, batch, indi) {
  model <- glm(g ~ batch + indi)
  model$coef[1] <- 0 # replace intercept with 0 to preserve reference batch.
  return(model$coef)
}
effects <- apply(
    logcounts(umi.qc), 
    1, 
    glm_fun, 
    batch = umi.qc$batch, 
    indi = umi.qc$individual
)
corrected <- logcounts(umi.qc) - t(effects[as.numeric(factor(umi.qc$batch)), ])
assay(umi.qc, "glm") <- corrected
```
##### 2.10.6 How to evaluate and compare confounder removal strategies
在考虑去除混杂因素的不同方法时，一个关键问题是如何定量确定哪一种方法最有效。比较具有挑战性的主要原因是因为通常很难知道什么对应于技术counfounders和什么是有趣的生物变异。在这里，我们根据我们对实验设计的了解，考虑三种不同的指标。根据您希望解决的生物学问题，选择一个允许您评估可能是特定情况下最关注的混杂因素的指标非常重要。
###### 2.10.6.1 Effectiveness 1
我们通过检查PCA图来评估标准化的有效性，其中颜色对应于技术重复，并且形状对应于不同的生物样品（个体）。生物样品和散布批次的分离表明技术变异已被消除。我们总是使用log2-cpm标准化数据来匹配PCA的假设。

```
for(n in assayNames(umi.qc)) {
    print(
        plotPCA(
            umi.qc[endog_genes, ],
            colour_by = "batch",
            size_by = "total_features",
            shape_by = "individual",
            exprs_values = n
        ) +
        ggtitle(n)
    )
}
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-10-1.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-10-2.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-10-3.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-10-4.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-10-5.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-10-6.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-10-7.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-10-8.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-10-9.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-10-10.png)
![image]https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-10-11.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-10-12.png)
###### 2.10.6.2 Effectiveness 2
我们还可以使用跨cell的相对对数表达式（RLE）检查校正的有效性，以确认已从数据集中删除技术噪声。注意RLE仅评估每个细胞的高于和低于平均值的基因数量是否相等 - 即系统性技术效应。 RLE可能无法检测到批次之间的随机技术噪音。

```
res <- list()
for(n in assayNames(umi.qc)) {
    res[[n]] <- suppressWarnings(calc_cell_RLE(assay(umi.qc, n), erccs))
}
par(mar=c(6,4,1,1))
boxplot(res, las=2)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-11-1.png)
###### 2.10.6.3 Effectiveness 3
我们可以重复前面的分析，检查批处理效果是否被删除

```
for(n in assayNames(umi.qc)) {
    print(
        plotQC(
            umi.qc[endog_genes, ],
            type = "expl",
            exprs_values = n,
            variables = c(
                "total_features",
                "total_counts",
                "batch",
                "individual",
                "pct_counts_ERCC",
                "pct_counts_MT"
            )
        ) +
        ggtitle(n)
    )
}
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/confound-cpm-1.png)![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/confound-cpm-2.png)![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/confound-cpm-3.png)![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/confound-cpm-4.png)![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/confound-cpm-5.png)![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/confound-cpm-6.png)![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/confound-cpm-7.png)![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/confound-cpm-8.png)![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/confound-cpm-9.png)![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/confound-cpm-10.png)![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/confound-cpm-11.png)![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/confound-cpm-12.png)
###### 2.10.6.4 Effectiveness 4
检查批次效应校正功效的另一种方法是考虑在数据的局部子样本中来自不同批次的点的混合。如果没有批次效应，则任何局部区域中每批次的细胞比例应等于每批次中细胞的全局比例。

kBET（Buttner等人，2017）在随机细胞周围采用kNN网络，并根据二项分布测试每批次的细胞数量。这些测试的拒绝率表明数据中仍存在批次效应的严重性（高拒绝率=强批次效应）。 kBET假设每批包含相同的生物组补充，因此如果使用完美平衡的设计，它只能应用于整个数据集。但是，如果将kBET分别应用于每个生物组，也可以将其应用于复制数据。在Tung数据的情况下，我们将独立地将kBET应用于每个个体以检查残余批次效应。然而，该方法不能识别与生物条件混淆的残留批次效应。此外，kBET不确定是否保留了生物信号。

```
compare_kBET_results <- function(sce){
    indiv <- unique(sce$individual)
    norms <- assayNames(sce) # Get all normalizations
    results <- list()
    for (i in indiv){ 
        for (j in norms){
            tmp <- kBET(
                df = t(assay(sce[,sce$individual== i], j)), 
                batch = sce$batch[sce$individual==i], 
                heuristic = TRUE, 
                verbose = FALSE, 
                addTest = FALSE, 
                plot = FALSE)
            results[[i]][[j]] <- tmp$summary$kBET.observed[1]
        }
    }
    return(as.data.frame(results))
}

eff_debatching <- compare_kBET_results(umi.qc)
```

```
require("reshape2")
require("RColorBrewer")
# Plot results
dod <- melt(as.matrix(eff_debatching),  value.name = "kBET")
colnames(dod)[1:2] <- c("Normalisation", "Individual")

colorset <- c('gray', brewer.pal(n = 9, "RdYlBu"))

ggplot(dod, aes(Normalisation, Individual, fill=kBET)) +  
    geom_tile() +
    scale_fill_gradient2(
        na.value = "gray",
        low = colorset[2],
        mid=colorset[6],
        high = colorset[10],
        midpoint = 0.5, limit = c(0,1)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    theme(
        axis.text.x = element_text(
            angle = 45, 
            vjust = 1, 
            size = 12, 
            hjust = 1
        )
    ) + 
    ggtitle("Effect of batch regression methods per individual")
```

![image](https://hemberg-lab.github.io/scRNA.seq.course/22-remove-conf_files/figure-html/unnamed-chunk-13-1.png)
##### 2.11 Dealing with confounders（Reads）

```
library(scRNA.seq.funcs)
library(RUVSeq)
library(scater)
library(SingleCellExperiment)
library(scran)
library(kBET)
library(sva) # Combat
library(edgeR)
set.seed(1234567)
options(stringsAsFactors = FALSE)
reads <- readRDS("tung/reads.rds")
reads.qc <- reads[rowData(reads)$use, colData(reads)$use]
endog_genes <- !rowData(reads.qc)$is_feature_control
erccs <- rowData(reads.qc)$is_feature_control

qclust <- quickCluster(reads.qc, min.size = 30)
reads.qc <- computeSumFactors(reads.qc, sizes = 15, clusters = qclust)
reads.qc <- normalize(reads.qc)
```

```
ruvg <- RUVg(counts(reads.qc), erccs, k = 1)
assay(reads.qc, "ruvg1") <- log2(
    t(t(ruvg$normalizedCounts) / colSums(ruvg$normalizedCounts) * 1e6) + 1
)
ruvg <- RUVg(counts(reads.qc), erccs, k = 10)
assay(reads.qc, "ruvg10") <- log2(
    t(t(ruvg$normalizedCounts) / colSums(ruvg$normalizedCounts) * 1e6) + 1
)
```

```
scIdx <- matrix(-1, ncol = max(table(reads.qc$individual)), nrow = 3)
tmp <- which(reads.qc$individual == "NA19098")
scIdx[1, 1:length(tmp)] <- tmp
tmp <- which(reads.qc$individual == "NA19101")
scIdx[2, 1:length(tmp)] <- tmp
tmp <- which(reads.qc$individual == "NA19239")
scIdx[3, 1:length(tmp)] <- tmp
cIdx <- rownames(reads.qc)
ruvs <- RUVs(counts(reads.qc), cIdx, k = 1, scIdx = scIdx, isLog = FALSE)
assay(reads.qc, "ruvs1") <- log2(
    t(t(ruvs$normalizedCounts) / colSums(ruvs$normalizedCounts) * 1e6) + 1
)
ruvs <- RUVs(counts(reads.qc), cIdx, k = 10, scIdx = scIdx, isLog = FALSE)
assay(reads.qc, "ruvs10") <- log2(
    t(t(ruvs$normalizedCounts) / colSums(ruvs$normalizedCounts) * 1e6) + 1
)
```

```
combat_data <- logcounts(reads.qc)
mod_data <- as.data.frame(t(combat_data))
# Basic batch removal
mod0 = model.matrix(~ 1, data = mod_data) 
# Preserve biological variability
mod1 = model.matrix(~ reads.qc$individual, data = mod_data) 
# adjust for total genes detected
mod2 = model.matrix(~ reads.qc$total_features, data = mod_data)
assay(reads.qc, "combat") <- ComBat(
    dat = t(mod_data), 
    batch = factor(reads.qc$batch), 
    mod = mod0,
    par.prior = TRUE,
    prior.plots = FALSE
)
```

```
## Standardizing Data across genes
```

```
do_mnn <- function(data.qc) {
    batch1 <- logcounts(data.qc[, data.qc$replicate == "r1"])
    batch2 <- logcounts(data.qc[, data.qc$replicate == "r2"])
    batch3 <- logcounts(data.qc[, data.qc$replicate == "r3"])
    
    if (ncol(batch2) > 0) {
        x = mnnCorrect(
          batch1, batch2, batch3,  
          k = 20,
          sigma = 0.1,
          cos.norm.in = TRUE,
          svd.dim = 2
        )
        res1 <- x$corrected[[1]]
        res2 <- x$corrected[[2]]
        res3 <- x$corrected[[3]]
        dimnames(res1) <- dimnames(batch1)
        dimnames(res2) <- dimnames(batch2)
        dimnames(res3) <- dimnames(batch3)
        return(cbind(res1, res2, res3))
    } else {
        x = mnnCorrect(
          batch1, batch3,  
          k = 20,
          sigma = 0.1,
          cos.norm.in = TRUE,
          svd.dim = 2
        )
        res1 <- x$corrected[[1]]
        res3 <- x$corrected[[2]]
        dimnames(res1) <- dimnames(batch1)
        dimnames(res3) <- dimnames(batch3)
        return(cbind(res1, res3))
    }
}

indi1 <- do_mnn(reads.qc[, reads.qc$individual == "NA19098"])
indi2 <- do_mnn(reads.qc[, reads.qc$individual == "NA19101"])
indi3 <- do_mnn(reads.qc[, reads.qc$individual == "NA19239"])

assay(reads.qc, "mnn") <- cbind(indi1, indi2, indi3)

# For a balanced design: 
#assay(reads.qc, "mnn") <- mnnCorrect(
#    list(B1 = logcounts(batch1), B2 = logcounts(batch2), B3 = logcounts(batch3)),  
#    k = 20,
#    sigma = 0.1,
#    cos.norm = TRUE,
#    svd.dim = 2
#)
```

```
glm_fun <- function(g, batch, indi) {
  model <- glm(g ~ batch + indi)
  model$coef[1] <- 0 # replace intercept with 0 to preserve reference batch.
  return(model$coef)
}
effects <- apply(
    logcounts(reads.qc), 
    1, 
    glm_fun, 
    batch = reads.qc$batch, 
    indi = reads.qc$individual
)
corrected <- logcounts(reads.qc) - t(effects[as.numeric(factor(reads.qc$batch)), ])
assay(reads.qc, "glm") <- corrected
```

```
for(n in assayNames(reads.qc)) {
    print(
        plotPCA(
            reads.qc[endog_genes, ],
            colour_by = "batch",
            size_by = "total_features",
            shape_by = "individual",
            exprs_values = n
        ) +
        ggtitle(n)
    )
}
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-10-1.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-10-2.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-10-3.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-10-4.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-10-5.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-10-6.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-10-7.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-10-8.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-10-9.png)![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-10-10.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-10-11.png)![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-10-12.png)

```
res <- list()
for(n in assayNames(reads.qc)) {
    res[[n]] <- suppressWarnings(calc_cell_RLE(assay(reads.qc, n), erccs))
}
par(mar=c(6,4,1,1))
boxplot(res, las=2)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-11-1.png)
```
for(n in assayNames(reads.qc)) {
    print(
        plotQC(
            reads.qc[endog_genes, ],
            type = "expl",
            exprs_values = n,
            variables = c(
                "total_features",
                "total_counts",
                "batch",
                "individual",
                "pct_counts_ERCC",
                "pct_counts_MT"
            )
        ) +
        ggtitle(n)
    )
}
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-12-1.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-12-2.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-12-3.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-12-4.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-12-5.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-12-6.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-12-7.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-12-8.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-12-9.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-12-10.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-12-11.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-12-12.png)

```
compare_kBET_results <- function(sce){
    indiv <- unique(sce$individual)
    norms <- assayNames(sce) # Get all normalizations
    results <- list()
    for (i in indiv){ 
        for (j in norms){
            tmp <- kBET(
                df = t(assay(sce[,sce$individual== i], j)), 
                batch = sce$batch[sce$individual==i], 
                heuristic = TRUE, 
                verbose = FALSE, 
                addTest = FALSE, 
                plot = FALSE)
            results[[i]][[j]] <- tmp$summary$kBET.observed[1]
        }
    }
    return(as.data.frame(results))
}

eff_debatching <- compare_kBET_results(reads.qc)
```

```
require("reshape2")
require("RColorBrewer")
# Plot results
dod <- melt(as.matrix(eff_debatching),  value.name = "kBET")
colnames(dod)[1:2] <- c("Normalisation", "Individual")

colorset <- c('gray', brewer.pal(n = 9, "RdYlBu"))

ggplot(dod, aes(Normalisation, Individual, fill=kBET)) +  
    geom_tile() +
    scale_fill_gradient2(
        na.value = "gray",
        low = colorset[2],
        mid=colorset[6],
        high = colorset[10],
        midpoint = 0.5, limit = c(0,1)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    theme(
        axis.text.x = element_text(
            angle = 45, 
            vjust = 1, 
            size = 12, 
            hjust = 1
        )
    ) + 
    ggtitle("Effect of batch regression methods per individual")
```

![image](https://hemberg-lab.github.io/scRNA.seq.course/23-remove-conf-reads_files/figure-html/unnamed-chunk-14-1.png)
### 3 Biological Analysis
#### 3.1 Clustering Introduction
一旦我们对数据进行了标准化并删除了混淆因素，我们就可以进行与手头的生物学问题相关的分析。分析的确切性质取决于数据集。然而，有一些方面在广泛的背景下是有用的，我们将在接下来的几章中讨论其中的一些方面。我们将从scRNA-seq数据的聚类开始。
##### 3.1.1 Introduction
scRNA-seq最有希望的应用之一是基于转录谱的de novo discovery和细胞类型的注释。在计算上，这是一个难题，因为它相当于无监督的聚类。也就是说，我们需要根据转录组的相似性识别细胞群，而无需事先了解标记。此外，在大多数情况下，我们甚至不知道先验集群的数量。由于高水平的噪音（技术和生物）和维度，问题变得更具挑战性。
##### 3.1.2 Dimensionality reductions
在处理大型数据集时，应用某种降维方法通常是有益的。通过将数据投影到较低维度的子空间上，人们通常能够显着减少噪声量。另一个好处是，在2维或3维子空间中可视化数据通常要容易得多。我们已经讨论过PCA和t-SNE。
##### 3.1.3 Clustering methods
无监督聚类在许多不同的应用中是有用的，并且已经在机器学习中被广泛研究。一些最流行的方法是层次聚类，k均值聚类和基于图的聚类。
###### 3.1.3.1 Hierarchical clustering（层次聚类）
在分层聚类中，可以使用自下而上或自上而下的方法。在前一种情况下，每个单元最初被分配到它自己的集群，随后合并成对的集合以创建一个hieararchy：
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/hierarchical_clustering1.png)  
原始的数据
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/hierarchical_clustering2.png)
使用自上而下策略，首先从一个集群中的所有观察开始，然后递归地拆分每个集群以形成层次结构。该策略的一个优点是该方法是确定性的。
###### 3.1.3.2 k-means(k均值聚类)
在k均值聚类中，目标是将N个单元划分为k个不同的聚类。以迭代方式，分配集群中心，并将每个单元分配给最近的集群：
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/k-means.png)
用于scRNA-seq分析的大多数方法在某些时候包括k均值步骤。
###### 3.1.3.3 Graph-based methods（基于图的聚类）
在过去的二十年中，人们对分析各个领域的网络产生了很大的兴趣。一个目标是识别网络中节点的组或模块。
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/graph_network.jpg)
通过构建图表，其中一些方法可以应用于scRNA-seq数据，其中每个节点代表一个细胞。请注意，构造图形并为边缘指定权重并非易事。基于图的方法的一个优点是它们中的一些非常有效并且可以应用于包含数百万个节点的网络。
##### 3.1.4 Challenge in clustering
- What is the number of clusters k?
- What is a cell type?
- Scalability: in the last few years the number of cells in scRNA-seq experiments has grown by several orders of magnitude 

```math
from ~ 
10^2to ~ 10^6
```
- Tools are not user-friendly
##### 3.1.5 Tools for scRNA-seq data
###### 3.1.5.1 SINCERA
- SINCERA（Guo等人，2015）基于层次聚类
- 在聚类之前将数据转换为z分数
- 通过查找层次结构中的第一个单例群集来识别k

###### 3.1.5.2 pcaReduce
pcaReduce (žurauskienė and Yau 2016) combines PCA, k-means and “iterative” hierarchical clustering. Starting from a large number of clusters pcaReduce iteratively merges similar clusters; after each merging event it removes the principle component explaning the least variance in the data.
###### 3.1.5.3 SC3
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/sc3.png)

- SC3（Kiselev等人，2017）基于PCA和光谱维数减少
- 利用k-means
- 另外执行一致性聚类
###### 3.1.5.4 tSNE+k-means
- 基于tSEN
- 利用k-means
###### 3.1.5.5 SNN-Cliq
SNN-Cliq（C. Xu和Su 2015）是一种基于图的方法。首先，该方法根据距离测量识别每个cell的k-最近邻居。这用于计算每对cell之间的共享最近邻居（SNN）的数量。如果它们具有至少一个SNN，则通过在两个cell之间放置边来构建图。集群被定义为使用“clique”方法在它们之间具有许多边缘的cell组。SNN-Cliq需要手动定义几个参数。
###### 3.1.5.6 Seurat clustering
Seurat聚类基于类似于SNN-Cliq的社区检测方法以及之前提出的用于分析CyTOF数据的方法（Levine等人2015）。由于Seurat变得更像是scRNA-seq数据分析的一体化工具，我们专门用一章来更详细地讨论它。
##### 3.1.6 Comparing clustering

为了比较两组聚类标签，我们可以使用调整后的Rand索引。该指数衡量两个数据集群之间的相似性。调整后的兰德指数的值在于[0;1]间隔，在哪里1意味着两个聚类是相同的0表示偶然预期的相似程度。
#### 3.2 Clustering example

```
library(pcaMethods)
library(pcaReduce)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
set.seed(1234567)
```

为了说明scRNA-seq数据的聚类，我们考虑了来自发育小鼠胚胎的细胞的Deng数据集（Deng等人，2014）。我们已经预处理了数据集并提前创建了一个SingleCellExperiment对象。我们还使用原始出版物中标识的cell类型注释了cell（它是colData中的cell_type2列）。
##### 3.2.1 Deng dataset
载入数据

```
deng <- readRDS("deng/deng-reads.rds")
deng
```

```
## class: SingleCellExperiment 
## dim: 22431 268 
## metadata(0):
## assays(2): counts logcounts
## rownames(22431): Hvcn1 Gbp7 ... Sox5 Alg11
## rowData names(10): feature_symbol is_feature_control ...
##   total_counts log10_total_counts
## colnames(268): 16cell 16cell.1 ... zy.2 zy.3
## colData names(30): cell_type2 cell_type1 ... pct_counts_ERCC
##   is_cell_control
## reducedDimNames(0):
## spikeNames(1): ERCC
```
查看细胞类型注释：

```
table(colData(deng)$cell_type2)
```

```
## 
##     16cell      4cell      8cell early2cell earlyblast  late2cell 
##         50         14         37          8         43         10 
##  lateblast   mid2cell   midblast         zy 
##         30         12         60          4
```
一个简单的PCA分析已经分离了一些强大的细胞类型，并提供了数据结构的一些见解：

```
plotPCA(deng, colour_by = "cell_type2")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/25-clustering_files/figure-html/unnamed-chunk-5-1.png)
正如你所看到的，早期的细胞类型分离得很好，但是三个胚泡的时间点更难以区分。
##### 3.2.2 SC3
让我们对deng数据进行SC3聚类。SC3的优点是它可以直接摄取SingleCellExperiment对象。
现在让我们想象一下，我们不知道集群k（cell类型）的数量。SC3可以为您估计一些集群：

```
deng <- sc3_estimate_k(deng)
## Estimating k...
metadata(deng)$sc3$k_estimation
## [1] 6
```

有趣的是，SC3预测的细胞类型数量小于原始数据注释。然而，在不同细胞类型的早期，中期和晚期阶段，我们将恰好有6种细胞类型。我们将合并的单元格类型存储在colData的cell_type1列中：

```
plotPCA(deng, colour_by = "cell_type1")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/25-clustering_files/figure-html/unnamed-chunk-7-1.png)
现在我们已经准备好运行SC3了（我们还要求它计算集群的生物属性）：

```
deng <- sc3(deng, ks = 10, biology = TRUE)
```

SC3的结果由几个不同的输出组成，这里我们展示了其中的一些：
共识度矩阵（Consensus matrix）:
![image](https://hemberg-lab.github.io/scRNA.seq.course/25-clustering_files/figure-html/unnamed-chunk-9-1.png)
轮廓图（silhouette plot）：Silhouette是指在数据集群内解释和验证一致性的方法。该技术提供了一个简洁的图形表示，表明每个对象在其集群中的位置。
轮廓值是对象与其自身群集（内聚力）相比与其他群集（分离）相似程度的度量。轮廓范围从-1到+1，其中高值表示对象与其自己的簇很好地匹配并且与相邻簇不匹配。如果大多数对象具有高值，则群集配置是合适的。如果许多点具有低值或负值，则群集配置可能具有太多或太少的群集。

```
sc3_plot_silhouette(deng, k = 10)
```

![image](https://hemberg-lab.github.io/scRNA.seq.course/25-clustering_files/figure-html/unnamed-chunk-10-1.png)  
表达矩阵的热图（Heatmap of the expression matrix:）：

```
sc3_plot_expression(deng, k = 10, show_pdata = "cell_type2")
```

![image](https://hemberg-lab.github.io/scRNA.seq.course/25-clustering_files/figure-html/unnamed-chunk-11-1.png)
确定marker基因：

```
sc3_plot_markers(deng, k = 10, show_pdata = "cell_type2")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/25-clustering_files/figure-html/unnamed-chunk-12-1.png)
PCA图，突出显示的SC3集群：

```
plotPCA(deng, colour_by = "sc3_10_clusters")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/25-clustering_files/figure-html/unnamed-chunk-13-1.png)
将SC3集群的结果与原始的cell类型标签进行比较：

```
adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$sc3_10_clusters)
```

```
## [1] 0.7705208
```
兰德指数（Rand index）在统计学中，特别是在数据聚类中，是两个数据聚类之间相似性的度量。可以定义Rand索引的形式，其针对元素的机会分组进行调整，这是调整后的Rand索引。从数学的角度来看，Rand指数与准确性有关，但即使不使用类标签也适用。
SC3也可以在一个交互式的Shiny会话中运行：

```
sc3_interactive(deng)
```

由于直接计算距离，SC3在细胞数量>5000时变得非常缓慢,对于包含的大型数据集达 
```math
10^5
```
我们用Seurat重新评论的细胞
##### 3.2.3 pcaReduce
pcaReduce直接在表达式矩阵上运行。建议在运行pcaReduce之前使用gene filter和log transfoemation。我们将使用默认的SC3基因过滤器（请注意，默认情况下，scater对象的exprs会进行对数转换）。

```
# use the same gene filter as in SC3
input <- logcounts(deng[rowData(deng)$sc3_gene_filter, ])
```

pcaReduce使用了几个参数：
* nbt定义了pcaReduce运行的次数（它是随机的，在不同的运行后可能有不同的解决方案）
* q定义了开始聚类的维数。输出将包含所有分区
ķ从2到q + 1。
* method定义了用于聚类的方法。
S- 执行基于sampleing的合并，
M- 基于最大概率执行合并。
我们运行pcaReduce一次：

```
# run pcaReduce 1 time creating hierarchies from 1 to 30 clusters
pca.red <- PCAreduce(t(input), nbt = 1, q = 30, method = 'S')[[1]]
```

```
colData(deng)$pcaReduce <- as.character(pca.red[,32 - 10])
plotPCA(deng, colour_by = "pcaReduce")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/25-clustering_files/figure-html/unnamed-chunk-18-1.png)

```
colData(deng)$pcaReducek2 <- as.character(pca.red[,32 - 2])
plotPCA(deng, colour_by = "pcaReducek2")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/25-clustering_files/figure-html/clust-pca-reduce2-1.png)
##### 3.2.4 tSEN+kmeans
我们之前看到的tSNE图，我们使用scater是由Rtsne和ggplot2包构建的。在这里我们将做同样的事情：

```
deng <- plotTSNE(deng, rand_seed = 1, return_SCE = TRUE)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/25-clustering_files/figure-html/clust-tsne-1.png)
请注意，上图中的所有点都是黑色的。这与我们之前看到的不同，当时细胞基于注释着色。这里我们没有任何注释，所有单元格都来自同一批次，因此所有点都是黑色的。
现在我们将k-means聚类算法应用于tSNE图上的点。你在中看到了多少组？

我们将从ķ=8开始

```
colData(deng)$tSNE_kmeans <- as.character(kmeans(deng@reducedDims$TSNE, centers = 8)$clust)
plotTSNE(deng, rand_seed = 1, colour_by = "tSNE_kmeans")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/25-clustering_files/figure-html/clust-tsne-kmeans2-1.png)
您可能已经注意到，pcaReduce和tSNE + kmeans都是随机的，每次运行时都会得到不同的结果。为了更好地了解解决方案，我们需要多次运行这些方法。 SC3也是随机的，但由于共识步骤，它更稳健，不太可能产生不同的结果。
##### 3.2.5 SNN-Clip
使用作者提供的SNN-cliq参数来运行：

```
distan <- "euclidean"
par.k <- 3
par.r <- 0.7
par.m <- 0.5
# construct a graph
scRNA.seq.funcs::SNN(
    data = t(input),
    outfile = "snn-cliq.txt",
    k = par.k,
    distance = distan
)
# find clusters in the graph
snn.res <- 
    system(
        paste0(
            "python utils/Cliq.py ", 
            "-i snn-cliq.txt ",
            "-o res-snn-cliq.txt ",
            "-r ", par.r,
            " -m ", par.m
        ),
        intern = TRUE
    )
cat(paste(snn.res, collapse = "\n"))
```

```
## input file snn-cliq.txt
## find 66 quasi-cliques
## merged into 29 clusters
## unique assign done
```

```
snn.res <- read.table("res-snn-cliq.txt")
# remove files that were created during the analysis
system("rm snn-cliq.txt res-snn-cliq.txt")

colData(deng)$SNNCliq <- as.character(snn.res[,1])
plotPCA(deng, colour_by = "SNNCliq")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/25-clustering_files/figure-html/unnamed-chunk-21-1.png)
##### 3.2.6 SINCREA
如前一章所述，SINCERA基于层次聚类。要记住的一件重要事情是它在进行聚类之前执行基因水平的z分数转换：

```
# perform gene-by-gene per-sample z-score transformation
dat <- apply(input, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
# hierarchical clustering
dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
hc <- hclust(dd, method = "average")
```
如果群集的数量未知，则SINCERA可以将k识别为生成不超过指定数量的单个群集（仅包含1个小区的群集）的分层树的最小高度

```
num.singleton <- 0
kk <- 1
for (i in 2:dim(dat)[2]) {
    clusters <- cutree(hc, k = i)
    clustersizes <- as.data.frame(table(clusters))
    singleton.clusters <- which(clustersizes$Freq < 2)
    if (length(singleton.clusters) <= num.singleton) {
        kk <- i
    } else {
        break;
    }
}
cat(kk)
```
查看SINCERA结果做的热图：


```
pheatmap(
    t(dat),
    cluster_cols = hc,
    cutree_cols = kk,
    kmeans_k = 100,
    show_rownames = FALSE
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/25-clustering_files/figure-html/clust-sincera-1.png)
#### 3.3 Feature Selection

```
library(scRNA.seq.funcs)
library(matrixStats)
library(M3Drop)
library(RColorBrewer)
library(SingleCellExperiment)
set.seed(1)
```
单细胞RNASeq能够测量每个细胞中数千个基因的表达。然而，在大多数情况下，只有一部分将显示对感兴趣的生物学状况的反应，例如，细胞类型的差异，分化的驱动因素，对环境刺激的反应。由于技术噪音，在scRNASeq实验中检测到的大多数基因仅在不同水平检测到。这样做的一个结果是技术噪声和批量效应可能使感兴趣的生物信号模糊不清。

因此，通常有利的是进行特征选择以去除那些仅表现出来自下游分析的技术噪音的基因。这通常不仅会增加数据中的信噪比;它还通过减少要处理的数据总量来降低分析的计算复杂性。

对于scRNASeq数据，我们将关注无需监督的特征选择方法，这些方法不需要任何先验信息，例如细胞类型标记或生物组，因为对于许多实验它们不可用，或者可能不可靠。相反，差异表达（第3.6章）可以被认为是监督特征选择的一种形式，因为它使用每个样本的已知生物标记来识别在不同组间表达的不同水平的特征（即基因）。
仍然使用Deng的数据：

```
deng <- readRDS("deng/deng-reads.rds")
cellLabels <- colData(deng)$cell_type2
```
可以使用M3Drop对该数据进行QCed并对文库大小进行标准化，从而去除检测到的基因很少的细胞，去除未检测到的基因，并将原始计数转换为CPM。

```
deng_list <- M3DropCleanData(
    counts(deng),
    labels = cellLabels,
    min_detected_genes = 100,
    is.counts = TRUE
)
expr_matrix <- deng_list$data # Normalized & filtered expression matrix
celltype_labs <- factor(deng_list$labels) # filtered cell-type labels
cell_colors <- brewer.pal(max(3,length(unique(celltype_labs))), "Set3")
```
##### 3.3.1 Identify Genes vs a Null model
无监督特征选择有两种主要方法。第一个是识别与仅描述数据集中预期的技术噪声的空模型不同的基因。

如果数据集包含spike-in RNA，它们可用于直接模拟技术噪声。然而，加spike-in的测量可能不会像内源性转录物那样经历相同的技术噪声（Svensson等，2017）。此外，scRNASeq实验通常只包含少量的加spike-in，这降低了我们对拟合模型参数的信心。
###### 3.3.1.1 Highly Variable Genes
提出用于鉴定scRNASeq数据集中的特征的第一种方法是鉴定高度可变的基因（HVG）。 HVG假设如果基因在细胞中的表达差异很大，那么这些差异中的一些差异是由于细胞之间的生物学差异而不是技术噪声。然而，由于计数数据的性质，基因的平均表达与细胞间读数的变化之间存在正相关关系。必须纠正这种关系才能正确识别HVG。
Brennecke等人提出了一种校正方差和平均表达之间关系的流行方法。为了使用Brennecke方法，我们首先对文库大小进行标准化，然后计算平均值和平方变异系数（变化除以平方均值）表达）。对于ERCC spike-in，二次曲线拟合这两个变量之间的关系，然后使用卡方检验找到显着高于曲线的基因。此方法作为Brennecke_getVariableGenes（count，spike）函数包含在M3Drop包中。但是，此数据集不包含spike-ins，因此我们将使用整个数据集来估算技术噪音。

```
plot(rowMeans(expr_matrix),rowVars(expr_matrix),log='xy')
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/26-dropouts_files/figure-html/unnamed-chunk-6-1.png)
在下图中，红色曲线是拟合的技术噪声模型，虚线是95％CI。粉红点是多次测试校正后具有显着生物变异性的基因。
```
Brennecke_HVG <- BrenneckeGetVariableGenes(
    expr_matrix,
    fdr = 0.01,
    minBiolDisp = 0.5
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/26-dropouts_files/figure-html/unnamed-chunk-7-1.png)
###### 3.3.1.2 High Dropout Genes
寻找HVG的另一种方法是鉴定具有意外大量零的基因。被称为“dropout rate”的零的频率与scRNASeq数据中的表达水平密切相关。零是单细胞RNASeq数据的主要特征，通常占最终表达矩阵中一半以上的条目。这些零点主要是由于mRNA未能被逆转录而失败（Andrews和Hemberg，2016）。逆转录是酶反应，因此可以使用Michaelis-Menten方程建模：
```math
P_dropout=1-S/(K+S)
```
S是细胞中的mRNA浓度（我们将其估计为平均表达）和ķ是Michaelis-Menten常数。
因为Michaelis-Menten方程是一个凸非线性函数，在我们的数据集中，两个或更多的细胞群的不同表达的基因将会被转移到Michaelis-Menten模型中.

```
K <- 49
S_sim <- 10^seq(from = -3, to = 4, by = 0.05) # range of expression values
MM <- 1 - S_sim / (K + S_sim)
plot(
    S_sim, 
    MM, 
    type = "l", 
    lwd = 3, 
    xlab = "Expression", 
    ylab = "Dropout Rate", 
    xlim = c(1,1000)
)
S1 <- 10
P1 <- 1 - S1 / (K + S1) # Expression & dropouts for cells in condition 1
S2 <- 750
P2 <- 1 - S2 / (K + S2) # Expression & dropouts for cells in condition 2
points(
    c(S1, S2),
    c(P1, P2), 
    pch = 16, 
    col = "grey85", 
    cex = 3
)
mix <- 0.5 # proportion of cells in condition 1
points(
    S1 * mix + S2 * (1 - mix), 
    P1 * mix + P2 * (1 - mix), 
    pch = 16, 
    col = "grey35", 
    cex = 3
)

```
我们使用M3Drop来识别在MM曲线右侧的重要异常值。我们也应用了1%的FDR多次测试修正：

```
M3Drop_genes <- M3DropFeatureSelection(
    expr_matrix,
    mt_method = "fdr",
    mt_threshold = 0.01
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/26-dropouts_files/figure-html/unnamed-chunk-10-1.png)
M3Drop包中包含一种替代方法，该方法专门针对UMI标记数据而定制，该数据通常包含由低序列覆盖率产生的许多零，这是由于逆转录不足产生的。 该模型是深度调整的负二项式（Depth-Adjusted Negative Binomial ，DANB）。 该方法将每个表达观察描述为负二项式模型，其具有与相应基因的平均表达和相应细胞的测序深度相关的平均值，以及与基因的平均表达相关的方差。

与Michaelis-Menten和HVG方法不同，对于该模型选择的特征没有可靠的统计检验，因此我们将考虑前1500个基因。

```
deng_int <- NBumiConvertToInteger(counts(deng))
DANB_fit <- NBumiFitModel(deng_int) # DANB is fit to the raw count matrix
# Perform DANB feature selection
DropFS <- NBumiFeatureSelectionCombinedDrop(DANB_fit)
DANB_genes <- names(DropFS[1:1500])
```
##### 3.3.2 Correlated Expression
一种完全不同的特征选择方法是使用基因 - 基因相关性。 该方法基于多个基因将在不同细胞类型或细胞状态之间差异表达的想法。 在相同细胞群中表达的基因将彼此正相关，其中在不同细胞群中表达的基因将彼此负相关。 因此，重要的基因可以通过它们与其他基因的相关程度来鉴定。
这种方法的局限性在于它假设技术噪声对于每个细胞是随机且独立的，因此不应产生基因 - 基因相关性，但是这种假设被批次效应所违反，批次效应通常在不同实验批次之间是系统的并且将产生基因 - 基因相关性。 因此，通过基因 - 基因相关性排序前几千个基因比考虑相关性的显着性更合适。

```
cor_mat <- cor(t(expr_matrix), method = "spearman") # Gene-gene correlations
diag(cor_mat) <- rep(0, times = nrow(expr_matrix))
score <- apply(cor_mat, 1, function(x) {max(abs(x))}) #Correlation of highest magnitude
names(score) <- rownames(expr_matrix);
score <- score[order(-score)]
Cor_genes <- names(score[1:1500])
```
最后，scRNASeq数据中另一种常用的特征选择方法是使用PCA loading。 具有高PCA负荷的基因可能是高度可变的并且与许多其他可变基因相关，因此可能与潜在的生物学相关。 然而，与基因 - 基因相关性一样，PCA载量往往易于检测由于批次效应引起的系统变异; 因此，建议绘制PCA结果以确定与生物变异相对应的那些组分而不是批次效应。

```
# PCA is typically performed on log-transformed expression data
pca <- prcomp(log(expr_matrix + 1) / log(2))

# plot projection
plot(
    pca$rotation[,1], 
    pca$rotation[,2], 
    pch = 16, 
    col = cell_colors[as.factor(celltype_labs)]
) 
```
![image]https://hemberg-lab.github.io/scRNA.seq.course/26-dropouts_files/figure-html/unnamed-chunk-13-1.png)
```
# calculate loadings for components 1 and 2
score <- rowSums(abs(pca$x[,c(1,2)])) 
names(score) <- rownames(expr_matrix)
score <- score[order(-score)]
PCA_genes <- names(score[1:1500])
```
##### 3.3.3 Comparing Methods
我们可以检查这些被识别的特征是否真的代表了在这个数据集中的细胞类型之间的不同表达。

```
MDrop：：M3DropExpressionHeatmap(
    M3Drop_genes,
    expr_matrix,
    cell_labels = celltype_labs
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/26-dropouts_files/figure-html/unnamed-chunk-15-1.png)
我们还可以考虑使用Jaccard索引的查看特性选择方法的一致性：

```
J <- sum(M3Drop_genes %in% HVG_genes)/length(unique(c(M3Drop_genes, HVG_genes)))
```
#### 3.4 Pseudotime analysis

```
library(SingleCellExperiment)
library(TSCAN)
library(M3Drop)
library(monocle)
library(destiny)
library(SLICER)
library(ouija)
library(scater)
library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
library(corrplot)
set.seed(1)
```
在许多情况下，人们正在研究细胞不断变化的过程。这包括，例如，在开发过程中发生的许多分化过程：在刺激之后，细胞将从一种细胞类型变为另一种细胞类型。理想情况下，我们希望监测单个细胞随时间的表达水平。不幸的是，用scRNA-seq进行这种监测是不可能的，因为当提取RNA时细胞被裂解（破坏）。

相反，我们必须在多个时间点进行采样并获得基因表达谱的快照。由于一些细胞将在分化过程中比其他细胞更快地进行，因此每个快照可以在沿着发育进程的不同点处包含细胞。我们使用统计方法沿着一个或多个轨迹对细胞进行排序，这些轨迹代表潜在的发育轨迹，这种排序被称为“pseudotime”。

在本章中，我们将考虑五种不同的工具：Monocle，TSCAN，destiny，SLICER和ouija，用于根据伪时间发展对细胞进行排序。为了说明我们将使用小鼠胚胎发育数据集的方法（Deng et al.2014）。该数据集由来自早期小鼠发育的10个不同时间点的268个细胞组成。在这种情况下，不需要伪时间对准，因为单元标签提供关于发展轨迹的信息。因此，标签允许我们建立基本事实，以便我们可以评估和比较不同的方法。

Cannoodt等人最近的综述提供了来自单细胞转录组学的轨迹推断的各种计算方法的详细总结（Cannoodt，Saelens和Saeys 2016）。他们讨论了几种工具，但不幸的是，出于我们的目的，许多这些工具没有完整或维护良好的实现，和/或没有在R中实现。
下面的图表总结了各种工具的一些特性：
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/cannoodt_pseudotime_properties.png)
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/cannoodt_pseudotime_methods.png)
##### 3.4.1 First look at Deng data
让我们先来看看Deng的数据，还没有应用复杂的pseudotime方法。 如下图所示，简单的PCA可以很好地显示这些数据中的结构。 只有在胚细胞类型（“earlyblast”，“midblat”，“latenlat”）时，PCA难以分离不同的细胞类型。

```
deng_SCE <- readRDS("deng/deng-reads.rds")
deng_SCE$cell_type2 <- factor(
    deng_SCE$cell_type2,
    levels = c("zy", "early2cell", "mid2cell", "late2cell",
                        "4cell", "8cell", "16cell", "earlyblast",
                        "midblast", "lateblast")
)
cellLabels <- deng_SCE$cell_type2
deng <- counts(deng_SCE)
colnames(deng) <- cellLabels
deng_SCE <- plotPCA(deng_SCE, colour_by = "cell_type2", 
                    return_SCE = TRUE) 
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/data-overview-1.png)

```
deng_SCE$PC1 <- reducedDim(deng_SCE, "PCA")[,1]
ggplot(as.data.frame(colData(deng_SCE)), aes(x = PC1, y = cell_type2, 
                              colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("First principal component") + ylab("Timepoint") +
    ggtitle("Cells ordered by first principal component")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/pca-pseudotime-1.png)
正如上面的图所显示的，PC1在发育时期的早期和晚期都很难正确地排列细胞，但是总的来说，在发育的时间点上对细胞进行排序是相对较好的。
##### 3.4.2 TSCAN
TSCAN将聚类与pseudotime分析相结合。 首先，它使用mclust聚类细胞，mclust基于正态分布的混合。 然后，它构建一个最小生成树来连接集群。 连接最大数量的簇的该树的分支是用于确定pseudotime的主分支。

```
procdeng <- TSCAN::preprocess(deng)
colnames(procdeng) <- 1:ncol(deng)
dengclust <- TSCAN::exprmclust(procdeng, clusternum = 10)
TSCAN::plotmclust(dengclust)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/tscan-all-genes-1.png)

```
dengorderTSCAN <- TSCAN::TSCANorder(dengclust, orderonly = FALSE)
pseudotime_order_tscan <- as.character(dengorderTSCAN$sample_name)
deng_SCE$pseudotime_order_tscan <- NA
deng_SCE$pseudotime_order_tscan[as.numeric(dengorderTSCAN$sample_name)] <- 
    dengorderTSCAN$Pseudotime
```
令人沮丧的是，TSCAN只提供268个细胞中的221个的pseudotime，默默地返回未分配cell的缺失值。

```
cellLabels[dengclust$clusterid == 10]
```

```
##  [1] late2cell late2cell late2cell late2cell late2cell late2cell late2cell
##  [8] late2cell late2cell late2cell
## 10 Levels: zy early2cell mid2cell late2cell 4cell 8cell ... lateblast
```

```
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_order_tscan, 
           y = cell_type2, colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("TSCAN pseudotime") + ylab("Timepoint") +
    ggtitle("Cells ordered by TSCAN pseudotime")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/tscan-vs-truth-1.png)
TSCAN的发展轨迹是“错误的方式”，因为后来的pseudotime 值对应于早期时间点，反之亦然。 这本身并不是一个问题（很容易扭转排序以获得伪时间的直观解释），但总体而言，建议TSCAN在此数据集上的性能优于PCA将是一个延伸。 （因为这是一种基于PCA的方法，也许这并不奇怪。）
##### 3.4.3 monocle
生成树以连接所有细胞。 Monocle然后识别该树中最长的路径以确定伪时间。 如果数据包含发散轨迹（即一种细胞类型区分为两种不同的细胞类型），则Monacle可以识别这些。 每个得到的分叉路径被定义为单独的cell状态。

不幸的是，Monocle在使用所有基因时都不起作用，因此我们必须进行特征选择。 首先，我们使用M3Drop：

```
m3dGenes <- as.character(
    M3DropFeatureSelection(deng)$Gene
)
```

```
## Warning in bg__calc_variables(expr_mat): Warning: Removing 1134 undetected
## genes.
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/m3d-select-genes-1.png)

```
d <- deng[which(rownames(deng) %in% m3dGenes), ]
d <- d[!duplicated(rownames(d)), ]
```

```
colnames(d) <- 1:ncol(d)
geneNames <- rownames(d)
rownames(d) <- 1:nrow(d)
pd <- data.frame(timepoint = cellLabels)
pd <- new("AnnotatedDataFrame", data=pd)
fd <- data.frame(gene_short_name = geneNames)
fd <- new("AnnotatedDataFrame", data=fd)

dCellData <- newCellDataSet(d, phenoData = pd, featureData = fd, expressionFamily = tobit())
dCellData <- setOrderingFilter(dCellData, which(geneNames %in% m3dGenes))
dCellData <- estimateSizeFactors(dCellData)
dCellDataSet <- reduceDimension(dCellData, pseudo_expr = 1)
dCellDataSet <- orderCells(dCellDataSet, reverse = FALSE)
plot_cell_trajectory(dCellDataSet)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/monocle-all-genes-1.png)

```
# Store the ordering
pseudotime_monocle <-
    data.frame(
        Timepoint = phenoData(dCellDataSet)$timepoint,
        pseudotime = phenoData(dCellDataSet)$Pseudotime,
        State = phenoData(dCellDataSet)$State
    )
rownames(pseudotime_monocle) <- 1:ncol(d)
pseudotime_order_monocle <-
    rownames(pseudotime_monocle[order(pseudotime_monocle$pseudotime), ])
```

```
deng_SCE$pseudotime_monocle <- pseudotime_monocle$pseudotime
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_monocle, 
           y = cell_type2, colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("monocle pseudotime") + ylab("Timepoint") +
    ggtitle("Cells ordered by monocle pseudotime")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/monocle-vs-truth-1.png)
Monocle - 至少是默认设置 - 对这些数据表现不佳。 “late2cell”组与“zy”，“early2cell”和“mid2cell”细胞完全分离（尽管这些是正确排序的），并且“4cell”，“8cell”，“16cell”或 任何blast胚细胞组之间没有分离
##### 3.4.4 Diffusion maps
扩散图由Ronald Coifman和Stephane Lafon引入，其基本思想是假设数据是来自扩散过程的样本。 该方法通过估计与数据相关的扩散算子的特征值和特征向量来推断低维流形。

Haghverdi等人已将扩散图概念应用于单细胞RNA-seq数据的分析，以创建称为destiny的R包。

我们将在这里将第一扩散映射分量中的cell的ranko prder作为“diffusion map pseudotime”。

```
deng <- logcounts(deng_SCE)
colnames(deng) <- cellLabels
dm <- DiffusionMap(t(deng))

tmp <- data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  Timepoint = deng_SCE$cell_type2)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
    geom_point() + scale_color_tableau() + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic()
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/destiny-deng-1.png)
```
deng_SCE$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_diffusionmap, 
           y = cell_type2, colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("Diffusion map pseudotime (first diffusion map component)") +
    ylab("Timepoint") +
    ggtitle("Cells ordered by diffusion map pseudotime")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/destiny-deng-2.png)
像其他方法一样，使用第一个扩散映射组件作为pseudotime，在排序早期的时间点上做得很好（如果我们在开发中采用“较早”的高值），但是它无法区分后面的时间点。
##### SLICER
SLICER方法是用于构建轨迹的算法，其描述在顺序生物过程期间基因表达变化，就像Monocle和TSCAN一样。 SLICER旨在捕获高度非线性的基因表达变化，自动选择与过程相关的基因，并检测轨迹中的多个分支和环路特征（Welch，Hartemink和Prins 2016）。 SLICER R软件包可从其GitHub存储库获得，可以使用devtools软件包从那里安装。

我们使用SLICER中的select_genes函数自动选择用于构建细胞轨迹的基因。 该函数使用“邻域方差”来识别在整组细胞中平滑变化而不是随机波动的基因。 在此之后，我们确定哪个“k”值（最近邻居的数量）产生最类似于轨迹的嵌入。 然后我们估计细胞的局部线性嵌入。

```
library("lle")
slicer_genes <- select_genes(t(deng))
k <- select_k(t(deng[slicer_genes,]), kmin = 30, kmax=60)
```

```
slicer_traj_lle <- lle(t(deng[slicer_genes,]), m = 2, k)$Y
```

```
reducedDim(deng_SCE, "LLE") <- slicer_traj_lle
plotReducedDim(deng_SCE, use_dimred = "LLE", colour_by = "cell_type2") +
    xlab("LLE component 1") + ylab("LLE component 2") +
    ggtitle("Locally linear embedding of cells from SLICER")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/slicer-analyis-1.png)
通过计算局部线性嵌入，我们可以构造完全连接的k-最近邻图。 该图显示每个单元格的（黄色）圆圈，单元格ID号码以蓝色覆盖。 这里我们显示使用10个最近邻居计算的图表。 在这里，SLICER似乎通过一个分支检测到一个主要轨迹。
```
slicer_traj_graph <- conn_knn_graph(slicer_traj_lle, 10)
plot(slicer_traj_graph, main = "Fully connected kNN graph from SLICER")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/slicer-build-graph-1.png)
从这张图中我们可以识别出“极端”细胞，它们是轨迹中的起始/结束细胞的候选细胞

```
ends <- find_extreme_cells(slicer_traj_graph, slicer_traj_lle)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/slicer-1.png)


```
start <- ends[1]
```
定义了一个开始cell后，我们可以在估计的pseudotime内对cell进行排序。

```
pseudotime_order_slicer <- cell_order(slicer_traj_graph, start)
branches <- assign_branches(slicer_traj_graph, start)

pseudotime_slicer <-
    data.frame(
        Timepoint = cellLabels,
        pseudotime = NA,
        State = branches
    )
pseudotime_slicer$pseudotime[pseudotime_order_slicer] <-
    1:length(pseudotime_order_slicer)
deng_SCE$pseudotime_slicer <- pseudotime_slicer$pseudotime
```

```
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_slicer, 
           y = cell_type2, colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("SLICER pseudotime (cell ordering)") +
    ylab("Timepoint") +
    theme_classic()
```
与前面的方法一样，SLICER为早期的时间点提供了良好的排序。它在“8cell”细胞之前放置“16cell”细胞，但比许多早期的方法提供了更好的blast细胞排序
##### 3.4.6 Ouija
Ouija（http://kieranrcampbell.github.io/ouija/）采用了与我们迄今为止所研究的伪时间估计方法不同的方法。 早期的方法都是“无人监督”的，也就是说，除了可能选择信息基因之外，我们不提供任何关于我们如何期望某些基因或整个轨迹表现的先验信息的方法。

相比之下，Ouija是一种概率框架，它允许仅使用一小组标记基因来解释单细胞pseudotime的学习。 这个方法：
- 从少数标记基因中推断出pseudotime，让您了解为什么pseudotime已经从这些基因中学习出来;
- 为可解释的基因调控行为（如peak time或upregulation time）提供参数估计（具有不确定性）;
- 进行了贝叶斯假设检验，以发现沿着轨迹在其他人之前受到调节的基因;
- 识别亚稳态，即沿连续轨迹的离散细胞类型。
我们将向Ouija提供以下标记基因（在时间点上，它们被认为是高度表达的）：
- Early timepoints: Dazl, Rnf17, Sycp3, Nanog, Pou5f1, Fgf8, Egfr, Bmp5, Bmp15
- Mid timepoints: Zscan4b, Foxa1, Prdm14, Sox21
- Late timepoints: Creb3, Gpx4, Krt8, Elf5, Eomes, Cdx2, Tdgf1, Gdf3
使用Ouija，我们可以将基因模型表现为单调上调或下调（称为开关样行为），或基因短暂达到峰值的瞬态行为。 默认情况下，Ouija假设所有基因都表现出类似开关的行为
在这里，我们可以“cheat”，并检查我们选择的标记基因是否确实识别了分化过程的不同时间点。

```
ouija_markers_down <- c("Dazl", "Rnf17", "Sycp3", "Fgf8", 
                        "Egfr", "Bmp5", "Bmp15", "Pou5f1")
ouija_markers_up <- c("Creb3", "Gpx4", "Krt8", "Elf5", "Cdx2", 
                      "Tdgf1", "Gdf3", "Eomes")
ouija_markers_transient <- c("Zscan4b", "Foxa1", "Prdm14", "Sox21")
ouija_markers <- c(ouija_markers_down, ouija_markers_up, 
                   ouija_markers_transient)
plotExpression(deng_SCE, ouija_markers, x = "cell_type2", colour_by = "cell_type2") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
```
![image](https://note.youdao.com/yws/public/resource/12beaf1b878bc1059015bf8c805c27d5/xmlnote/33E31DDE55D94CE1AEEF02639EABBC86/2019)
为了适应pseudotime，只需调用ouija，传入预期的响应类型。 请注意，如果没有提供响应类型，则默认情况下它们都被假定为switch-like，我们将在此处执行此操作。 Ouija的输入可以是非负表达式值的cell-by-gene矩阵，也可以是ExpressionSet对象，或者，幸运的是，通过从SingleCellExperiment对象中选择logcounts值。
我们可以应用先前的信息，关于基因在分化过程中是否向上或向下调节，并提供关于何时表达的转换或表达的峰值可能发生的先验信息。
- Hamiltonian Monte Carlo（HMC） - 完整的MCMC推理，其中log-posterior的梯度信息用于“引导”随机遍历参数空间，或者
- 自动微分变分贝叶斯（ADVI或简称VI） - 近似推断，其中KL偏差到近似分布被最小化。
一般而言，HMC将为所有参数提供更准确的推理和近似正确的后验方差。 然而，VB比HMC快几个数量级，虽然它可能低估了后验方差，但是Ouija的作者认为，通常情况下，它通常与HMC一样发现后部pseudotime。
为了帮助Ouija模型，我们为它提供了关于上调和下调基因的开关强度的先验信息。 通过将下调基因的开关强度设置为-10，将上调基因的开关强度设置为10，先验强度标准差为0.5，我们告诉模型我们对这些基因在分化过程中的预期行为充满信心。

```
options(mc.cores = parallel::detectCores())
response_type <- c(rep("switch", length(ouija_markers_down) + 
                           length(ouija_markers_up)), 
                   rep("transient", length(ouija_markers_transient)))
switch_strengths <- c(rep(-10, length(ouija_markers_down)),
                      rep(10, length(ouija_markers_up)))
switch_strength_sd <- c(rep(0.5, length(ouija_markers_down)),
                      rep(0.5, length(ouija_markers_up)))
garbage <- capture.output(
    oui_vb <- ouija(deng_SCE[ouija_markers,],
                    single_cell_experiment_assay = "logcounts", 
                    response_type = response_type,
                    switch_strengths = switch_strengths,
                    switch_strength_sd = switch_strength_sd,
                    inference_type = "vb")
)

print(oui_vb)
```

```
## A Ouija fit with 268 cells and 20 marker genes 
## Inference type:  Variational Bayes 
## (Gene behaviour) Switch/transient: 16 / 4
```
我们可以使用plot_expression函数绘制伪照片上的基因表达以及平均函数（S形或高斯瞬态函数）的最大后验（MAP）估计。
```
plot_expression(oui_vb)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/ouija-plot-exprs-1.png)
我们还可以使用plot_switch_times和plot_transient_times函数可视化基因调控行为在轨迹中何时以切换时间或峰值时间（对于开关样或瞬态基因）的形式发生：
```
plot_switch_times(oui_vb)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/ouija-plot-switch-times-1.png)

```
plot_peak_times(oui_vb)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/ouija-plot-switch-times-2.png)
使用一致性矩阵确定亚稳态。
```

cmo <- consistency_matrix(oui_vb)
plot_consistency(oui_vb)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/ouija-consistency-1.png)

```
cell_classifications <- cluster_consistency(cmo)
```

```
map_pst <- map_pseudotime(oui_vb)
ouija_pseudotime <- data.frame(map_pst, cell_classifications)

ggplot(ouija_pseudotime, aes(x = map_pst, y = cell_classifications)) +
  geom_point() +
  xlab("MAP pseudotime") +
  ylab("Cell classification")
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/ouija-pseudotime-1.png)

```
deng_SCE$pseudotime_ouija <- ouija_pseudotime$map_pst
deng_SCE$ouija_cell_class <- ouija_pseudotime$cell_classifications

ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_ouija, 
           y = cell_type2, colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("Ouija pseudotime") +
    ylab("Timepoint") +
    theme_classic()
```

![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/ouija-pseudotime-2.png)
虽然它可以对标记基因的选择和提供的先前信息敏感，但Ouija在细胞的排序方面表现相当不错。 如果选择不同的标记基因或更改先验，结果会如何变化？

Ouija在这里确定了四个亚稳状态，我们可以将其注释为“zygote / 2cell”，“4/8/16 cell”，“blast1”和“blast2”。

```
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = as.factor(ouija_cell_class), 
           y = pseudotime_ouija, colour = cell_type2)) +
    geom_boxplot() + 
    coord_flip() +
    scale_color_tableau() + theme_classic() +
    xlab("Ouija cell classification") +
    ylab("Ouija pseudotime") +
    theme_classic()
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/ouija-states-1.png)

```
gene_regs <- gene_regulation(oui_vb)
head(gene_regs)
```
一个常见的分析是制定基因的调节顺序。 例如，基因A在基因B之前是否被上调？ 基因C在基因D下调前是否达到峰值？ Ouija根据贝叶斯假设检验来回答这些问题，该检验是关于调节时间的差异（切换时间或峰值时间）是否显着不同于0.这是使用gene_regulation函数进行整理。
```
## # A tibble: 6 x 7
## # Groups:   label, gene_A [6]
##   label        gene_A gene_B mean_difference lower_95 upper_95 significant
##   <chr>        <chr>  <chr>            <dbl>    <dbl>    <dbl> <lgl>      
## 1 Bmp15 - Cdx2 Bmp15  Cdx2           -0.0434 -0.0912   0.0110  FALSE      
## 2 Bmp15 - Cre… Bmp15  Creb3           0.278   0.220    0.330   TRUE       
## 3 Bmp15 - Elf5 Bmp15  Elf5           -0.656  -0.688   -0.618   TRUE       
## 4 Bmp15 - Eom… Bmp15  Eomes           0.0766  0.00433  0.153   TRUE       
## 5 Bmp15 - Fox… Bmp15  Foxa1          -0.0266 -0.0611   0.00287 FALSE      
## 6 Bmp15 - Gdf3 Bmp15  Gdf3            0.0704  0.00963  0.134   TRUE
```
##### 3.4.7 Comparison of the methods
TSCAN，Monocle，Diffusion Map，SLICER和Ouija推断出的轨迹哪一个更好？

TSCAN和Diffusion Map方法使得轨迹“错误”，因此我们将对这些比较进行调整。

```
df_pseudotime <- as.data.frame(
    colData(deng_SCE)[, grep("pseudotime", colnames(colData(deng_SCE)))]
)
colnames(df_pseudotime) <- gsub("pseudotime_", "", 
                                colnames(df_pseudotime))
df_pseudotime$PC1 <- deng_SCE$PC1
df_pseudotime$order_tscan <- -df_pseudotime$order_tscan
df_pseudotime$diffusionmap <- -df_pseudotime$diffusionmap

corrplot.mixed(cor(df_pseudotime, use = "na.or.complete"), 
               order = "hclust", tl.col = "black",
               main = "Correlation matrix for pseudotime results",
               mar = c(0, 0, 3.1, 0))
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/compare-results-1.png)
##### 3.4.8 Expression of genes through time
每个包还可以通过pseudotime显示表达式。 遵循单个基因非常有助于鉴定在分化过程中起重要作用的基因。 我们用Rhoa基因说明了这个过程。

我们已将使用所有方法计算的pseudotime添加到SCE对象的colData中。 完成后，scater包的完整绘图功能可用于研究基因表达，细胞群和假型的关系。 这对于不提供绘图功能的SLICER等软件包特别有用。
Principal components
```
plotExpression(deng_SCE, "Rhoa", x = "PC1", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)
```

![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/Rhoa-pc1-1.png)
TSCAN

```
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_order_tscan", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/Rhoa-tscan-1.png)
Monocle

```
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_monocle", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/Rhoa-monocle-1.png)
Diffusion Map

```
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_diffusionmap", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/Rhoa-diff-map-1.png)
SLICER

```
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_slicer", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)
```

![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/Rhoa-slicer-1.png)
Ouija

```
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_ouija", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/27-pseudotime_files/figure-html/Rhoa-ouija-1.png)
##### 3.5 Imputation
如前所述，分析scRNA-seq数据时的主要挑战之一是存在零或丢失。假设丢失出现了三个可能的原因：

- 该基因不在细胞中表达，因此没有序列的转录物
- 该基因被表达，但由于某种原因，转录物在测序前丢失了某处
- 表达基因并捕获转录物并转化为cDNA，但测序深度不足以产生任何读数。
因此，丢失可能是实验缺陷的结果，如果是这种情况，那么我们想提供计算校正。一种可能的解决方案是将丢失量计算在表达式矩阵中。为了能够克服基因表达值，必须有一个基础模型。然而，由于我们不知道哪些丢失事件是技术人为因素，哪些本来就存在的，因此归咎是一项艰巨的挑战。

据我们所知，目前有两种不同的插补方法：MAGIC（Dijk等人2017）和scImpute（W. V. Li和Li 2017）。 MAGIC仅适用于Python或Matlab，但我们将从R内部运行它。
##### 3.5.1 scImpute
为了测试scImpute，我们使用默认参数，并将其应用到我们以前处理过的Deng数据集。scImpute将一个.csv或.txt文件作为输入：

```
deng <- readRDS("deng/deng-reads.rds")
write.csv(counts(deng), "deng.csv")
scimpute(
    count_path = "deng.csv",
    infile = "csv",
    outfile = "txt", 
    out_dir = "./",
    Kcluster = 10,
    ncores = 2
)
```
现在我们可以通过考虑一个PCA图来比较结果和原始数据

```
res <- read.table("scimpute_count.txt")
colnames(res) <- NULL
res <- SingleCellExperiment(
    assays = list(logcounts = log2(as.matrix(res) + 1)), 
    colData = colData(deng)
)
rowData(res)$feature_symbol <- rowData(deng)$feature_symbol
plotPCA(
    res, 
    colour_by = "cell_type2"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/28-imputation_files/figure-html/unnamed-chunk-4-1.png)
为了评估imputation的影响，我们使用SC3来对估算矩阵进行分组

```
res <- sc3_estimate_k(res)
metadata(res)$sc3$k_estimation
```

```
## [1] 6
```

```
res <- sc3(res, ks = 10, gene_filter = FALSE)
```

```
adjustedRandIndex(colData(deng)$cell_type2, colData(res)$sc3_10_clusters)
```

```
## [1] 0.4687332
```

```
plotPCA(
    res, 
    colour_by = "sc3_10_clusters"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/28-imputation_files/figure-html/unnamed-chunk-5-1.png)
##### 3.5.2 MAGIC

```
system("python3 utils/MAGIC.py -d deng.csv -o MAGIC_count.csv --cell-axis columns -l 1 --pca-non-random csv")
```

```
res <- read.csv("MAGIC_count.csv", header = TRUE)
rownames(res) <- res[,1]
res <- res[,-1]
res <- t(res)
res <- SingleCellExperiment(
    assays = list(logcounts = res), 
    colData = colData(deng)
)
rowData(res)$feature_symbol <- rownames(res)
plotPCA(
    res, 
    colour_by = "cell_type2"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/28-imputation_files/figure-html/unnamed-chunk-7-1.png)
为了评估imputation的影响，我们使用SC3来对估算矩阵进行分组

```
res <- sc3_estimate_k(res)
metadata(res)$sc3$k_estimation
```

```
## [1] 4
```

```
res <- sc3(res, ks = 10, gene_filter = FALSE)
```

```
adjustedRandIndex(colData(deng)$cell_type2, colData(res)$sc3_10_clusters)
```

```
## [1] 0.3752866
```

```
plotPCA(
    res, 
    colour_by = "sc3_10_clusters"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/28-imputation_files/figure-html/unnamed-chunk-8-1.png)
#### 3.6 Differential Expression(DE) analysis
##### 3.6.1 Bulk RNA-seq
在处理大量RNA-seq数据时，最常见的分析类型之一是识别不同表达的基因。通过比较两个条件之间的变化，例如突变和野生型，或者刺激和不受刺激，有可能描述这种变化背后的分子机制。
几种不同的方法，例如 已经开发出DESeq2和edgeR用于大量RNA-seq。 此外，还有广泛的数据集，其中RNA-seq数据已经使用RT-qPCR验证。 这些数据可用于对DE查找算法进行基准测试，现有证据表明算法表现良好。
##### 3.6.2 Single cell  RNA-seq
与大量RNA-seq相反，在scRNA-seq中，我们通常没有一组确定的实验条件。 相反，如前一章（3.2）所示，我们可以通过使用无监督的聚类方法来识别细胞群。 一旦鉴定了组，就可以通过比较各组之间的方差差异（如SC3中实施的Kruskal-Wallis检验），或通过以成对方式比较簇之间的基因表达，找到差异表达的基因。 在下一章中，我们将主要考虑为成对比较开发的工具。
与大量RNA-seq不同，我们通常在单细胞实验中比较的每组中有大量样品（即细胞）。因此，我们可以利用每组中表达值的完整分布来识别组之间的差异，而不是如bulk RNA-seq一样仅仅比较平均表达的估计。

比较分布有两种主要方法。首先，我们可以使用现有的统计模型/分布，并将相同类型的模型拟合到每个组中的表达式，然后测试每个模型的参数差异，或者测试如果允许特定参数不同，模型是否更合适据小组说。例如，在第2.10章中，我们使用edgeR来测试在不同批次中允许平均表达是否不同显着改善了数据的负二项模型的拟合。

或者，我们可以使用非参数测试，该测试不假设表达值遵循任何特定分布，例如， Kolmogorov-Smirnov检验（KS检验）。非参数测试通常将观察到的表达值转换为等级并测试一组的等级分布是否明显不同于另一组的等级分布。然而，一些非参数方法在存在大量绑定值时失败，例如单细胞RNA-seq表达数据中的丢失（零）的情况。此外，如果参数测试的条件成立，那么它通常比非参数测试更强大。
下面用一个测试数据来评价一下不同的算法的表现。处理同样的表达矩阵得到差异结果跟已知的差异结果进行比较看看overlap怎么样。评价指标主要是：
- 召回率，True positive rate (TPR), TP/(TP + FN)
- 准确率，False positive rate (FPR), FP/(FP+TP)
- receiver-operating-characteristic curve (ROC)
- area under this curve (AUC)
3.6.4 Models of single-cell RNAseq Data
RNAseq数据最常见的模型是负二项模型

```
set.seed(1)
hist(
    rnbinom(
        1000, 
        mu = 10, 
        size = 100), 
    col = "grey50", 
    xlab = "Read Counts", 
    main = "Negative Binomial"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/29-de-intro_files/figure-html/nb-plot-1.png)
Negative Binomial distribution of read counts for a single gene across 1000 cells

```math
Mean：μ=mu
```


```math
Variance:σ^2=mu+mu^2/size
```
它通过平均表达式（mu）和分散（size）参数化，其与方差成反比。负二项式模型很好地拟合了大量RNA-seq数据，并且它被用于为这些数据设计的大多数统计方法。此外，它已被证明适合从独特的分子标识符（UMI）标记的数据获得的分子计数的分布相当好.
然而，由于相对于非零读取计数的高丢失率，原始负二项式模型也不适合全长转录本数据。对于这种类型的数据，已经提出了各种zero-inflated负二项模型
```
d <- 0.5;
counts <- rnbinom(
    1000, 
    mu = 10, 
    size = 100
)
counts[runif(1000) < d] <- 0
hist(
    counts, 
    col = "grey50", 
    xlab = "Read Counts", 
    main = "Zero-inflated NB"
)
```

![image](https://hemberg-lab.github.io/scRNA.seq.course/29-de-intro_files/figure-html/zero-inflation-plot-1.png)
 Zero-inflated Negative Binomial distribution
 ```math
Mean：μ=mu*(1-d)
```


```math
Variance:σ^2=mu(1-d)+(1+d*μ+μ/size)
```
由于丢失率，这个模型在负二项模型中引入了一个新的参数d，正如我们在19章中所见到的，一个基因的丢失率与这个基因的平均表达值密切相关。不同的zero-inflated负二项模型利用mu和d之间的不同关系以及一些可能适合的μ和d独立的表达每一个基因。
最后，几种方法使用Poisson-Beta分布，该分部基于转录破裂的机制模型。这个模型有很强的实验支持，但是他比负二项模型更加不容易使用，

```
a <- 0.1
b <- 0.1
g <- 100
lambdas <- rbeta(1000, a, b)
counts <- sapply(g*lambdas, function(l) {rpois(1, lambda = l)})
hist(
    counts, 
    col = "grey50", 
    xlab = "Read Counts", 
    main = "Poisson-Beta"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/29-de-intro_files/figure-html/pois-beta-plot-1.png)
```math
Mean：μ=g*a/(a+b)
```


```math
Variance:σ^2=g^2*a*b/((a+b+1)*(a+b)^2)

```
这个模型使用3个参数：a表示转率激活率；b表示转率抑制率，g转录产生的速率，即转录在该位点是活跃的。差异表达方法可以测试每个参数的组间差异或仅一个（通常是g）。
所有这些模型可以进一步扩展以明确地解释基因表达差异的其他来源，例如批次效应或库深度，这取决于特定的DE算法。
#### 3.7 DE in a real dataset

```
library(scRNA.seq.funcs)
library(edgeR)
library(monocle)
library(MAST)
library(ROCR)
set.seed(1)
```
###### 3.7.1 Introdution
为了测试不同单细胞差异表达方法我们将使用7-17章中的Blischak数据集。在这个实验中除了产生单细胞数据，每一个细胞系还产生了bulk RNA-seq数据。我们将使用在各自的批量数据上使用标准方法确定的不同表达的基因作为评估每个单细胞方法的准确性的基本事实。为了节省时间，我们预先计算了这些。您可以运行下面的命令来加载这些数据。

```
DE <- read.table("tung/TPs.txt")
notDE <- read.table("tung/TNs.txt")
GroundTruth <- list(
    DE = as.character(unlist(DE)), 
    notDE = as.character(unlist(notDE))
)
```

这对NA19101和NA19239的比较的结果。现在加载各自的单细胞数据：

```
molecules <- read.table("tung/molecules.txt", sep = "\t")
anno <- read.table("tung/annotation.txt", sep = "\t", header = TRUE)
keep <- anno[,1] == "NA19101" | anno[,1] == "NA19239"
data <- molecules[,keep]
group <- anno[keep,1]
batch <- anno[keep,4]
# remove genes that aren't expressed in at least 6 cells
gkeep <- rowSums(data > 0) > 5;
counts <- data[gkeep,]
# Library size normalization
lib_size = colSums(counts)
norm <- t(t(counts)/lib_size * median(lib_size)) 
# Variant of CPM for datasets with library sizes of fewer than 1 mil molecules
```
##### 3.7.2 Kolmogorov-Smirnov test
最容易使用的测试类型是非参数化的。最常用的非参数化测试是Kolmogorov-Smirnov test（KS-test），我们可以用它来对比两个个体中每一个基因的分布。
KS-test量化了这两组中每个基因表达的经验累积分布之间的距离。它对平均表达的变化和方差的变化很敏感。然而，它假设数据是连续的，并且当数据包含大量相同的值时可能表现不佳（例如。0)。KS-test的另一个问题是，它对大的样本容量非常敏感，因此它可能会变得非常重要，即使差异的大小非常小。
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/KS2_Example.png)
Illustration of the two-sample Kolmogorov–Smirnov statistic. Red and blue lines each correspond to an empirical distribution function, and the black arrow is the two-sample KS statistic.

```
pVals <- apply(
    norm, 1, function(x) {
        ks.test(
            x[group == "NA19101"], 
            x[group == "NA19239"]
        )$p.value
    }
)
# multiple testing correction
pVals <- p.adjust(pVals, method = "fdr")
```
“applies”函数应用于表达矩阵的每一行（1）。在这个函数中我们仅仅返回；额ks。test输出中的p.value.我们现在可以考虑，KS-test检测到多少个正和负的基因：
###### 3.7.2.1 Evalusting Accuracy

```
sigDE <- names(pVals)[pVals < 0.05]
length(sigDE) 
## [1] 5095
# Number of KS-DE genes
sum(GroundTruth$DE %in% sigDE) 
## [1] 792
# Number of KS-DE genes that are true DE genes
sum(GroundTruth$notDE %in% sigDE)
## [1] 3190
# Number of KS-DE genes that are truly not-DE
```
正如您所看到的，通过KS检验（假阳性）确定了更多假阳性基因而不是真阳性基因（真阳性），但这可能是由于大量的NOTDE基因，因此我们通常 将这些计数归一化为真阳性率（TPR），TP /（TP + FN）和假阳性率（FPR），FP /（FP + TP）。

```
tp <- sum(GroundTruth$DE %in% sigDE)
fp <- sum(GroundTruth$notDE %in% sigDE)
tn <- sum(GroundTruth$notDE %in% names(pVals)[pVals >= 0.05])
fn <- sum(GroundTruth$DE %in% names(pVals)[pVals >= 0.05])
tpr <- tp/(tp + fn)
fpr <- fp/(fp + tn)
cat(c(tpr, fpr))
## 0.7346939 0.2944706
```

现在我们可以看到，TPR比FPR要高得多，表明KS测试可以识别DE基因。
到目前为止，我们只评估了单个显着性阈值的性能。 通常情况下，改变阈值并评估一系列值的性能是有益的。 然后将其绘制为receiver-operating-characteristic curve（ROC），并且可以将一般精度统计量计算为该曲线下的面积（AUC）。 我们将使用ROCR包来促进这种绘图。
可以看到KS检验判断的显著差异基因实在是太多了，高达5095个。所以它能找回来792个真正的差异基因。但是却找到了3190个假阳性。所以计算得到召回率73.46%，但是准确率只有29.44%，这个表现不佳。

```
# Only consider genes for which we know the ground truth
pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
               names(pVals) %in% GroundTruth$notDE] 
truth <- rep(1, times = length(pVals));
truth[names(pVals) %in% GroundTruth$DE] = 0;
pred <- ROCR::prediction(pVals, truth)
perf <- ROCR::performance(pred, "tpr", "fpr")
ROCR::plot(perf)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/30-de-real_files/figure-html/ks-roc-plot-1.png)

```
aucObj <- ROCR::performance(pred, "auc")
aucObj@y.values[[1]] # AUC
```
Finally to facilitate the comparisons of other DE methods let’s put this code into a function so we don’t need to repeat it:

```
DE_Quality_AUC <- function(pVals) {
    pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
                   names(pVals) %in% GroundTruth$notDE]
    truth <- rep(1, times = length(pVals));
    truth[names(pVals) %in% GroundTruth$DE] = 0;
    pred <- ROCR::prediction(pVals, truth)
    perf <- ROCR::performance(pred, "tpr", "fpr")
    ROCR::plot(perf)
    aucObj <- ROCR::performance(pred, "auc")
    return(aucObj@y.values[[1]])
}
```
##### 3.7.3 Wilcox/Mann-Whitney-U Test
Wilcox-rank-sum test是另一种非参数检验，但是这个test对于一组中的值是否小于、大于另一组中的值更加有针对性。因此，它通常被认为是两组之间中值差异的检验，而KS-test对于表达值分布的任何变化都很敏感。

```
pVals <- apply(
    norm, 1, function(x) {
        wilcox.test(
            x[group == "NA19101"], 
            x[group == "NA19239"]
        )$p.value
    }
)
# multiple testing correction
pVals <- p.adjust(pVals, method = "fdr")
DE_Quality_AUC(pVals)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/30-de-real_files/figure-html/wilcox-plot-1.png)

```
## [1] 0.8320326
```
##### 3.7.4 edgeR
我们已经在第2.10章中使用了edgeR来表示差异表达式。 edgeR基于基因表达的负二项式模型，并使用广义线性模型（GLM）框架，使我们能够将其他因素（例如批次）包含到模型中。

```
dge <- DGEList(
    counts = counts, 
    norm.factors = rep(1, length(counts[1,])), 
    group = group
)
group_edgeR <- factor(group)
design <- model.matrix(~ group_edgeR)
dge <- estimateDisp(dge, design = design, trend.method = "none")
fit <- glmFit(dge, design)
res <- glmLRT(fit)
pVals <- res$table[,4]
names(pVals) <- rownames(res$table)

pVals <- p.adjust(pVals, method = "fdr")
DE_Quality_AUC(pVals)

```
![image](https://hemberg-lab.github.io/scRNA.seq.course/30-de-real_files/figure-html/edger-plot-1.png)

```
## [1] 0.8466764
```
##### 3.7.5 Monocle
Monocle可以使用几种不同的DE模型。 对于count data，它建议使用负二项模型（negbinomial.size）。 对于标准化数据，它建议使用正态分布（gaussianff）对其进行对数变换。 与edgeR类似，此方法使用GLM框架，因此理论上可以考虑批次，但实际上，如果包含批次，则模型将无法适用于此数据集。

```
pd <- data.frame(group = group, batch = batch)
rownames(pd) <- colnames(counts)
pd <- new("AnnotatedDataFrame", data = pd)

Obj <- newCellDataSet(
    as.matrix(counts), 
    phenoData = pd, 
    expressionFamily = negbinomial.size()
)
Obj <- estimateSizeFactors(Obj)
Obj <- estimateDispersions(Obj)
res <- differentialGeneTest(Obj, fullModelFormulaStr = "~group")

pVals <- res[,3]
names(pVals) <- rownames(res)
pVals <- p.adjust(pVals, method = "fdr")
DE_Quality_AUC(pVals)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/30-de-real_files/figure-html/Monocle-plot-1.png)
##### 3.7.6 MAST
MAST是基于zero-inflated的负二项模型。它通过一个障碍模型测试差异表达式，将离散（0 vs非零）和连续（非零值）的基因表达的测试结合起来。同样，它使用一个线性建模框架来考虑复杂的模型。

```
log_counts <- log(counts + 1) / log(2)
fData <- data.frame(names = rownames(log_counts))
rownames(fData) <- rownames(log_counts);
cData <- data.frame(cond = group)
rownames(cData) <- colnames(log_counts)

obj <- FromMatrix(as.matrix(log_counts), cData, fData)
colData(obj)$cngeneson <- scale(colSums(assay(obj) > 0))
cond <- factor(colData(obj)$cond)

# Model expression as function of condition & number of detected genes
zlmCond <- zlm.SingleCellAssay(~ cond + cngeneson, obj) 
## Warning: 'zlm.SingleCellAssay' is deprecated.
## Use 'zlm' instead.
## See help("Deprecated")
## Warning in .nextMethod(object = object, value = value): Coefficients
## condNA19239 are never estimible and will be dropped.
summaryCond <- summary(zlmCond, doLRT = "condNA19101")
summaryDt <- summaryCond$datatable

summaryDt <- as.data.frame(summaryDt)
pVals <- unlist(summaryDt[summaryDt$component == "H",4]) # H = hurdle model
names(pVals) <- unlist(summaryDt[summaryDt$component == "H",1])
pVals <- p.adjust(pVals, method = "fdr")
DE_Quality_AUC(pVals)

```
![image](https://hemberg-lab.github.io/scRNA.seq.course/30-de-real_files/figure-html/MAST-plot-1.png)

```
## [1] 0.8284046
```
##### 3.7.7 Slow Methods
这个方法运行时间很长
#####  3.7.8 BPSC
BPSC使用了单细胞基因表达的Poisson-beta模型，我们在前一章中讨论过，并将它与我们在使用edgeR时已经遇到的广义线性模型结合在一起。BPSC将一个或多个组与一个参照组（“control”）进行比较，并且可以包括其他因素，比如模型中的批次。

```
library(BPSC)
bpsc_data <- norm[,batch=="NA19101.r1" | batch=="NA19239.r1"]
bpsc_group = group[batch=="NA19101.r1" | batch=="NA19239.r1"]

control_cells <- which(bpsc_group == "NA19101")
design <- model.matrix(~bpsc_group)
coef=2 # group label
res=BPglm(data=bpsc_data, controlIds=control_cells, design=design, coef=coef, 
                estIntPar=FALSE, useParallel = FALSE)
pVals = res$PVAL
pVals <- p.adjust(pVals, method = "fdr")
DE_Quality_AUC(pVals)
```
##### 3.7.9 SCDE
SCDE是第一个单细胞特异的DE方法。它适用于使用贝叶斯统计数据的零膨胀负二项模型。下面的用法测试了不同群体中单个基因的平均表达差异，但最近的版本包括了测试平均表达或基因组分散的方法，通常代表一条途径。

```
library(scde)
cnts <- apply(
    counts,
    2,
    function(x) {
        storage.mode(x) <- 'integer'
        return(x)
    }
)
names(group) <- 1:length(group)
colnames(cnts) <- 1:length(group)
o.ifm <- scde::scde.error.models(
    counts = cnts,
    groups = group,
    n.cores = 1,
    threshold.segmentation = TRUE,
    save.crossfit.plots = FALSE,
    save.model.plots = FALSE,
    verbose = 0,
    min.size.entries = 2
)
priors <- scde::scde.expression.prior(
    models = o.ifm,
    counts = cnts,
    length.out = 400,
    show.plot = FALSE
)
resSCDE <- scde::scde.expression.difference(
    o.ifm,
    cnts,
    priors,
    groups = group,
    n.randomizations = 100,
    n.cores = 1,
    verbose = 0
)
# Convert Z-scores into 2-tailed p-values
pVals <- pnorm(abs(resSCDE$cZ), lower.tail = FALSE) * 2
DE_Quality_AUC(pVals)
```
#### 3.8 Comparing/Combining scRNAseq datasets
随着越来越多的scRNA-seq数据集变得可用，将merged-seurat进行比较是关键。比较scRNASeq数据集有两种主要方法。第一种方法是“label-centric”，它专注于通过比较单个细胞或细胞群来识别不同数据集的等效细胞类型/状态。另一种方法是“交叉数据标准化(cross-dataset normalization)”，它试图通过计算来消除特定于实验的技术/生物效应，从而使来自多个实验的数据能够被组合和联合分析。
以标签为中心的方法可以与具有高置信度细胞注释的数据集一起使用，例如，人类细胞图谱（HCA）（Regev等人2017）或Tabula Muris（???）一旦完成，将细胞或簇从新样本投射到该参考上以考虑组织成分和/或识别细胞小说/未知身份。从概念上讲，这种预测类似于流行的BLAST方法（Altschul等人，1990），这使得可以在数据库中快速找到最接近的匹配，用于新鉴定的核苷酸或氨基酸序列。以标签为中心的方法还可用于比较由不同实验室收集的相似生物来源的数据集，以确保注释和分析是一致的
cross-dataset标准化方法也可以用来比较类似的生物起源的数据集，不像label-centric，它支持多个数据集的连接分析，以方便识别罕见的细胞类型，这些类型可能在每个单独的数据集中被很稀疏地取样，以便可靠地检测到。然而，cross-dataset规范化并不适用于非常大的和多样化的引用，因为它假定每个数据集的大部分生物变量都与其他数据集重叠。
###### 3.8.2 Datasets
我们将在两个人类胰腺数据集上运行这些方法。由于胰腺已经被广泛研究，这些数据集都有很好的注释。

```
muraro <- readRDS("pancreas/muraro.rds")
segerstolpe <- readRDS("pancreas/segerstolpe.rds")
```

此数据已针对scmap格式化。 cell类型标签必须存储在colData的cell_type1列中，并且两个数据集之间一致的基因ID必须存储在rowData的feature_symbol列中。
首先，让我们检查两个数据集的基因id匹配：

```
sum(rowData(muraro)$feature_symbol %in% rowData(segerstolpe)$feature_symbol)/nrow(muraro)
```

```
## [1] 0.9599519
```

```
sum(rowData(segerstolpe)$feature_symbol %in% rowData(muraro)$feature_symbol)/nrow(segerstolpe)
```

```
## [1] 0.719334
```
在这里，我们可以看到muraro中存在的96％基因与segerstople中的基因匹配，而segerstolpe中72％的基因是muraro中的匹配基因。 这是预期的，因为segerstolpe数据集比muraro数据集更深入地排序。 然而，它强调了比较scRNASeq数据集的一些困难。
我们可以通过检查这两个数据集的总体大小来确认这一点。

```
dim(muraro)
## [1] 19127  2126
dim(segerstolpe)
## [1] 25525  3514
```
此外，我们还可以使用下面的命令检查每个数据集的cell类型注解：

```
summary(factor(colData(muraro)$cell_type1))
##      acinar       alpha        beta       delta      ductal endothelial 
##         219         812         448         193         245          21 
##     epsilon       gamma mesenchymal     unclear 
##           3         101          80           4
summary(factor(colData(segerstolpe)$cell_type1))
##                 acinar                  alpha                   beta 
##                    185                    886                    270 
##          co-expression                  delta                 ductal 
##                     39                    114                    386 
##            endothelial                epsilon                  gamma 
##                     16                      7                    197 
##                   mast           MHC class II         not applicable 
##                      7                      5                   1305 
##                    PSC           unclassified unclassified endocrine 
##                     54                      2                     41
```
在这里我们可以看到，尽管两个数据集都将两个数据集视为相同的生物组织，但它们已经注释了一组略有不同的细胞类型。 如果您熟悉胰腺生物学，您可能会认识到segerstolpe中的胰腺星状细胞（PSCs）是一种间充质干细胞，在muraro中属于“间充质”类型。 但是，目前尚不清楚这两个注释是否应该被视为同义词。 我们可以使用以标签为中心的比较方法来确定这两个细胞类型的注释是否确实相同。
或者，我们可能有兴趣了解那些“未分类的内分泌”细胞的功能，或者通过利用数据集中的形成来理解每个数据集中原始聚类的质量太差（“不适用”）。 我们可以尝试使用以标签为中心的方法来推断它们最可能属于哪些现有注释，或者我们可以尝试使用跨数据集规范化来揭示其中的新颖cell类型（或现有注释中的子类型）。
为了简化我们的演示分析，我们将删除少量未分配的cell和质量较差的cell。 我们将保留“未分类的内分泌”，看看这些方法中是否有任何一种可以阐明它们属于哪种细胞类型。

```
segerstolpe <- segerstolpe[,colData(segerstolpe)$cell_type1 != "unclassified"]
segerstolpe <- segerstolpe[,colData(segerstolpe)$cell_type1 != "not applicable",]
muraro <- muraro[,colData(muraro)$cell_type1 != "unclear"]
```
##### 3.8.3 Projecting cells onto annotated cell-types（scmap）
我们最近开发了scmap（kis列夫和Hemberg，2017）——一种将细胞从scRNA-seq实验投射到其他实验中确定的细胞类型的方法。
###### 3.8.3.1 Feature Selection
一旦我们有了SingleCellExperiment对象，我们就可以运行scmap了。 首先，我们必须构建参考集群的“索引”。 由于我们想知道PSC和间充质细胞是否是同义词，我们将每个数据集投射到另一个数据集，因此我们将为每个数据集构建索引。 这需要首先为参考数据集选择信息量最大的特征。

```
muraro <- selectFeatures(muraro, suppress_plot = FALSE)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/31-projection_files/figure-html/unnamed-chunk-12-1.png)

```
segerstolpe <- selectFeatures(segerstolpe, suppress_plot = FALSE)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/31-projection_files/figure-html/unnamed-chunk-13-1.png)
从这些图的y轴我们可以看到scmap使用基于dropmerged-seserat的特征选择方法。
计算cell-type索引：

```
muraro <- indexCluster(muraro)
segerstolpe <- indexCluster(segerstolpe)
```
可视化索引：

```
heatmap(as.matrix(metadata(muraro)$scmap_cluster_index))
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/31-projection_files/figure-html/unnamed-chunk-15-1.png)
如果功能过于集中在少数单元类型中，您可能需要使用setFeatures函数调整功能。 在这种情况下，基于dropmerged_seurat的功能看起来很好
###### 3.8.3.2 Projecting
scmap计算从每个cell到参考索引中每个cell类型的距离，然后应用经验导出的阈值来确定将哪些cell分配给最近的参考cell类型以及哪些cell不分配。为了解释测序深度的差异，使用spearman相关和余弦距离计算距离，并且仅将具有两个距离的一致分配的单元返回为指定的。
我们将把segerstolpe数据集投射到muraro数据集：

```
seger_to_muraro <- scmapCluster(
  projection = segerstolpe,
  index_list = list(
    muraro = metadata(muraro)$scmap_cluster_index
  )
)
```
muraro 投射到 segerstolpe

```
muraro_to_seger <- scmapCluster(
  projection = muraro,
  index_list = list(
    seger = metadata(segerstolpe)$scmap_cluster_index
  )
)
```
现在，让我们将原始的细胞类型标签与投影标签进行比较：

```
table(colData(muraro)$cell_type1, muraro_to_seger$scmap_cluster_labs)
```

```
##              
##               acinar alpha beta co-expression delta ductal endothelial
##   acinar         211     0    0             0     0      0           0
##   alpha            1   763    0            18     0      2           0
##   beta             2     1  397             7     2      2           0
##   delta            0     0    2             1   173      0           0
##   ductal           7     0    0             0     0    208           0
##   endothelial      0     0    0             0     0      0          15
##   epsilon          0     0    0             0     0      0           0
##   gamma            2     0    0             0     0      0           0
##   mesenchymal      0     0    0             0     0      1           0
##              
##               epsilon gamma MHC class II PSC unassigned
##   acinar            0     0            0   0          8
##   alpha             0     2            0   0         26
##   beta              0     5            1   2         29
##   delta             0     0            0   0         17
##   ductal            0     0            5   3         22
##   endothelial       0     0            0   1          5
##   epsilon           3     0            0   0          0
##   gamma             0    95            0   0          4
##   mesenchymal       0     0            0  77          2
```
在这里，我们可以看到细胞类型在segerstolpe中映射到它们的等价物，而且重要的是我们看到所有的“mesenchymal”细胞都被分配到“PSC”类中。

```
table(colData(segerstolpe)$cell_type1, seger_to_muraro$scmap_cluster_labs)
```

```
##                         
##                          acinar alpha beta delta ductal endothelial
##   acinar                    181     0    0     0      4           0
##   alpha                       0   869    1     0      0           0
##   beta                        0     0  260     0      0           0
##   co-expression               0     7   31     0      0           0
##   delta                       0     0    1   111      0           0
##   ductal                      0     0    0     0    383           0
##   endothelial                 0     0    0     0      0          14
##   epsilon                     0     0    0     0      0           0
##   gamma                       0     2    0     0      0           0
##   mast                        0     0    0     0      0           0
##   MHC class II                0     0    0     0      0           0
##   PSC                         0     0    1     0      0           0
##   unclassified endocrine      0     0    0     0      0           0
##                         
##                          epsilon gamma mesenchymal unassigned
##   acinar                       0     0           0          0
##   alpha                        0     0           0         16
##   beta                         0     0           0         10
##   co-expression                0     0           0          1
##   delta                        0     0           0          2
##   ductal                       0     0           0          3
##   endothelial                  0     0           0          2
##   epsilon                      6     0           0          1
##   gamma                        0   192           0          3
##   mast                         0     0           0          7
##   MHC class II                 0     0           0          5
##   PSC                          0     0          53          0
##   unclassified endocrine       0     0           0         41
```
我们再次看到细胞类型相互匹配，其中“PSCs”匹配“mesenchymal”细胞，提供强有力的证据证明这两个注释应该被认为是同义的。
使用Sankey diagram进行可视化：
```
plot(getSankey(colData(muraro)$cell_type1,  muraro_to_seger$scmap_cluster_labs[,1], plot_height=400))
```
##### 3.8.4 Cell-to-Cell mapping
scmap还可以将一个数据集中的每个cell投射到参考数据集中的近似最近的相邻cell中。它使用了一个高度优化的搜索算法，允许将其扩展到非常大的引用（在理论上是10万-数百万个cell）。然而，这个过程是随机的，所以我们必须修正随机的种子，以确保我们能够重现我们的结果。
我们已经为这个数据集执行了特性选择，这样我们就可以直接构建索引了。

```
set.seed(193047)
segerstolpe <- indexCell(segerstolpe)
muraro <- indexCell(muraro)
```
在这种情况下，索引是每个cell的一系列集群，使用不同的特性集，参数k和M是集群的数量和每个子集群中使用的特性的数量。在每个子集群中，将新cell分配给最近的集群，以生成独特的集群分配模式。然后我们在参考资料集中找到具有相同或最相似的集群分配模式的cell。
我们可以使用以下方法检查参考数据集的集群分配模式：

```
metadata(muraro)$scmap_cell_index$subclusters[1:5,1:5]
```


```
##      D28.1_1 D28.1_13 D28.1_15 D28.1_17 D28.1_2
## [1,]       4       42       27       43      10
## [2,]       5        8        2       33      37
## [3,]      11       32       35       17      26
## [4,]       2        4       32        2      18
## [5,]      31       18       21       40       1
```
为了找到w最近的邻居，我们使用了一个与之前类似的命令：

```
muraro_to_seger <- scmapCell(
  projection = muraro,
  index_list = list(
    seger = metadata(segerstolpe)$scmap_cell_index
  ),
  w = 5
)
```
查看结果：

```
muraro_to_seger$seger[[1]][,1:5]
##      D28.1_1 D28.1_13 D28.1_15 D28.1_17 D28.1_2
## [1,]    2201     1288     1117     1623    1078
## [2,]    1229     1724     2104     1448    1593
## [3,]    1793     1854     2201     2039    1553
## [4,]    1882     1737     1081     1202    1890
## [5,]    1731      976     1903     1834    1437
```
这显示了segerstolpe中5列的最近邻居到muraro中的每个cell。然后，我们可以通过从segerstolpe数据集的colData中选择适当的数据来计算伪时间估计，分支分配或其他单元级数据。 作为演示，我们将找到每个cell的最近邻居的cell类型。
##### 3.8.5 Metaneighbour
Metaneighbour专门用于询问细胞类型标签是否在数据集之间保持一致。它有两个版本。首先是一种完全监督的方法，它假设所有数据集中都知道细胞类型，并计算这些细胞类型标记的“好”。 （下面将描述“好”的确切含义）。或者，metaneighbour可以估计数据集内部和跨数据集的所有细胞类型彼此之间的相似程度。我们将只使用无监督版本，因为它具有更广泛的适用性，并且更容易解释结果。

Metaneighbour通过构建细胞 - 细胞spearman相关网络来比较数据集中的细胞类型。然后该方法尝试通过其最近邻居的加权“投票”来预测每个cell的标签。然后将两个聚类之间的整体相似性评分为AUROC，用于基于这些加权投票将类型A的cell分配给类型B. AUROC为1表示在任何其他细胞出现之前，所有类型A的细胞都被分配到typeB，如果细胞被随机分配，则AUROC为0.5。

Metanighbour只是几个R函数而不是一个完整的包，所以我们必须使用source加载它们

```
source("2017-08-28-runMN-US.R")
```
###### 3.8.5.1 Prepare Data
Metaneighbour要求在运行之前将所有数据集组合成single expression 矩阵：

```
is.common <- rowData(muraro)$feature_symbol %in% rowData(segerstolpe)$feature_symbol
muraro <- muraro[is.common,]
segerstolpe <- segerstolpe[match(rowData(muraro)$feature_symbol, rowData(segerstolpe)$feature_symbol),]
rownames(segerstolpe) <- rowData(segerstolpe)$feature_symbol
rownames(muraro) <- rowData(muraro)$feature_symbol
identical(rownames(segerstolpe), rownames(muraro))
## [1] TRUE
combined_logcounts <- cbind(logcounts(muraro), logcounts(segerstolpe))
dataset_labels <- rep(c("m", "s"), times=c(ncol(muraro), ncol(segerstolpe)))
cell_type_labels <- c(colData(muraro)$cell_type1, colData(segerstolpe)$cell_type1)

pheno <- data.frame(Sample_ID = colnames(combined_logcounts),
                Study_ID=dataset_labels,
                Celltype=paste(cell_type_labels, dataset_labels, sep="-"))
rownames(pheno) <- colnames(combined_logcounts)
```
Metaneighbor包含一种特征选择方法来识别高度可变的基因。

```
var.genes = get_variable_genes(combined_logcounts, pheno)
```
由于Metaneighbor比scmap慢得多，我们将对这些数据集进行抽样。
```
subset <- sample(1:nrow(pheno), 2000)
combined_logcounts <- combined_logcounts[,subset]
pheno <- pheno[subset,]
cell_type_labels <- cell_type_labels[subset]
dataset_labels <- dataset_labels[subset]
```
现在我们已经准备好运行metanbor。首先，我们将运行无人监督的版本，它将让我们看到两个数据集之间最相似的cell类型。

```
unsup <- run_MetaNeighbor_US(var.genes, combined_logcounts, unique(pheno$Celltype), pheno)
heatmap(unsup)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/31-projection_files/figure-html/unnamed-chunk-32-1.png)
##### 3.8.6 mnnCorrect
mnnCorrect纠正数据集以促进联合分析。它为了解释两个重复或两个不同实验之间的组成差异，它首先匹配实验中的单个细胞以找到重叠的生物结构。使用该重叠，它学习表达的哪个维度对应于生物状态以及哪个维度对应于批次/实验效果; mnnCorrect假设这些维度在高维表达空间中是彼此直系的。最后，它从整个表达式矩阵中删除批次/实验效果，以返回校正后的矩阵。
为了跨数据集匹配单个单元格，mnnCorrect使用余弦距离来避免库大小效应，然后识别跨数据集的互补最近邻居（k确定邻域大小）。只有重叠的生物群应该有相互最近的邻居（见下面的图b）。然而，这假设k被设置为大约数据集中最小生物群的大小，但是k太低将识别太少的相互最近邻居对以获得我们想要移除的批量效果的良好估计。
学习生物/技术效果的方法是使用奇异值分解，类似于我们在批量校正部分遇到的RUV，或者使用优化的irlba包进行主成分分析，这应该比SVD快。参数svd.dim指定应该保留多少维度来总结数据的生物结构，我们将它设置为三，因为我们发现使用上面的Metaneighbor的三个主要组。可以通过平滑（sigma）和/或方差调整（var.adj）进一步调整这些估计。

mnnCorrect还假设你已经将你的表达式矩阵子集化，以便它们以相同的顺序包含相同的基因，幸运的是，当我们为Metaneighbor设置数据时，我们已经为我们的数据集做了一些事情。
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/mnnCorrectDiagramCropped.png)

```
require("scran")
# mnnCorrect will take several minutes to run
corrected <- mnnCorrect(logcounts(muraro), logcounts(segerstolpe), k=20, sigma=1, pc.approx=TRUE, subset.row=var.genes, svd.dim=3)
```
首先让我们检查一下我们找到了足够数量的mnn对，mnnCorrect返回一个数据帧列表，其中包含每个数据集的mnn对。

```
dim(corrected$pairs[[1]]) # muraro -> others
## [1] 0 3
dim(corrected$pairs[[2]]) # seger -> others
## [1] 2533    3
```
第一列和第二列包含cell列ID，第三列包含指示第2列单元属于哪个数据集/批次的数字。 在我们的例子中，我们只比较两个数据集，因此所有mnn对都已分配给第二个表，第三个列只包含一个。

```
head(corrected$pairs[[2]])
## DataFrame with 6 rows and 3 columns
##   current.cell other.cell other.batch
##      <integer>      <Rle>       <Rle>
## 1         1553          5           1
## 2         1078          5           1
## 3         1437          5           1
## 4         1890          5           1
## 5         1569          5           1
## 6          373          5           1
total_pairs <- nrow(corrected$pairs[[2]])
n_unique_seger <- length(unique((corrected$pairs[[2]][,1])))
n_unique_muraro <- length(unique((corrected$pairs[[2]][,2])))
```

mnnCorrect在n_unique_seger segerstolpe cell和n_unique_muraro muraro cells之间找到了2533组相互最近邻居。 这应该是足够数量的对，但每个数据集中的独特细胞数量较少表明我们可能没有捕获每个数据集中的完整生物信号。
现在我们可以创建一个组合数据集来共同分析这些数据。 但是，更正的数据不再是count，通常包含负表达值，因此某些分析工具可能不再合适。 为简单起见，我们只绘制一个联合TSNE。

```
require("Rtsne")
## Loading required package: Rtsne
```

```
joint_expression_matrix <- cbind(corrected$corrected[[1]], corrected$corrected[[2]])

# Tsne will take some time to run on the full dataset
joint_tsne <- Rtsne(t(joint_expression_matrix[rownames(joint_expression_matrix) %in% var.genes,]), initial_dims=10, theta=0.75,
                        check_duplicates=FALSE, max_iter=200, stop_lying_iter=50, mom_switch_iter=50)
dataset_labels <- factor(rep(c("m", "s"), times=c(ncol(muraro), ncol(segerstolpe))))
cell_type_labels <- factor(c(colData(muraro)$cell_type1, colData(segerstolpe)$cell_type1))
plot(joint_tsne$Y[,1], joint_tsne$Y[,2], pch=c(16,1)[dataset_labels], col=rainbow(length(levels(cell_type_labels)))[cell_type_labels])
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/31-projection_files/figure-html/unnamed-chunk-38-1.png)
##### 3.8.7 Cannonical Correlation Analysis（Seurat）
Seurat包中包含另一种用于组合多个数据集的校正方法，称为CCA。但是，与mnnCorrect不同，它不直接更正表达式矩阵本身。相反，Seurat为每个数据集找到一个较低维度的子空间，然后校正这些子空间。与mnnCorrect不同，Seurat一次只组合一对数据集。

Seurat使用基因 - 基因相关性通过称为典型相关分析（CCA）的方法识别数据集中的生物结构。 Seurat了解基因 - 基因相关性的共享结构，然后评估每个细胞适合这种结构的程度。假设必须通过数据特定维数减少方法而不是共享相关结构更好地描述的单元表示数据集特定的单元类型/状态，并且在对齐两个数据集之前被丢弃。最后，使用“扭曲”算法对两个数据集进行对齐，该算法以对人口密度差异稳健的方式对每个数据集的低维表示进行归一化。

请注意，因为Seurat占用了大量的库空间，您必须重新启动R会话才能加载它，并且不会在此页面上自动生成plots / merged_seuratput。
载入数据：

```
muraro <- readRDS("pancreas/muraro.rds")
segerstolpe <- readRDS("pancreas/segerstolpe.rds")
segerstolpe <- segerstolpe[,colData(segerstolpe)$cell_type1 != "unclassified"]
segerstolpe <- segerstolpe[,colData(segerstolpe)$cell_type1 != "not applicable",]
muraro <- muraro[,colData(muraro)$cell_type1 != "unclear"]
is.common <- rowData(muraro)$feature_symbol %in% rowData(segerstolpe)$feature_symbol
muraro <- muraro[is.common,]
segerstolpe <- segerstolpe[match(rowData(muraro)$feature_symbol, rowData(segerstolpe)$feature_symbol),]
rownames(segerstolpe) <- rowData(segerstolpe)$feature_symbol
rownames(muraro) <- rowData(muraro)$feature_symbol
identical(rownames(segerstolpe), rownames(muraro))
```
首先，我们将把数据重新格式化为Seurat对象：
```
require("Seurat")
set.seed(4719364)
muraro_seurat <- CreateSeuratObject(raw.data=assays(muraro)[["normcounts"]]) # raw counts aren't available for muraro
muraro_seurat@meta.data[, "dataset"] <- 1
muraro_seurat@meta.data[, "celltype"] <- paste("m",colData(muraro)$cell_type1, sep="-")

seger_seurat <- CreateSeuratObject(raw.data=assays(segerstolpe)[["counts"]])
seger_seurat@meta.data[, "dataset"] <- 2
seger_seurat@meta.data[, "celltype"] <- paste("s",colData(segerstolpe)$cell_type1, sep="-")
```
接下来，我们必须对每个数据集进行标准化、缩放和识别高度可变的基因：

```
muraro_seurat <- NormalizeData(object=muraro_seurat)
muraro_seurat <- ScaleData(object=muraro_seurat)
muraro_seurat <- FindVariableGenes(object=muraro_seurat, do.plot=TRUE)

seger_seurat <- NormalizeData(object=seger_seurat)
seger_seurat <- ScaleData(object=seger_seurat)
seger_seurat <- FindVariableGenes(object=seger_seurat, do.plot=TRUE)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/muraro_seurat_hvg.png)
         muraro variable genes
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/seger_seurat_hvg.png)
          segerstolpe variable genes
虽然Seurat校正了dispersion和平均表达之间的关系，但在排序特征时它不会使用校正值。 将以下命令的结果与上图中的结果进行比较：

```
head(muraro_seurat@hvg.info, 50)
head(seger_seurat@hvg.info, 50)
```
但是我们将按照他们的例子使用2000个最分散withmerged_seurat的基因来校正每个数据集的平均表达式。

```
gene.use <- union(rownames(x = head(x = muraro_seurat@hvg.info, n = 2000)),
                  rownames(x = head(x = seger_seurat@hvg.info, n = 2000)))
```
现在我们将运行CCA来查找这两个数据集的共享相关结构：

请注意，为了加快计算速度，我们将仅使用前5个维度，但理想情况下，您会考虑更多维度，然后使用DimHeatmap选择信息量最大的维度。

```
merged_seurat <- RunCCA(object=muraro_seurat, object2=seger_seurat, genes.use=gene.use, add.cell.id1="m", add.cell.id2="s", num.cc = 5)
DimPlot(object = merged_seurat, reduction.use = "cca", group.by = "dataset", pt.size = 0.5) # Before correcting
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/cca_before.png)
为了识别特定于数据集的细胞类型，我们比较了CCA与数据集特定主成分分析“解释”细胞的程度。

```
merged_seurat <- CalcVarExpRatio(object = merged_seurat, reduction.type = "pca", grouping.var = "dataset", dims.use = 1:5)
merged.all <- merged_seurat
merged_seurat <- SubsetData(object=merged_seurat, subset.name="var.ratio.pca", accept.low = 0.5) # CCA > 1/2 as good as PCA
merged.discard <- SubsetData(object=merged.all, subset.name="var.ratio.pca", accept.high = 0.5)

summary(factor(merged.discard@meta.data$celltype)) # check the cell-type of the discarded cells.
```
在这里我们可以看到，尽管两个数据集都包含内皮细胞，但几乎所有数据集都被丢弃为“数据集特定”。 现在我们可以对齐数据集：

```
merged_seurat <- AlignSubspace(object = merged_seurat, reduction.type = "cca", grouping.var = "dataset", dims.align = 1:5)
DimPlot(object = merged_seurat, reduction.use = "cca.aligned", group.by = "dataset", pt.size = 0.5) # After aligning subspaces
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/cca_after.png)
#### 3.9 Search scRNA-Seq data
scfind是一种工具，它允许人们使用基因列表搜索单个细胞RNA-Seq集合（Atlas），例如搜索特定的基因集合的细胞和细胞类型。
![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/scfind.png)
##### 3.9.2 Dataset

```
muraro <- readRDS("pancreas/muraro.rds")
```
##### 3.9.3 Gene INdex
使用刚刚加载的数据创建索引：

```
cellIndex <- buildCellIndex(muraro)
```
基因索引包含表达它的细胞的每个基因索引。 这类似于表达矩阵的稀疏化。 除此之外，索引也以一种可以非常快速地访问的方式进行压缩。 我们估计用这种方法可以达到5-10的压缩系数。

默认情况下，SingleCellExperiment对象的colData槽的cell_type1列用于定义单元类型，但也可以使用buildCellTypeIndex函数的cell_type_column参数手动定义
##### 3.9.4 Marker genes
现在让我们定义细胞类型特异性标记基因的列表。我们将使用在原始出版物中确定的标记基因：

```
# these genes are taken from fig. 1
muraro_alpha <- c("GCG", "LOXL4", "PLCE1", "IRX2", "GC", "KLHL41", 
                  "CRYBA2", "TTR", "TM4SF4", "RGS4")
muraro_beta <- c("INS", "IAPP", "MAFA", "NPTX2", "DLK1", "ADCYAP1", 
                 "PFKFB2", "PDX1", "TGFBR3", "SYT13")
muraro_gamma <- c("PPY", "SERTM1", "CARTPT", "SLITRK6", "ETV1", 
                  "THSD7A", "AQP3", "ENTPD2", "PTGFR", "CHN2")
muraro_delta <- c("SST", "PRG4", "LEPR", "RBP4", "BCHE", "HHEX", 
                  "FRZB", "PCSK1", "RGS2", "GABRG2")
```
##### 3.9.5 Search cells by a gene list
findCell函数返回与给定数据集中所有cell类型对应的p值列表。 它还输出一组细胞，其中来自给定基因列表的基因共表达。 我们将在上面定义的所有标记基因列表上运行它：

```
res <- findCell(cellIndex, muraro_alpha)
barplot(-log10(res$p_values), ylab = "-log10(pval)", las = 2)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/32-search_files/figure-html/unnamed-chunk-7-1.png)
### 4 Seurat
Seurat最初是作为scRNA-seq数据的聚类工具开发的，但是在过去几年中，该软件包的重点变得不那么具体，目前Seurat是一个流行的R软件包，可以执行QC，分析和scRNA的探索 -seq数据，即本课程涵盖的许多任务。 虽然作者提供了几个教程，但在这里我们按照Seurat的作者（2,800个外周血单核细胞）创建的例子提供了一个简要的概述。 我们主要在各种函数调用中使用默认值，有关详细信息，请参阅文档和作者。 出于课程目的，将使用前面章节中描述的小邓数据集：

```
deng <- readRDS("deng/deng-reads.rds")
```
#### 4.1 Seurat object class
Seurat没有集成上面描述的SingleCellExperiment Bioconductor类，而是引入了自己的对象类 - seurat。 本章中的所有计算都是在此类的对象上执行的。 要开始分析，我们首先需要使用原始（非规范化）数据初始化对象。 我们将保留所有表达>=3的基因和所有检测到至少200个基因的细胞。

```
library(SingleCellExperiment)
library(Seurat)
library(mclust)
library(dplyr)
seuset <- CreateSeuratObject(
    raw.data = counts(deng),
    min.cells = 3, 
    min.genes = 200
)
```
#### 4.2 Expression QC
Seurat允许您根据任何用户定义的标准轻松地探索QC度量和过滤cell。我们可以想象基因和分子的数量并绘制它们的关系：

```
VlnPlot(
    object = seuset, 
    features.plot = c("nGene", "nUMI"), 
    nCol = 2
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/33-seurat_files/figure-html/unnamed-chunk-5-1.png)

```
GenePlot(
    object = seuset, 
    gene1 = "nUMI", 
    gene2 = "nGene"
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/33-seurat_files/figure-html/unnamed-chunk-5-2.png)
现在，我们将排除具有明显离群值的细胞：

```
seuset <- FilterCells(
    object = seuset, 
    subset.names = c("nUMI"), 
    high.thresholds = c(2e7)
)
```
#### 4.3 Normalization
从数据集中删除不需要的单元格后，下一步是规范化数据。 默认情况下，我们采用全局缩放规范化方法LogNormalize，通过总表达式对每个单元格的基因表达式测量值进行标准化，将其乘以比例因子（默认为10,000），并对结果进行对数转换：

```
seuset <- NormalizeData(
    object = seuset, 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
)
```
#### 4.4 Highly Variable genes
Seurat计算高度可变的基因，并将重点放在这些基因上进行下游分析。findvariable基因计算每个基因的平均表达和色散，将这些基因放入bin中，然后计算每个bin内的分散度。这有助于控制可变性和平均表达之间的关系：
```
seuset <- FindVariableGenes(
    object = seuset,
    mean.function = ExpMean, 
    dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, 
    x.high.cutoff = 3, 
    y.cutoff = 0.5
)
```

```
length(x = seuset@var.genes)
## [1] 6127
```
#### 4.5 Dealing with confounders
为了减轻混杂因素的影响，Seurat构建线性模型以基于用户定义的变量预测基因表达。 这些模型的缩放z得分残差存储在scale.data槽中，用于降维和聚类。

Seurat可以消除由批次，细胞排列率（由Drop-seq数据的Drop-seq工具提供），检测到的分子数量，线粒体基因表达和细胞周期驱动的基因表达中的细胞 - 细胞变异。 在这里，我们回归每个细胞检测到的分子数量。

```
seuset <- ScaleData(
    object = seuset, 
    vars.to.regress = c("nUMI")
)
```
#### 4.6 Linear dimensionality reduction
接下来，我们对缩放数据执行PCA。 默认情况下，object@var.genes中的基因用作输入，但也可以使用pc.genes定义。 在高度可变的基因上运行降维可以提高性能。 然而，对于某些类型的数据（UMI） - 特别是在回归技术变量之后，当在更大的基因子集（包括整个转录组）上运行时，PCA返回相似（尽管更慢）的结果。

```
seuset <- RunPCA(
    object = seuset, 
    pc.genes = seuset@var.genes, 
    do.print = TRUE, 
    pcs.print = 1:5, 
    genes.print = 5
)
```
Seurat提供了几种有用的方法来可视化定义PCA的细胞和基因：
```
PrintPCA(object = seuset, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
```

```
VizPCA(object = seuset, pcs.use = 1:2)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/33-seurat_files/figure-html/unnamed-chunk-11-1.png)

```
PCAPlot(object = seuset, dim.1 = 1, dim.2 = 2)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/33-seurat_files/figure-html/unnamed-chunk-11-2.png)
特别是，PCHeatmap允许轻松探索数据集中异质性的主要来源，并且在尝试确定要包含哪些PC以进行进一步的下游分析时非常有用。 细胞和基因都根据其PCA分数排序。 将cells.use设置为数字可以绘制光谱两端的极端单元格，从而大大加快了大型数据集的绘制速度：
```
PCHeatmap(
    object = seuset, 
    pc.use = 1:6, 
    cells.use = 500, 
    do.balanced = TRUE, 
    label.columns = FALSE,
    use.full = FALSE
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/33-seurat_files/figure-html/unnamed-chunk-12-1.png)
#### 4.7 Significant PCs
为了克服任何单个基因中scRNA-seq数据的广泛技术噪音，Seurat根据其PCA评分对细胞进行聚类，每个PC基本上代表一种结合了相关基因集的信息的metagene。 因此，确定下游包含多少PC是重要的一步。 Seurat随机置换数据的一个子集（默认为1％）并重新运行PCA，通过重复此过程构建基因分数的空分布。 我们将重要的PC识别为具有低p值基因强烈富集的PC：
```
seuset <- JackStraw(
    object = seuset, 
    num.replicate = 100, 
    do.print = FALSE
)
```
JackStrawPlot函数提供了一个可视化工具，用于比较每个具有均匀分布（虚线）的PC的p值分布。 重要的PC将显示具有低p值的基因的强烈富集（虚线上方的实线）。 在这种情况下，PC 1-8似乎很重要。
```
JackStrawPlot(object = seuset, PCs = 1:9)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/33-seurat_files/figure-html/unnamed-chunk-14-1.png)
一种更特别的方法来确定要使用哪些PCs，看一看主成分的标准偏差的图，并在图中有一个清晰的弯头的地方画出你的截止点。这可以用PCElbowPlot来完成。在这个例子中，看起来拐点会落在PC 5周围。

```
PCElbowPlot(object = seuset)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/33-seurat_files/figure-html/unnamed-chunk-15-1.png)
#### 4.8 Clustering cells
Seurat实现了基于图形的聚类方法。基于先前识别的PC计算cell之间的距离。 Seurat方法受到最近manuscript的启发，这些manuscript将基于图形的聚类方法应用于scRNA-seq数据--SNN-Cliq（（C。Xu和Su 2015））和CyTOF数据 - PhenoGraph（（Levine等人，2015））。简而言之，这些方法将细胞嵌入图结构中 - 例如K-最近邻（KNN）图，在具有相似基因表达模式的细胞之间绘制边缘，然后尝试将该图划分为高度互连的准集团或群落。与PhenoGraph一样，我们首先根据PCA空间中的欧氏距离构建KNN图，并根据其局部邻域中的共享重叠（Jaccard距离）细化任意两个cell之间的边权重。为了聚类细胞，我们应用模块化优化技术 - SLM（（Blondel等人，2008）），将细胞迭代分组在一起，目的是优化标准模块化功能。

FindClusters函数实现该过程，并包含一个分辨率参数，用于设置下游群集的粒度，增加的值会导致更多的群集。我们发现在这之间设置此参数
0.6- 1.2通常会为大约包含3000细胞的单细胞数据集返回良好的结果。较大的数据集通常会增加最佳分辨率。群集保存在object@ ident中。

```
seuset <- FindClusters(
    object = seuset, 
    reduction.type = "pca", 
    dims.use = 1:8, 
    resolution = 1.0, 
    print.output = 0, 
    save.SNN = TRUE
)
```
Seurat的一个有用特性是能够回忆起在最新的函数调用中使用的参数。对于findcluster来说，有一个函数PrintFindClustersParams来打印一个格式良好的参数总结，这些参数是选择的：

```
PrintFindClustersParams(object = seuset)
## Parameters used in latest FindClusters calculation run on: 2018-05-29 15:27:42
## =============================================================================
## Resolution: 1
## -----------------------------------------------------------------------------
## Modularity Function    Algorithm         n.start         n.iter
##      1                   1                 100             10
## -----------------------------------------------------------------------------
## Reduction used          k.param          k.scale          prune.SNN
##      pca                 30                25              0.0667
## -----------------------------------------------------------------------------
## Dims used in calculation
## =============================================================================
## 1 2 3 4 5 6 7 8
```
我们可以查看集群结果，并将它们与原始的细胞标签进行比较：
```
table(seuset@ident)
## 
##  0  1  2  3 
## 85 75 59 34
adjustedRandIndex(colData(deng)[seuset@cell.names, ]$cell_type2, seuset@ident)
## [1] 0.3981315
```
Seurat还利用tSNE图来识别聚类结果。 作为对tSNE的输入，我们建议使用作为聚类分析的输入的相同PC，尽管使用genes.use参数也支持基于缩放基因表达计算tSNE。
```
seuset <- RunTSNE(
    object = seuset,
    dims.use = 1:8,
    do.fast = TRUE
)
TSNEPlot(object = seuset)

```
![image](https://hemberg-lab.github.io/scRNA.seq.course/33-seurat_files/figure-html/unnamed-chunk-19-1.png)
#### 4.9 Marker genes
Seurat可以帮助您找到通过差异表达定义聚类的标记。 默认情况下，它标识单个群集的正负标记，与所有其他单元格进行比较。 您可以测试群集彼此之间或所有cell。 例如，要找到第2组的标记基因，我们可以运行：

```
markers2 <- FindMarkers(seuset, 2)
```
然后，标记基因可以被可视化：

```
VlnPlot(object = seuset, features.plot = rownames(markers2)[1:2])
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/33-seurat_files/figure-html/unnamed-chunk-21-1.png)
```
FeaturePlot(
    seuset, 
    head(rownames(markers2)), 
    cols.use = c("lightgrey", "blue"), 
    nCol = 3
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/33-seurat_files/figure-html/unnamed-chunk-21-2.png)
FindAllMarkers自动执行此过程并查找所有群集的标记
```
markers <- FindAllMarkers(
    object = seuset, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    thresh.use = 0.25
)
```
DoHeatmap为给定的细胞和基因生成表达式热图。 在这种情况下，我们为每个群集绘制前10个标记（或所有标记，如果小于20）：
```
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(
    object = seuset, 
    genes.use = top10$gene, 
    slim.col.label = TRUE, 
    remove.key = TRUE
)
```
![image](https://hemberg-lab.github.io/scRNA.seq.course/33-seurat_files/figure-html/unnamed-chunk-23-1.png)
### 5 “Ideal” scRNA-seq pipline
#### 5.1 Experimental Design
- 避免混淆生物和批处理效果
  - 如果可能的话，应该在相同的芯片上捕获多个条件
  - 对每个条件执行多个复制，如果可能的话，应该在一起执行不同条件的复制
  - 统计数字不能纠正一个完全令人困惑的实验！
- 独特的分子标识符
  - 大大减少数据中的噪声
  - 可以降低基因检测率（不清楚是UMIs还是其他操作差异）
  - 失去连接信息
  - 使用更长的umi(~ 10 bp)
  - 使用UMI-tools对UMIs中的排序错误进行校正
- Spike-ins
  - 有利于质量控制
  - 可能对标准化读取计数很有用
  - 可用于近似细胞大小/ RNA含量（如果与生物学问题相关）
  - 通常表现出比内源基因更高的噪音（移液错误，混合物质量）
  - 需要更多的测序才能获得每个细胞足够的内源性读数
- 细胞数量 vs read depth
  - 基因检测平台从每个细胞1百万个读数开始
  - 转录因子检测（调节网络）需要高读取深度和最敏感的方案（即Fluidigm C1）
  - 细胞聚类和细胞类型鉴定受益于大量细胞，并且不需要高测序深度（每个细胞约100,000个读数）。
 ![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/Pipeline-batches.png)
#### 5.2 Processing Reads
- Read QC & Trimming
  - FASTQC，cutadapt
- Mapping
  - Small datasets or UMI datasets：align to genome/transcriptome using STAR
  - Large datasets :pseudo-alignment with Salmon or kallisto
- Quantification
  - Small dataset,no UMI:featureCounts
  - Large dataset,no UMI:Salmon,kallisto
  - UMI dataset:UMI-tools+featureCounts
#### 5.3 Preparing Expression Matrix
- Cell QC
  - scater
  - consider:mtRNA,rRNA,spikr-ins(if available),number of detected genes per cell,total reads/molecules per cell
- Library Size Normalization
  - scran
- Batch correction(if appropriate)
  - Replicates/Confounded RUVs
  - Unknown or unbalanced biological groups mnnCorrect
  - Balanced sedign Combat
#### 5.4 Biological Interpretation
- Feature Selection
  - M3Drop
- Clustering and Marker Gene Identification
  - <=5000 cells:SC3
  - > 5000 cells:Seurat
- Pseudotime 
  - distinct timepoints:TSCAN
  - small dataset/unknown number of branches:Monocle2
  - large continuous dataset:destiny
- Differential Expression
  - Small number of cells and few group:scde
  - Replicates with bath effects:mixture/linear models
  - Balanced batchs:edgeR or MAST
  - Large dataset:Kruskal-Wallix test(all groups a once),or Wilcox-test(compare 2-groups at a time)
 