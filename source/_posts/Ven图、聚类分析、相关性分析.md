---
title: Ven图、聚类分析、相关性分析
tag: [Ven图,聚类分析,相关性分析]
categories: R
---


##### 绘制Ven图
```R
  #载入相应的安装包
  library(VennDiagram)
  #进入到需要绘制Ven图的文件夹
  setwd("/data/YJ/YanNa/stringtie/trs/Pt/diff10/LPt")
  #filelist=list("Pt15","Pt17","Pt19","Pt2-","Pt20","Pt22","Pt3","Pt5","Pt6","Pt8","Pt9")
  #将需要绘制Ven图的样本的rowname汇集到一个list中
  filelist=list("LPt16","LPt24")
  for (Pt in filelist){
    files=list.files(path="/data/YJ/YanNa/stringtie/trs/Pt/diff10/LPt",pattern = Pt)
    trslist=list()
    image_name=paste(Pt,"_Ven",".tiff",sep="")
    for (i in files){
      base=strsplit(i,"_")[[1]][1]
      trslist[base]=read.table(i,header = T,colClasses = c("character","NULL","NULL"))
    }
    #设置颜色
    num=148
    col=c()
    for (i in 1:length(trslist)){
      col[i]=colors()[num]
      num=num+150
    }
    #绘制Ven图
    venn.diagram(trslist,filename = image_name,col="transparent",fil=col,alpha=0.5,resolution = 500)
    
```
<!--more-->
##### 聚类分析
###### K均值聚类
K均值聚类又称为动态聚类，它的计算方法较为简单，也不需要输入距离矩阵。首先要指定聚类的分类个数N，随机取N个样本作为初始类的中心，计算各样本与类中心的距离并进行归类，所有样本划分完成后重新计算类中心，重复这个过程直到类中心不再变化。 
###### 二、层次聚类
层次聚类又称为系统聚类，首先要定义样本之间的距离关系，距离较近的归为一类，较远的则属于不同的类。可用于定义“距离”的统计量包括了欧氏距离(euclidean)、马氏距离(manhattan)、 两项距离(binary)、明氏距离(minkowski)。还包括相关系数和夹角余弦。

层次聚类首先将每个样本单独作为一类，然后将不同类之间距离最近的进行合并，合并后重新计算类间距离。这个过程一直持续到将所有样本归为一类为止。在计算类间距离时则有六种不同的方法，分别是最短距离法、最长距离法、类平均法、重心法、中间距离法、离差平方和法。 

```
==!!注意==
K均值聚类能处理比层次聚类更大的数据集。由于K均值聚类在开始要随机选择k个中心点，在每次调用函数时可能获得不同的方案。使用
set.seed() 函数可以保证结果是可复制的。此外，聚类方法对初始中心值的选择也很敏感。

```
==聚类分析的一般步骤==  

1. 选择合适的变量，这是第一步，也可能是最重要的一步，再高级的聚类方法也不能弥补聚类变量选择不好的问题。因此，首先选择你感觉可能对于识别和理解数据中不同观测值分组有重要影响的变量。  
2. 缩放数据，由于不同变量可能有不同的变化范围，以免那些变化范围大的变量对结果有不成比例的影响，常常需要在分析之前缩放数据。通常，将每个变量标准化为均值为0和标准差为1的变量。还有另外两种替代方法，比如每个变量被其最大值相除，或该变量减去它的平均值并除以变量的平均绝对偏差。实现代码如下：
```
df1 <- apply(mydata, 2, function(x){(x-mean(x))/sd(x)})
df2 <- apply(mydata, 2, function(x){x/max(x)})
df3 <- apply(mydata, 2, function(x){(x – mean(x))/mad(x)})
```
在R中可以通过使用scale()函数实现df1代码片段的功能。
3. 寻找异常点，通常聚类方法对于异常值比较敏感，对此有几种解决方法，通过outliers包中的函数来筛选异常单变量离群点，mvoutlier包中 有识别多元变量的离群点的函数。  
4. 计算距离，虽然不同的聚类算法差异很大，但是通常需要计算被聚类的实体之间的距离。两个观测之间最常用的距离度量是欧几里得距离，其它常用的还有曼哈顿距离、兰氏距离、非对称二元距离、最大距离和闵可夫斯基距离。在下文的论述当中默认采用欧几里得距离。  
两个观测值之间的欧几里得距离定义如下：
![image](https://www.zhihu.com/equation?tex=d_%7Bij%7D+%3D%5Csqrt%7B%5Csum_%7Bp%3D1%7D%5E%7Bp%7D%7B+%5Cleft%28+x_%7Bip%7D+-x_%7Bjp%7D+%5Cright%29%5E%7B2%7D+%7D+%7D+)  
用R中自带的dist()函数可以计算矩阵或数据框中所有行（观测值）之间的距离，dist()函数的格式为dist(x, method =)，默认为欧几里得距离。

观测值之间的距离越大，则异质性越大，距离也越远。

通常，欧几里得距离作为连续型数据的距离度量，对于其它类型的数据，可以使用cluster包中的daisy()函数获得包含任意二元、名义、有序或者连续属性组合的相异矩阵，并且cluster包中的其它函数可以用这种异质性进行聚类分析。

另外，当一个观测中的某一个变量变换范围太大，缩放数据有利于均衡各变量的影响。
5. 选择聚类算法，对于不同类型和不同样本量的数据，针对性地选择不同的聚类算法，从大的方面讲有层次聚类和划分聚类，通常前者对于小样本来说很实用（如150个观测值或更少），而且这种情况下嵌套聚类更实用；后者能够处理更大的数据量，但是需要事先确定聚类的个数。然后，不管是层次聚类还是划分聚类，都需要一个特定的聚类算法，不同算法有着各自的优缺点，需要根据工程实践进行综合比选确定。  
6. 获 得一种或多种聚类方法，这是步骤5的延续。  
7. 确定类的数目，在确定最终聚类方案时必须确定类的数目。通常尝试不同的类型（比如2~K）并比较解的质量。NbClust包中的NbClust()函数30个不同的指标来帮助你进行选择，指标如此之多，可见这在聚类分析当中是个难题。
8. 获得最终的聚类解决方案，类的个数确定以后，就可以提取出子群，形成最终的聚类方案。  
8. 结果可视化，可视化可以帮助判定聚类方案的意义和用处，层次聚类的结果通常表示为一个树状图。划分的结果通常利用可视化双变量聚类图来表示。  
9. 解读类，聚类方案确定以后，结果也出来了，必须对其进行解读。比如，一个类中的观测值有何相似之处，不同类之间有何不同，这一步通常通过获得类中每个变量的汇总统计来完成。对于连续数据，每一类中变量的均值和中位数会被计算出来。对于混合数据（数据中包含分类变量），结果中将返回各类的众数或类别分布。  
20. 验证结果，对于聚类方案的验证相当于确认下，如果采用不同的聚类方法或者不同样本，是否会产生相同的类。fpc、clv和clValid包包含了评估聚类解的稳定性的函数。


```R
#安装并加载包
  #使用k-means聚类所需要的包：factoextra和cluster
  site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
  pakeags_list=c("factoextra","cluster","NbClust")
  for(p in pakeags_list){
    if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
      install.packages(p, repos=site)
      suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
    }
  }
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("ConsensusClusterPlus", version = "3.8")
  #载入相关的程序包
  library(ConsensusClusterPlus)
  library(factoextra)
  library(cluster)
  library(NbClust)
  library(reshape2)
```

```

  #载入数据
  data<-read.table("Pt_RPt_LW.txt",header=T)
  #根据rownames去掉数据中的重复行
  index<-duplicated(data$trs)
  data<-data[!index,]
  rownames(data)<-data$trs
  data<-data[,-1]
  #t()函数对数据进行转置，形成，行为基因名，列为样本名称的数据集
  tdata<-as.matrix(t(data))
  #去掉数据集中的NA函数
  my_data<-na.omit(tdata)
  
  #在聚类之前，进行不要的数据检查即数据描述性统计，如平均值、标准差
  stats=data.frame(Min=apply(data,2,min),
                   Med=apply(data,2,median),
                   Mean=apply(data,2,mean),
                   SD=apply(data,2,sd),
                   Max=apply(data,2,max)
    
  )
  stats=round(stats,1)#保留小数点最后一位
  
  #数据标准化和评估
  #变量有很大的方差及均值时需要进行标准化
  df=scale(my_data)
  #数据集群性评估，使用get_clust_tendency()计算Hopkins统计量
  res=get_clust_tendency(df,40,graph=TRUE)
  res$hopkins_stat #Hopkins统计量的值<0.5,表明数据是高度可聚合的
  res$plot
```


###### 使用Kmean函数
算法的简单描述如下：
1. 选择K个中心点（随机选择K行）；  
1. 把每个数据点分配到离它最近的中心点； 
1. 重新计算每类中的点到该类中心点距离的平均值；
1. 分配每个数据到它最近的中心点； 
1. 重复步骤3和步骤4，直到所有的观测值不再被分配或是达到最大的迭代次数。  

从算法上可以看出，初始中心值的选择对于结果也非常敏感，kmeans()函数有一个nstart选项可以尝试多种初始配置并输出最好的一个结果。举个例子，加上nstart=25会生成25个初始配置，一般来说，我们推荐使用这种方法。
```R
#计算最适合的K-mean聚类的分类数
  set.seed(123)
  gap_stat=clusGap(df,FUN=kmeans,nstart=25,K.max = 10,B=500)
  fviz_gap_stat(gap_stat)
  km.res=kmeans(df,2,nstart = 25)
  fviz_cluster(km.res,df)
```
###### eclust():增强的聚类分析  
与其他聚类分析包相比，eclust()有以下优点：  
- 简化了聚类分析的工作流程，
- 可以用于计算层次聚类和分区聚类，
- eclust()自动计算最佳聚类簇数。 
- 自动提供Silhouette plot，可以结合ggplot2绘制优美的图形
```R
#使用eclust()的K均值聚类
res.km<-eclust(df,"kmeans")
fviz_gap_stat(res.km$gap_stat)
#利用eclust()做层次聚类
res.hc=eclust(df,"hclust")
fviz_dend(res.hc,rect = TRUE)
```
###### 利用hclust()做层次聚类
1. 层次聚类的算法如下：
1. 定义每一个观测值（行或单元）为一类；
1. 计算每类和其他各类的距离；
1. 把距离最短的两类合并成一类，减少一个类的个数；
1. 重复步骤2和3，直到包含所有观测值的类合并成单个类为止。  
 
不同的层次聚类算法主要区别在于它们对类的定义不同，常见的五种聚类方法及其距离定义如下
![image](https://pic3.zhimg.com/80/v2-e6a84d52ffeaa71cfc2947d8141a0932_hd.png)  
层次聚类方法的实现格式如下：  

```
hclust(d, method=)
```

其中d是通过dist() 函数产生的距离矩阵， 并且方法包括 "single" 、"complete" 、"average" 、"centroid"和"ward"。
```
#先求样本之间两两相似性
result<-dist(df,method="euclidean")
#产生层次结构
result_hc<-hclust(d=result,method = "ward.D2")
#初步结果
fviz_dend(result_hc,cex=0.3)
#可视化展示
fviz_dend(result_hc,k=2,cex=0.3,k_colors = c("#E7B800", "#00AFBB"),color_labels_by_k = TRUE,rect = TRUE)
```

###### 利用NbCluster()函数首先选择最适合的分类k
```
##利用NbCluster()函数首先选择最适合的分类k
  res_nbclust<-NbClust(df,distance = "euclidean",min.nc = 2,max.nc = 10,method = "complete",index = "all")
```


##### 相关性分析
相关系数的显著性水平
使用Hmisc包，计算矩阵相关系数及其对应的显著性水平

```R
library(Hmisc)
res<-rcorr(as.matrix(mtcars))
res
       mpg   cyl  disp    hp  drat    wt  qsec    vs    am  gear  carb
mpg   1.00 -0.85 -0.85 -0.78  0.68 -0.87  0.42  0.66  0.60  0.48 -0.55
cyl  -0.85  1.00  0.90  0.83 -0.70  0.78 -0.59 -0.81 -0.52 -0.49  0.53
disp -0.85  0.90  1.00  0.79 -0.71  0.89 -0.43 -0.71 -0.59 -0.56  0.39
hp   -0.78  0.83  0.79  1.00 -0.45  0.66 -0.71 -0.72 -0.24 -0.13  0.75
drat  0.68 -0.70 -0.71 -0.45  1.00 -0.71  0.09  0.44  0.71  0.70 -0.09
wt   -0.87  0.78  0.89  0.66 -0.71  1.00 -0.17 -0.55 -0.69 -0.58  0.43
qsec  0.42 -0.59 -0.43 -0.71  0.09 -0.17  1.00  0.74 -0.23 -0.21 -0.66
vs    0.66 -0.81 -0.71 -0.72  0.44 -0.55  0.74  1.00  0.17  0.21 -0.57
am    0.60 -0.52 -0.59 -0.24  0.71 -0.69 -0.23  0.17  1.00  0.79  0.06
gear  0.48 -0.49 -0.56 -0.13  0.70 -0.58 -0.21  0.21  0.79  1.00  0.27
carb -0.55  0.53  0.39  0.75 -0.09  0.43 -0.66 -0.57  0.06  0.27  1.00

n= 32 


P
     mpg    cyl    disp   hp     drat   wt     qsec   vs     am     gear  
mpg         0.0000 0.0000 0.0000 0.0000 0.0000 0.0171 0.0000 0.0003 0.0054
cyl  0.0000        0.0000 0.0000 0.0000 0.0000 0.0004 0.0000 0.0022 0.0042
disp 0.0000 0.0000        0.0000 0.0000 0.0000 0.0131 0.0000 0.0004 0.0010
hp   0.0000 0.0000 0.0000        0.0100 0.0000 0.0000 0.0000 0.1798 0.4930
drat 0.0000 0.0000 0.0000 0.0100        0.0000 0.6196 0.0117 0.0000 0.0000
wt   0.0000 0.0000 0.0000 0.0000 0.0000        0.3389 0.0010 0.0000 0.0005
qsec 0.0171 0.0004 0.0131 0.0000 0.6196 0.3389        0.0000 0.2057 0.2425
vs   0.0000 0.0000 0.0000 0.0000 0.0117 0.0010 0.0000        0.3570 0.2579
am   0.0003 0.0022 0.0004 0.1798 0.0000 0.0000 0.2057 0.3570        0.0000
gear 0.0054 0.0042 0.0010 0.4930 0.0000 0.0005 0.2425 0.2579 0.0000       
carb 0.0011 0.0019 0.0253 0.0000 0.6212 0.0146 0.0000 0.0007 0.7545 0.1290
     carb  
mpg  0.0011
cyl  0.0019
disp 0.0253
hp   0.0000
drat 0.6212
wt   0.0146
qsec 0.0000
vs   0.0007
am   0.7545
gear 0.1290
carb

```
提取矩阵相关及其P值

```
CorMatrix <- function(cor,p) {
+                               ut <- upper.tri(cor) 
+                               data.frame(row = rownames(cor)[row(cor)[ut]] ,
+                               column = rownames(cor)[col(cor)[ut]], 
+                               cor =(cor)[ut], 
+                               p = p[ut] )
+ }
res <- rcorr(as.matrix(mtcars))
> CorMatrix (res$r, res$P)
    row column         cor            p
1   mpg    cyl -0.85216196 6.112688e-10
2   mpg   disp -0.84755138 9.380328e-10
3   cyl   disp  0.90203287 1.803002e-12
4   mpg     hp -0.77616837 1.787835e-07
5   cyl     hp  0.83244745 3.477861e-09
6  disp     hp  0.79094859 7.142679e-08
7   mpg   drat  0.68117191 1.776240e-05
8   cyl   drat -0.69993811 8.244636e-06
9  disp   drat -0.71021393 5.282022e-06
10   hp   drat -0.44875912 9.988772e-03
11  mpg     wt -0.86765938 1.293958e-10
12  cyl     wt  0.78249579 1.217567e-07
13 disp     wt  0.88797992 1.222311e-11
14   hp     wt  0.65874789 4.145827e-05
15 drat     wt -0.71244065 4.784260e-06
16  mpg   qsec  0.41868403 1.708199e-02
17  cyl   qsec -0.59124207 3.660533e-04
18 disp   qsec -0.43369788 1.314404e-02
19   hp   qsec -0.70822339 5.766253e-06
20 drat   qsec  0.09120476 6.195826e-01
21   wt   qsec -0.17471588 3.388683e-01
22  mpg     vs  0.66403892 3.415937e-05
23  cyl     vs -0.81081180 1.843018e-08
24 disp     vs -0.71041589 5.235012e-06
25   hp     vs -0.72309674 2.940896e-06
26 drat     vs  0.44027846 1.167553e-02
27   wt     vs -0.55491568 9.798492e-04
28 qsec     vs  0.74453544 1.029669e-06
29  mpg     am  0.59983243 2.850207e-04
30  cyl     am -0.52260705 2.151207e-03
31 disp     am -0.59122704 3.662114e-04
32   hp     am -0.24320426 1.798309e-01
33 drat     am  0.71271113 4.726790e-06
34   wt     am -0.69249526 1.125440e-05
35 qsec     am -0.22986086 2.056621e-01
36   vs     am  0.16834512 3.570439e-01
37  mpg   gear  0.48028476 5.400948e-03
38  cyl   gear -0.49268660 4.173297e-03
39 disp   gear -0.55556920 9.635921e-04
40   hp   gear -0.12570426 4.930119e-01
41 drat   gear  0.69961013 8.360110e-06
42   wt   gear -0.58328700 4.586601e-04
43 qsec   gear -0.21268223 2.425344e-01
44   vs   gear  0.20602335 2.579439e-01
45   am   gear  0.79405876 5.834043e-08
46  mpg   carb -0.55092507 1.084446e-03
47  cyl   carb  0.52698829 1.942340e-03
48 disp   carb  0.39497686 2.526789e-02
49   hp   carb  0.74981247 7.827810e-07
50 drat   carb -0.09078980 6.211834e-01
51   wt   carb  0.42760594 1.463861e-02
52 qsec   carb -0.65624923 4.536949e-05
53   vs   carb -0.56960714 6.670496e-04
54   am   carb  0.05753435 7.544526e-01
55 gear   carb  0.27407284 1.290291e-01

```

corrplot
```
  install.packages("corrplot")
  library(corrplot)
  cor<-cor(data,method = "pearson")
  corrplot(cor,order = "hclust",method="shade",addrect = 2,tl.cex=0.25,tl.col = "blue")
  write.table(cor,file="Pt_RLW_FPKM_cor.txt",sep=" ",row.names = TRUE,col.names = TRUE,quote=F)
```
scatter plots
```
library(PerformanceAnalytics)
chart.Correlation(mtcars,histogram = TRUE,pch=19)
```
![image](https://upload-images.jianshu.io/upload_images/9218360-772f8e43d7d95533.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1000)