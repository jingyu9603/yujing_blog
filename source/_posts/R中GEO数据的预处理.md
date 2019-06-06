---
title: R中GEO数据的预处理
tag: GEO
categories: R
---

获得芯片探针与gene的对应关系的三种方法：  
- 在基因芯片的厂商的官网直接下载
- 从NCBI里面下载文件来解析
- 直接使用Bioconductor包  
<!--more-->
以下是直接使用Biocondutor包的流程：

- 获取平台与Bioconductor包之间的对应关系

```
library(BiocManager)
BiocManager::install("GEOmetadb")
install.packages("RSQLite")
devtools::install_github("ggrothendieck/sqldf")
library(GEOmetadb)
library(RSQLite)
library(sqldf)
getSQLiteFile()
con<-dbConnect(RSQLite::SQLite(),"GEOmetadb.sqlite")
gplToBioc<-dbGetQuery(con,'select gpl,bioc_package,title from gpl where bioc_package is not null')
```
- 下载所需要的包
```
for(i in 1:nrow(gplToBioc)){
  platform=gplToBioc[i,1]
  platformDB=paste(platform,".db",sep = "")
  if(platformDB %in% rownames(installed.packages()) ==FALSE){
    BiocManager::install(platformDB)
  }
}
```
- 批量获得芯片探针与gene的对应关系
```
library(hgu133a.db)
library(GEOquery)
#下载GEO数据集
GEO<-"GSE24673"
if(!file.exists(GEO)){
  gset<-getGEO(GEO,destdir = ".",getGPL = F,AnnotGPL = F)
  save(gset,file=GEO)
}
data<-gset[[1]]
data<-exprs(data)
ids<-toTable(hgu133aSYMBOL)
data<-data[ids$probe_id,]
ids$median<-apply(data,1,median)
ids<-ids[order(ids$symbol,ids$median,decreasing = T),]
ids<-ids[!duplicated(ids$symbol),]
data<-data[ids$probe_id,]
rownames(data)<-ids$symbol
```