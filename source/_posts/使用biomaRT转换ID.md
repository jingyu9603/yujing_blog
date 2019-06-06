---
title: 使用biomaRT转换ID
tag: biomaRT
categories: R
---

```R
library(biomaRt)
#显示包含的数据库及其版本
```
<!--more-->
显示包含的数据库及其版本,
```
> listMarts()
               biomart               version
1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 94
2   ENSEMBL_MART_MOUSE      Mouse strains 94
3     ENSEMBL_MART_SNP  Ensembl Variation 94
4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 94
```
本次使用选用ENSEMBL_MART_ENSEMBL这个数据库

```
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
```
显示该数据库包含的子数据库

```
head(listDatasets(ensembl))
```

```
47            ggallus_gene_ensembl
48           ggorilla_gene_ensembl
49            gmorhua_gene_ensembl
50           hburtoni_gene_ensembl
51             hcomes_gene_ensembl
52            hfemale_gene_ensembl
53              hmale_gene_ensembl
54           hsapiens_gene_ensembl
```
选择查询的数据库及其相关数据集

```
ensembl<-useDataset("hsapiens_gene_ensembl",mart=ensembl)
```
查询filter函数包含的属性，这里filter函数代表的输入（即已知信息）的属性

```
filters<-listFilters(ensembl)
head(filters)
```
显示attributes函数的属性，这里attributes函数要选择需要查询的属性

```
attributes<-listAttributes(ensembl)
head(attributes)
```

```
                           name                  description         page
1               ensembl_gene_id               Gene stable ID feature_page
2       ensembl_gene_id_version       Gene stable ID version feature_page
3         ensembl_transcript_id         Transcript stable ID feature_page
4 ensembl_transcript_id_version Transcript stable ID version feature_page
5            ensembl_peptide_id            Protein stable ID feature_page
6    ensembl_peptide_id_version    Protein stable ID version feature_page
```

```R
rownames_old=rownames(Pt_union_20)
##去掉id后面的小数点，也就是版本号
ens=c()
for(i in rownames_old){
  ens[i]=strsplit(i,"\\.")[[1]][1]
}
names(ens)<-NULL
#利用getBM()函数进行ID转换
gene_symbol=getBM(attributes = c("ensembl_transcript_id","hgnc_symbol"),filters = "ensembl_transcript_id",values = ens,mart=ensembl)
```

```R
##将ID转换过程中，在数据库中没有找到对应ID的空白值转换为“none字符”
for(j in which(gene_symbol[,2]=="")){
  gene_symbol[j,2]="none"
}
```

```R
#将ID与genesymbol连起来，形成新的ID
rownames_new=c()
for(i in rownames_old){
  trs_id=strsplit(i,"\\.")[[1]][1]
  print(trs_id)
  index=match(trs_id,gene_symbol$ensembl_transcript_id)
  print(index)
  rownames_new[i]<-paste(i,gene_symbol[index,2],sep=",")
}
names(rownames_new)<-NULL
```
当需要转换的id中既有gene_id又有transcrption_id时，可以采用以下的方法：

```
rownames_old=rownames(Pt_union_top100)
g_id=c()
t_id=c()
for(i in rownames_old){
  id=strsplit(i,"\\.")[[1]][1]
  if(grepl("G",id)){
    g_id[i]=id
  }
  else{
    t_id[i]=id
  }
}
names(g_id)<-NULL
names(t_id)<-NULL

gene_symbol_t=getBM(attributes = c("ensembl_transcript_id","hgnc_symbol"),filters = "ensembl_transcript_id",values = t_id,mart=ensembl)
gene_symbol_g=getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),filters="ensembl_gene_id",values = g_id,mart=ensembl)

for(j in which(gene_symbol_t[,2]=="")){
  gene_symbol_t[j,2]="none"
}
for(j in which(gene_symbol_g[,2]=="")){
  gene_symbol_g[j,2]="none"
}
rownames_new=c()
for(i in rownames_old){
  trs_id=strsplit(i,"\\.")[[1]][1]
  print(trs_id)
  if(grepl("G",trs_id)){
    index=match(trs_id,gene_symbol_g$ensembl_gene_id)
    rownames_new[i]<-paste(i,gene_symbol_g[index,2],sep=";")
  }
  else{
    index=match(trs_id,gene_symbol_t$ensembl_transcript_id)
    rownames_new[i]<-paste(i,gene_symbol_t[index,2],sep=";")
  }
}
names(rownames_new)<-NULL
rownames(Pt_union_top100)<-rownames_new
```

