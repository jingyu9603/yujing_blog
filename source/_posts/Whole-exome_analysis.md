---
title: Whole-exome analysis
tag: exome
categories: Bioinformatic_project
---
#####  软件准备
在相关软件官网下载最新版本的软件  
- bwa-0.7.17.tar.bz2   
- gatk-4.1.2.0.zip  
- picard-2.20.0-0.tar.bz2   
- samtools-1.9.tar.bz2 
<!--more-->
##### 原始数据质量控制
分析数据质量的几个方面：  
- read各个位置的碱基质量分布
- 碱基的总体质量分布
- read各个位置上碱基分布比例看，目的是为了分析碱基的分离程度
- GC含量分布
- read各个位置的N含量
- read是否还包含测序的接头序列
- read重复率，这个是实验的扩增所引入的  

**需要的软件**
- FastQC，一个java程序，可以很好的帮助我们理解测序数据的质量情况，唯一不足的就是图太丑了
- MultiQC,一个基于Python的小工具，能够将测序数据的多个QC结果整合成一个HTLM网页交互式报告，同时也能导出pdf
- Trimmomatic（cutadapt、sickle、seqtk）等可用于切除接头序列和read的低质量序列

##### 全外显子组分析（WES）
需要的工具：  
- bwa
- samtools
- bgzip
- GATK 4.0
- GATK 3.x
- Picard

######  数据预处理
**序列比对**
- 为参考基因组构建索引
```
bwa index hg38.fa
```
- 将read比对到参考基因组
```
bwa mem -t 4 -R '@RG\tID:foo_lane\tPL:illumina\tLB:library\tSM:sample_name' /path/to/human.fasta read_1.fq.gz read_2.fq.gz > sample_name.sam
```
-t，线程数，我们在这里使用4个线程  ；  
-R 接的是Read Group的字符串信息，这是一个非常重要的信息，以@RG开头，它是用来将比对的read进行分组的。这个信息对于我们后续对比对数据进行错误率分析和Mark duplicate时非常重要。在Read Group中，有如下几个信息非常重要：  
- 1）ID，这是Read Group的分组ID，一般设置为测序的lane ID（不同lane之间的测序过程认为是独立的），下机数据中我们都能看到这个信息的，一般都是包含在fastq的文件名中；
- 2）PL，指的是所用的测序平台，这个信息不要随便写！特别是当我们需要使用GATK进行后续分析的时候，更是如此！这是一个很多新手都容易忽视的一个地方，在GATK中，PL只允许被设置为：ILLUMINA,SLX,SOLEXA,SOLID,454,LS454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS或UNKNOWN这几个信息。
- 3）SM，样本ID，同样非常重要，有时候我们测序的数据比较多的时候，那么可能会分成多个不同的lane分布测出来，这个时候SM名字就是可以用于区分这些样本；
- 4）LB，测序文库的名字，这个重要性稍微低一些，主要也是为了协助区分不同的group而存在。文库名字一般可以在下机的fq文件名中找到，如果上面的lane ID足够用于区分的话，也可以不用设置LB。  

为了方便后续的分析，一般在输出文件时，将其转换为bam文件
```
bwa mem -t 4 -R '@RG\tID:foo_lane\tPL:illumina\tLB:library\tSM:sample_name' /path/to/human.fasta read_1.fq.gz read_2.fq.gz | samtools view -S -b - > sample_name.bam
```
**排序（sort）**  
第一步的比对是按照FASTQ文件的顺序把read逐一定位到参考基因组上之后，随即就输出了，它不会也不可能在这一步里面能够自动识别比对位置的先后位置重排比对结果。因此，比对后得到的结果文件中，每一条记录之间位置的先后顺序是乱的，我们后续去重复等步骤都需要在比对记录按照顺序从小到大排序下来才能进行，所以这才是需要进行排序的原因。
```
time samtools sort -@ 4 -m 4G -O bam -o sample_name.sorted.bam sample_name.bam
```
-@，用于设定排序时的线程数，我们设为4；  
-m，限制排序时最大的内存消耗，这里设为4GB；  
-O 指定输出为bam格式；  
-o 是输出文件的名字，

**标记重复序列（删除重复序列）**  
我们使用Picard来完成这个事情
```
##标记重复序列
java -jar picard.jar MarkDuplicates \ 
  I=sample_name.sorted.bam \
  O=sample_name.sorted.markdup.bam \
  M=sample_name.markdup_metrics.txt
```
这里只把重复序列在输出的新结果中标记出来，但不删除。如果我们非要把这些序列完全删除的话可以这样做：
```
java -jar picard.jar MarkDuplicates \ 
  REMOVE_DUPLICATES=true \
  I=sample_name.sorted.bam \
  O=sample_name.sorted.markdup.bam \
  M=sample_name.markdup_metrics.txt
```
在后面的分析中，在`Somatic SNVs+Indels`的过程中使用标记的bam的文件，并且使用GATK4.0的主流pipline；在`Detect MSI（微卫星不稳定性）`的过程中使用GATK3.x的过程

###### Somatic SNVs+Indels  
**重新矫正碱基质量值（BQSR）**
分为两步:
```{R}
 gatk BaseRecalibrator \
   -I my_reads.bam \
   -R reference.fasta \
   --known-sites sites_of_variation.vcf \
   --known-sites another/optional/setOfSitesToMask.vcf \
   -O recal_data.table
## 用到的已知变异位点
common_all_20180418_dbsnp_151_hg38.vcf
Homo_sapiens_assembly38.known_indels.vcf
Mills_and_1000G_gold_standard.indels.hg38.vcf
下载地址（GATK的bundle以及1000genomic官网）：ftp://ftp.broadinstitute.org/bundle/hg38/
https://ftp.ncbi.nlm.nih.gov/snp/organisms/
```

```
 gatk ApplyBQSR \
   -R reference.fasta \
   -I input.bam \
   --bqsr-recal-file recalibration.table \
   -O output.bam
```

**查看测序深度**
```
samtools depth  bamfile  |  awk '{sum+=$3} END { print "Average = "sum/NR}'
```

**使用qualimap来查看比对好的bam文件的质量**  
在[官网](http://qualimap.bioinfo.cipf.es/)下载使用qualimap
```
wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v.2.12.zip
unzip qualimap_v.2.12.zip
cd qualimap_v.2.12
./qualimap -h
# 使用qualimap时 一般需要添加内存参数,对于外显子数据，需要添加关于外显子区域的bed文件
qualimap --java-mem-size=8G bamqc -bam in.bam -gff exome.bed
```
###### 使用旧版的mutects来call somatic mutation（现在GATK网站已经不推荐了）
** CreateSomaticPanelOfNormals **
panel of normals (PoN) containing germline and artifactual sites for use with Mutect2.  

Step 1. Run Mutect2 in tumor-only mode for each normal sample.
```
~/biosoft/GATK4.0/gatk-4.0.5.1/gatk Mutect2 \
-R ~/reference/genome/gatk_hg38/Homo_sapiens_assembly38.fasta \
-I HG00190.bam \
-tumor HG00190 \
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
-L interval_list \
--max-map-distance 0 \ ## 防止在产生的vcf文件中出现MNPS现象
-O 3_HG00190.vcf.gz
```
其中的interval_list来自于安捷伦公司的官网，使用的具体bed文件为：S07604514_Padded.bed   
或者是可以通过CCDS文件自己制作interval_list   
- 首先在CCDS官网（ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human）下载文件
- 将txt文件转换为bed形式
```
cat CCDS.20160908.txt |grep -w -v Withdrawn|perl -alne '{/\[(.*?)\]/;next unless $1;$gene=$F[2];$exons=$1;$exons=~s/\s//g;$exons=~s/-/\t/g;print "$F[0]\t$_\t$gene" foreach split/,/,$exons;}'|sort -u |bedtools sort -i >exon_probe.hg38.gene.bed
```
- 制作interval_list
```
cat exon_probe.hg38.gene.bed | awk '{print "chr"$0}' >hg38.chr.bed

picard BedToIntervalList \
       I=hg38.chr.bed \
       O=hg38.interval_list\
       SD=/data/reference/GRCh38/GRCh38.primary_assembly.genome.dict
       
gatk4 PreprocessIntervals 
      -L hg38.interval_list  
      --sequence-dictionary /data/reference/GRCh38/GRCh38.primary_assembly.genome.dict 
      --reference /data/reference/GRCh38/GRCh38.primary_assembly.genome.fa 
      --padding 100 
      --bin-length 0 
      --interval-merging-rule OVERLAPPING_ONLY
      --output hg38.Padded.interval_list
```

Step 2. Create a GenomicsDB from the normal Mutect2 calls.
由于在第一步产生的vcf文件中有多行对应一个位点，所以在进行Step2之前要先将这些同一位点的行合并起来
```
bcftools norm \
    -m +any --do-not-normalize \
    delta_af-only-gnomad_Hg19toGRCh38.vcf.gz \
    -Oz -o zeta_af-only-gnomad_Hg19toGRCh38.vcf.gz
## 为norm之后的新vcf文件创建索引
gatk IndexFeatureFile \
     -F cohort.vcf.gz
```
```R
## step.2
 gatk GenomicsDBImport -R reference.fasta -L intervals.interval_list \
       --genomicsdb-workspace-path pon_db \
       -V normal1.vcf.gz \
       -V normal2.vcf.gz \
       -V normal3.vcf.gz
```
Step 3. Combine the normal calls using CreateSomaticPanelOfNormals.
```
gatk CreateSomaticPanelOfNormals -R reference.fasta -V gendb://pon_db -O pon.vcf.gz
```
** 对tumor和matched normal进行calls somatic variants **
```
for sample in *C*.bam;  
do base=${sample%-*};  
   gatk4 Mutect2   
   -R /data/reference/GRCh38/GRCh38.primary_assembly.genome.fa   
   -I "$base"-N.sorted.markdup.recal.bam  
   -I $sample   
   -tumor "$base"-C   
   -normal "$base"-N   
   -pon ./Mutect2_interval/pon.vcf.gz  
   --germline-resource ../../reference/somatic-hg38_af-only-gnomad.hg38.vcf.gz  
   --af-of-alleles-not-in-resource 0.0000025  
   -disable-read-filter MateOnSameContigOrNoMappedMateReadFilter  
   -L ../../reference/S07604514_hs_hg38/S07604514_hs_hg38/S07604514_Padded.bed   
   -O ./Mutect2_interval/Mutect2_output/"$base"_somatic_m2.vcf.gz   
   -bamout ./Mutect2_interval/Mutect2_output/"$base"_tumor_normal_m2.bam;done
```
** 使用`GetPileupSummaries`工具来计算肿瘤样本的污染情况 ** 
先使用Selectvariants工具生成包含commom biallelic variants的.vcf文件
```R
 
gatk SelectVariants \
-R ~/reference/genome/gatk_hg38/Homo_sapiens_assembly38.fasta \
-V af-only-gnomad.hg38.vcf.gz \
--select-type-to-include SNP \
--restrict-alleles-to BIALLELIC \
-O af-only-gnomad.hg38.SNP_biallelic.vcf.gz
```
接着用GetPileupSummaries计算resource site位点上的count数，这样并不是计算所有的位点，而是用AF参数过滤后的（大于0.01并小于0.2），同时也可以用 -L 参数指定区域，分别对Normal和Tumor样本计算
```
~/biosoft/GATK4.0/gatk-4.0.5.1/gatk GetPileupSummaries \
-I Normal_blood.allready.bam \
-L ../exon_150bp.list \
-V ~/annotation/GATK/resources/bundle/hg38/af-only-gnomad.hg38.SNP_biallelic.vcf.gz \
-O Normal_blood.pileups.table
```
总之最后通过计算出的污染比例（XXX.calculatecontamination.table文件中）来过滤掉somatic variant中的一些可能是污染导致的假阳性突变
```
gatk CalculateContamination \
-I Primary-IIIG.pileups.table \
-matched Normal_blood.pileups.table \
-O Primary-IIIG.calculatecontamination.table
```

###### 使用Muect2新流程进行calling somatic mutation
首先需要创建PON，步骤与旧版的相同，接下买的操作如下
```
## First, run Mutect2 with the --f1r2-tar-gz argument. This creates an output with raw data
if [ 1 -eq 1 ];
then
	workdir=/media/data1/YJ/CRC/tag_duplicates/tag-dup_rec/;
outdir=/media/data1/YJ/CRC/tag_duplicates/tag-dup_rec/Mutect2_interval/Mutect2_new_withnewPON_matchedmodel/;
refdir=/media/data1/YJ/CRC/reference/;
cd $workdir;
for sample in *C*.bam;
do
	base=${sample%%-*};
	gatk4 Mutect2 -R /data/reference/GRCh38/GRCh38.primary_assembly.genome.fa \
		-L /media/data1/YJ/CRC/reference/S07604514_hs_hg38/S07604514_hs_hg38/S07604514_Padded.bed \
		-I $sample \
		-I "$base"-N.sorted.markdup.recal.bam \
		-tumor "$base"-C \
		-normal "$base"-N \
		-germline-resource "$refdir"somatic-hg38_af-only-gnomad.hg38.vcf.gz \
		-pon "$workdir"Mutect2_interval/pon_withgnomAD.vcf.gz \
		--f1r2-tar-gz "$outdir""$base".f1r2.tar.gz \
		-O "$outdir""$base".unfiltered.vcf;
done;
fi

## Next, pass this raw data to LearnReadOrientationModel:
if [ 1 -eq 1 ] ;
then
	workdir=/media/data1/YJ/CRC/tag_duplicates/tag-dup_rec/Mutect2_interval/Mutect2_new_withnewPON_matchedmodel;
	cd $workdir;
	for sample in *.tar.gz;
	do
		base=${sample%%.*};
		gatk4 LearnReadOrientationModel -I $sample -O "$base"_read-orientation-model.tar.gz;
	done;
fi

## Finally, pass the learned read orientation model to FilterMutectCallswith the -ob-priors argument:
if [ 1 -eq 1 ];
then
	workdir=/media/data1/YJ/CRC/tag_duplicates/tag-dup_rec/Mutect2_interval/Mutect2_new_withnewPON_matchedmodel;
	refdir=/data/reference/GRCh38;
	cd $workdir
	for sample in *.unfiltered.vcf;
	do
		base=${sample%%.*};
		gatk4 FilterMutectCalls -V $sample \
			--ob-priors "$base"_read-orientation-model.tar.gz \
			-R "$refdir"/GRCh38.primary_assembly.genome.fa \
			-O "$base"_filtered.vcf
	done;
fi
```
最终获得利用`FilterMutectCalls`过滤过的filtered.vcf文件然后得到通过过滤的vcf文件
```R
##filter the .filtered.vcf and annotation it with gene_based model
if [ 1 -eq 1 ];
then
	for sample in *_filtered.vcf;
	do
	
		base=${sample%%_*};
		grep -v "#" $sample | awk '$7=="PASS"' > "$base"_filtered_PASS.cvf;
	done;
fi
```

可以通过`VariantAnnotator`对vcf文件进行dbsnp注释
```
##Annotation vcf with dbsnp                                                                                                   
if [ 1 -eq 0 ];                                                                                                                
then                                                                                                   
        for sample in *_filtered.vcf;                                                               
        do
                base=${sample%%_*};  
                gatk4 VariantAnnotator -R /data/reference/GRCh38/GRCh38.primary_assembly.genome.fa \ 
                        -I ../../"$base"-C.sorted.markdup.recal.bam \  
                        -V $sample \        
                        --output "$base"_filtered_dbsnp.vcf \ 
                        -A Coverage \           
                        --dbsnp /media/data1/YJ/CRC/reference/common_all_20180418_dbsnp_151_hg38.vcf ;
        done                   
fi 
```

注释完的vcf文件中多了一列ID的信息，里面有snp号，可以通过[SNPedia网站](https://www.snpedia.com/index.php/SNPedia)查看snp的具体信息

###### 使用不同的方法进行annotation  

**ANNOVAR**
在ANNOVAR官网下载最新版本的软件
```R
## 下载适合的数据库(下载的数据储存在GRCH38文件夹中)
 annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene GRCH38/
## 利用 ANNOVAR进行基于gene的annotation
if [ 1 -eq 1 ];
then
	sofadir=/media/data1/YJ/CRC/softwave/annovar;
	for sample in *filtered_PASS.vcf;
	do
		base=${sample%%.*};
		
		$sofadir/convert2annovar.pl -format vcf4 $sample -outfile "$base".avinput;
		$sofadir/annotate_variation.pl -out "$base" \
			-build hg38 \
			"$base".avinput \
			--geneanno -dbtype refGene \
			$sofadir/GRCH38/
	done;
fi

```
结果会得到两个文件：
- 在外显子位点的注释：s0020487998_filtered_PASS.exonic_variant_function
![image](https://note.youdao.com/yws/public/resource/bed499f2133fd1cc1e1ced485447f8c1/xmlnote/WEBRESOURCE22c444411c67fe341deac1f544ccb7ae/6557)
- 在所有位点的注释：
s0020487998_filtered_PASS.variant_function
![image](https://note.youdao.com/yws/public/resource/bed499f2133fd1cc1e1ced485447f8c1/xmlnote/WEBRESOURCEe3710742911832dc4bc8e26297d630c7/6567)
XX.variant_function文件一般关心前两列，后面几列均是变异位点的一些信息
第一列是变异所在基因组的位置，如：exonic、splicing、ncRNA等；其优先级：exonic = splicing > ncRNA> > UTR5/UTR3 > intron > upstream/downstream > intergenic
第二列信息则给突变位点所在的基因名称（如果突变在exonic/intronic/ncRNA），或者给出临近基因的名称
**ClinEff**  
在ClinEff的[官网](http://www.dnaminer.com/clineff.html)下载最新版本的软件已经相关的database
```
tar -xvf clinEff.1.0h.tgz
tar -xvf clinEff_db38_v1.0.tgz
```
**SnpEff**  
在[SnpEff的官网](http://snpeff.sourceforge.net/download.html)下载最新版本的软件
```R
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
#查看homo_sapiens相关的数据库
java -jar snpEff.jar databases | grep "Home_sapiens" | cut -f 1,2

GRCh37.75       Homo_sapiens
GRCh38.86       Homo_sapiens
hg19            Homo_sapiens (UCSC)
hg19kg          Homo_sapiens (UCSC KnownGenes)
hg38            Homo_sapiens (UCSC)
hg38kg          Homo_sapiens (UCSC KnownGenes)
testHg19ChrM    Homo_sapiens (UCSC)

#由于我们之前是使用的GRCh38的基因组，所以我们就下载GRCh38.86
java -jar snpEff.jar download GRCh38.86
#注释
java -Xmx4g -jar /media/data1/YJ/CRC/softwave/snpEff/snpEff.jar GRCh38.86 s0020347988_filtered_PASS.vcf > s0020347988_filtered_PASS.snpEff.vcf
```
输出的vcf文件，在输入vcf文件的基础上添加了一些tag：ANN、LOF、NMD  
因此我可以将这个vcf格式文件稍微处理下，保留原来的vcf文件的前5列，再加上ANNtag形成一个新文件来查看
```
for sample in *PASS*;
do 
    base=${sample%%_*};
    perl -alne 'next if $_ =~ /^#/;$F[7] =~ /(ANN=\S+)/;print "$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$1"' $sample > "$base".anntag.vcf
done
```
vcf文件的内容如下：
```
29275630        .       C       T       .       ANN=T|stop_gained|HIGH|PTPRU|ENSG00000060656|transcript|ENST00000345512.7|protein_coding|8/31|c.1327C>T|p.Gln443*|1456/4470|1327/4341|443/1446||, 
                                                    T|stop_gained|HIGH|PTPRU|ENSG00000060656|transcript|ENST00000373779.7|protein_coding|8/30|c.1327C>T|p.Gln443*|1456/5579|1327/4311|443/1436||,
                                                    T|stop_gained|HIGH|PTPRU|ENSG00000060656|transcript|ENST00000428026.6|protein_coding|8/30|c.1327C>T|p.Gln443*|1423/5550|1327/4302|443/1433||,
                                                    T|stop_gained|HIGH|PTPRU|ENSG00000060656|transcript|ENST00000460170.2|protein_coding|8/31|c.1327C>T|p.Gln443*|1333/5337|1327/4323|443/1440||,
                                                    T|upstream_gene_variant|MODIFIER|PTPRU|ENSG00000060656|transcript|ENST00000415600.6|processed_transcript||n.-3921C>T|||||3921|,
                                                    T|non_coding_transcript_exon_variant|MODIFIER|PTPRU|ENSG00000060656|transcript|ENST00000527027.1|processed_transcript|3/3|n.564C>T||||||;
                                                    LOF=(PTPRU|ENSG00000060656|11|0.36);
                                                    NMD=(PTPRU|ENSG00000060656|11|0.36) 
```
往往一个突变条目对应着多个annotation报告，主要有一下几种原因：
- 一个突变位点可以影响多个基因，比如一个突变可能是一个基因的DOWNSTREAM同时可能是另一个基因的UPSTREAM
- 在复杂的基因组中，一个基因可能对应着多个转录组，snpeff会指出受该突变影响的每一个转录本
- 一条vcf记录可能包含不止一种突变
```
#CHROM  POS      ID    REF  ALT    QUAL FILTER    INFO
1       889455   .     G    A,T    .    .         ANN=A|stop_gained|HIGH|NOC2L|ENSG00000188976|transcript|ENST00000327044|protein_coding|7/19|c.706C>T|p.Gln236*|756/2790|706/2250|236/749||
                                                     ,T|missense_variant|MODERATE|NOC2L|ENSG00000188976|transcript|ENST00000327044|protein_coding|7/19|c.706C>A|p.Gln236Lys|756/2790|706/2250|236/749||
                                                     ,A|downstream_gene_variant|MODIFIER|NOC2L|ENSG00000188976|transcript|ENST00000487214|processed_transcript||n.*865C>T|||||351|
                                                     ,T|downstream_gene_variant|MODIFIER|NOC2L|ENSG00000188976|transcript|ENST00000487214|processed_transcript||n.*865C>A|||||351|
                                                     ,A|downstream_gene_variant|MODIFIER|NOC2L|ENSG00000188976|transcript|ENST00000469563|retained_intron||n.*878C>T|||||4171|
                                                     ,T|downstream_gene_variant|MODIFIER|NOC2L|ENSG00000188976|transcript|ENST00000469563|retained_intron||n.*878C>A|||||4171|
                                                     ,A|non_coding_exon_variant|MODIFIER|NOC2L|ENSG00000188976|transcript|ENST00000477976|retained_intron|5/17|n.2153C>T||||||
                                                     ,T|non_coding_exon_variant|MODIFIER|NOC2L|ENSG00000188976|transcript|ENST00000477976|retained_intron|5/17|n.2153C>A||||||;LOF=(NOC2L|ENSG00000188976|6|0.17);NMD=(NOC2L|ENSG00000188976|6|0.17)
```

当存在多种effect时，effect排序的依据：
- 根据Putative impact：具有更高等级的putative impact将会排在前面
- Effect type：有害的effect会排在前面
- Canonical trancript before non-canonical.
- Marker genomic coordinates (e.g. genes starting before first).　　


`ANN`tag将注释信息以“|”分割，每个filed有其对应的信息：
- Allele（or ALT）：突变的碱基
- Annotation：表示突变的类型，当有多种类型存在时，用@连接，可以在 [The Sequence ontology官网](http://sequenceontology.org/browser/current_svn/term/SO:0001792)查询每一种突变类型的具体含义
- Putative_impact：对突变的影响进行的预测，有4个程度
1.HIGH：突变对蛋白的影响很大，可能导致蛋白质截短、功能丧失或引发无意义的介导衰退。eg.stop_gained, frameshift_variant
2.MODERATE：一种可能改变蛋白质有效性的无破坏性突变  eg.missense_variant, inframe_deletion  
3.LOW：几乎是没害的，不会改变蛋白质的性能，eg.synonymous_variant  
4.MODIFIER:通常是影响非编码基因的突变或非编码区的突变，对其很难预测它的影响。eg.exon_variant, downstream_gene_variant
- Gene ID：使用ENSEMBL id
- Feature type：表示突变所在区域的类型，比如transcript, motif, miRNA等
- Feature ID ：表示Feature type对应的id
- Transcript biotype ：关于transcript的类型，根据ENSEMBL biotypes来定义，
- Rank/total：1/1表示Exon or Intron rank / total number of exons or introns，前面的1表示这个突变是在第1个exon上（因为annatation已经给出了这个是突变是在exon上），后面的11表示这个突变所在的transcript总共有1个exon
- HGVS.c ：n.321T>C表示Variant using HGVS notation (DNA level)
- HGVS.p：如果突变发生在ｃｏｄｉｎｇ区域，根据HGVS标记法在蛋白质水平表示突变类型。p.Gln443*表示443位的Gln氨基酸突变成了终止密码子
- cDNA_position/cDNA_len:cDNA的位置以及，cDNA的长度
- CDS_position / CDS_len: Position and number of coding bases (one based includes START and STOP codons).
- Protein_position / Protein_len: Position and number of AA (one based, including START, but not STOP).

至于LOF和NMD标签这是表示：Loss of funcation(LOF) and nonsense-mediated   decay(NMD) 预测  ，分析相关的effect是否可以产生LOF或者是NMD影响  

NMD：无义介导的mRNA降解作为真核细胞中重要RNA监控机制，识别并降解开放阅读框中含有提前终止密码子（premature termination codon，PTC）的mRNA，以避免因截短的蛋白产物积累对细胞造成毒害. NMD还调控正常生理基因的表达，暗示其在真核细胞中扮演重要角色.  
 
 LOF和NMD标签的形式如下：
 - Gene：gene的名称
 - ID：gene id，通常是ENSEMBL id
 - Num_transcripts：gene中的转录本的数量
 - percent_affected：被改突变多影响的转录本的比例
 

###### 使用VarScan2 Somatic来call somatic（没有写完）

###### 环境搭建
1 创建文件夹
```
mkdir -p input reference scripts temp
```
2 安装conda并使用conda创建env
```R
# 添加信号源
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
# 创建env
conda create --name altwes gatk R samtools trimmomatic picard bam-readcount
varscan
# 激活env
conda activate altwes
```
3 安装GATK3.8
```
WES=/
cd $WES/scripts
wget -r -np -nd \
-O gatk3.tar.bz2 \
’https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-0-ge9d806836’
bunzip2 gatk3.tar.bz2
tar -xf gatk3.tar
```

##### Detect MSI  

###### 数据预处理

**局部重比对**  
接下来是局部区域重比对，通常也叫Indel局部区域重比对。有时在进行这一步骤之前还有一个merge的操作，将同个样本的所有比对结果合并成唯一一个大的BAM文件，merge的例子如下：
```
samtools merge <out.bam> <in1.bam> [<in2.bam>...<inN.bam>]
```
局部重比对的目的是将BWA比对过程中所发现有潜在序列插入或者序列删除（insertion和deletion，简称Indel）的区域进行重新校正。这个过程往往还会把一些已知的Indel区域一并作为重比对的区域  
GATK4.0中没有分离出单独的局部比对脚本。只能使用3.X版本中的脚本,这里包含了两个步骤：
- 第一步，RealignerTargetCreator ，目的是定位出所有需要进行序列重比对的目标区域
- 第二步，IndelRealigner，对所有在第一步中找到的目标区域运用算法进行序列重比对，最后得到捋顺了的新结果。
```{Linux}
java -jar /path/to/GenomeAnalysisTK.jar \
 -T RealignerTargetCreator \
 -R /path/to/human.fasta \
 -I sample_name.sorted.markdup.bam \
 -known /path/to/gatk/bundle/1000G_phase1.indels.b37.vcf \
 -known /path/to/gatk/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf \
 -o sample_name.IndelRealigner.intervals

java -jar /path/to/GenomeAnalysisTK.jar \
 -T IndelRealigner \
 -R /path/to/human.fasta \
 -I sample_name.sorted.markdup.bam \
 -known /path/to/gatk/bundle/1000G_phase1.indels.b37.vcf \
 -known /path/to/gatk/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf \
 -o sample_name.sorted.markdup.realign.bam \
 --targetIntervals sample_name.IndelRealigner.intervals
 
 
```
这一步中我们实际使用到的known indel集为：
 Homo_sapiens_assembly38.known_indels.vcf
 Mills_and_1000G_gold_standard.indels.hg38.vcf
 这些文件可以方便的在GATK bundle里面下载（ftp://ftp.broadinstitute.org/bundle/hg38/）  
 
 后面的变异检测使用GATK，而且使用GATK的HaplotypeCaller模块或者是Mutect2模块时，当这个时候才可以减少这个Indel局部重比对的步骤。 
 
 **重新矫正碱基质量值**
 分为两步:
```{R}
 gatk BaseRecalibrator \
   -I my_reads.bam \
   -R reference.fasta \
   --known-sites sites_of_variation.vcf \
   --known-sites another/optional/setOfSitesToMask.vcf \
   -O recal_data.table
## 用到的已知编译位点
common_all_20180418_dbsnp_151_hg38.vcf
Homo_sapiens_assembly38.known_indels.vcf
Mills_and_1000G_gold_standard.indels.hg38.vcf
下载地址（GATK的bundle以及1000genomic官网）：ftp://ftp.broadinstitute.org/bundle/hg38/
https://ftp.ncbi.nlm.nih.gov/snp/organisms/
```

```
 gatk ApplyBQSR \
   -R reference.fasta \
   -I input.bam \
   --bqsr-recal-file recalibration.table \
   -O output.bam
```
###### MSI鉴定  
通过以上步骤我么就获得了干净的bam文件可以用于下一步的MSI鉴定分析。  
在此次鉴定中一共使用到了三种MSI鉴定软件用于鉴定我们的tumer-normal paired样本

** MSIsensor **  
MSIsensor通过分别判断每个微卫星位点的稳定性，然后以不稳定微卫星位点的比例作为MSI得分。MSIsensor需要基于配对的肿瘤-正常样本进行MSI的判定。
 首先, 对于在肿瘤和正常样本中测序深度均大于等于20的微卫星位点, 计算其等位基因的分布信息; 其次, 通过卡方检验比较肿瘤和正常样本的相同微卫星位点的等位基因分布, 若显著不同, 则认为该微卫星位点不稳定; 最后统计不稳定位点的比例, 若该比例超过阈值, 则判定为MSI-H,这个阈值一般是3.5%
 ```R
 #重参考基因组中获取微卫星位点
 msisensor scan -d reference.fa -o microsatellites.list
 #MSI扫描
 msisensor msi -d microsatellites.list -n normal.bam -t tumor.bam -e bed.file -o output.prefix
 ```
 结果显示3个样本都是MSS  
** MANTIS **
`MANTIS`软件使用Python编辑，最新的版本2.7.1,但是Python3.0版本也是支持的
下载地址：https://github.com/OSU-SRLab/MANTIS  
类似于MSIsensor, MANTIS也获得了肿瘤-正常配对样本在每个微卫星位点的等位基因分布信息; 与MSIsensor不同的是, 对于每个微卫星位点, MANTIS把上述两组数据看作两个向量, 定义这两个向量的 L1范数为样本中该位点的稳定程度, 对所有位点的L1范数求平均值即为样本的MSI得分
```R
 #从参考基因组中获得微卫星位点
 ./RepeatFinder –i /path/to/genome.fasta -o /path/to/loci.bed
 #对上一步获得的微卫星位点根据exon测序实验室中的target位置信息进行过滤
 bedtools intersect -a MANTIS.loc.bed -b hg38_Regions.bed > MANTIS_filter.bed
 #进行MSI鉴定
 python mantis.py --bedfile /path/to/loci.bed --genome /path/to/genome.fasta -n /path/to/normal.bam -t /path/to/tumor.bam -o /path/to/output/file.txt
 ```
 对于外显子测序，可使用以下不太严格的参数设置：
 ```R
 -mrq 20.0
 -mlq 25
 -mlc 20
 -mrr 1
 ```
 结果显示3个样本都是MSS 
 ** VisualMSI **  
 通过模仿PCR来检测MSI。VisualMSI从参考基因组中提取PCR adapter，并将其map到测序read中。加入adapter可以成功的map到raed上，将计算这些adapter的插入长度。VisualMSI统计计算这些长度的分布情况。
 ```
 visualmsi -i tumor.sorted.bam -n normal.sorted.bam -r hg19.fasta -t targets/msi.tsv
 ```
 关于target file
 VisualMSI官网提供的参考文件为：
```
#CHROM  POSITION  NAME
chr4  55598216  BAT25
chr2  47641568  BAT26
chr14 23652365  NR-21
chr11 102193518 NR-27
chr2  95849372  NR-24
```
但是他是根据hg19参考基因组来做的，我们需要将它转换为hg38参考基因组的,可以利用UCSC的liftover工具
- 下载注释文件 https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/liftOver/	中的hg19ToHg38.over.chain.gz文件
- 输入文件需要是标准的bed格式

```
必须的字段：
- 染色体
- 起始位置（第一个氨基酸为0）
- 终止位置
chr1  213941196  213942363
chr1  213942363  213943530
chr1  213943530  213944697
chr2  158364697  158365864
chr2  158365864  158367031
可选字段：name、score、strand....
```
- 进行转换
```
../softwave/liftOver CRC_cancer-Specific_marker.bed hg19ToHg38.over.chain.gz hg38_CRC_cancer-Specific_marker.tsv unmap
```
不知道为什么大多数的position都显示能通过过滤标准的reads特别的少，不能通过质量控制，这个软件不适合？还是说参数设置不对？  
** mSINGS **
1) Catalog microsatellites present in your host genome.  There are a number of algorithms available to to this, but one we have found particularly easy to use is MISA (http://pgrc.ipk-gatersleben.de/misa/). 
2) Limit the list of microsatellites to those which are present in your capture design.  This can be done by converting the location of identified microsatellites from step #1 to BED file format, converting the coordinates of your capture design to BED file format, and using the BEDTOOLS intersect function to find where these two files overlap (http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html).  
3) Generate support files and a baseline file from no fewer than 10-20 specimens, as detailed in the readme.
4) Examine the baseline file.  Identify any loci which have an average peak count of 1.0 / standard deviation of 0.0, as these loci may cause artifacts.  These should be removed from the baseline file and the original locus BED file.  The intervals file will need to be re-generated with the edited input BED file.
5) Perform mSINGS analysis of the individual files that you used to generate your baseline with the edited files.  Any remaining loci which are called as unstable in a measurable fraction (~10% or more) of known negative specimens should also be dropped due to their potential to cause artifacts.  Remove these loci and edit files as in step #4
6) Optional - if you have access to a number of known MSI positive specimens (10 or more), it can be helpful to analyze them at this point.  Discriminatory power of mSINGS can be improved if uninformative loci (those which are not unstable in 1 or more known MSI specimens) are removed from the panel.  You may remove uninformative loci as in step #4.  If these specimens are not available, be aware that an empirically-determined cutoff for discriminating MSS form MSI-H specimens may need to be established for your assay.  

mSINGs的使用流程：
```R
##由于mSINGs是基于Python 2.7写的，所以首先要将shell切换到Python 2.7环境
source activate py2.7
git clone https://bitbucket.org/uwlabmed/msings.git
cd msings
bash ./dev/bootstrap.sh
##切换到msings-env环境
source msings-env/bin/activate
##利用MANTIS软件中的获得的filter_bed（capture.bed的交集）文件进行下一步分析
./scripts/create_intervals.sh ../../reference/mSINGS/MANTIS_filtered.bed
##创建msi_bed文件
./scripts/create_baseline.sh ../../reference/mSINGS/normal_bamlist ../../reference/mSINGS/MANTIS_filtered.msi_intervals ../../reference/mSINGS/MANTIS_filtered.bed /data/reference/GRCh38/GRCh38.primary_assembly.genome.fa
```
上面的到`MSI_BASELINE.txt`文件，筛选到其中`Standard_Deviation=0`的Position，这些Position被认为是没有意义的误差，然后与之前获得的`bed`文件取交集，获得有意义的bed文件  
对normal样本进行MSI鉴定，对bed文件进行进一步的验证筛选
```
./scripts/run_msings.sh
          ../../reference/mSINGS/cancer_bamlist ../../reference/mSINGS/MANTIS_filtered.filterSD.msi_intervals
          ../../reference/mSINGS/MANTIS_filtered.filterSD.bed 
          /data/reference/GRCh38/GRCh38.primary_assembly.genome.fa 
          ../../reference/mSINGS/MSI_BASELINE_filterSD.txt
```
根据获的`Combined_MSI.txt`结果，再从`MSI_BASELIEN_filterSD.txt`文件中删除掉在instable_num/total_num > 10% 的position，同样的`bed`文件和`interval`也做相同的处理,得到最终的`MSI_BASELINE_filterSD_checknorm.txt`  
最后进行cancer样本的MSI鉴定
```
./scripts/run_msings.sh
          ../../reference/mSINGS/cancer_bamlist
          ../../reference/mSINGS/MANTIS_filtered.filterSD.checknorm.msi_intervals
          ../../reference/mSINGS/MANTIS_filtered.filterSD.checknorm.bed
          /data/reference/GRCh38/GRCh38.primary_assembly.genome.fa 
          ../../reference/mSINGS/MSI_BASELINE_filterSD_checknorm.txt
```  
结果显示3个样本都是MSS的