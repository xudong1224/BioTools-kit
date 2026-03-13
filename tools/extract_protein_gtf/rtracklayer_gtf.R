#genecode下载gtf文件https://www.gencodegenes.org/human/
##安装包BiocManager::install("rtracklayer")
library(rtracklayer)
library(dplyr)

setwd("C:/Users/DELL/Desktop/")

gtf <- rtracklayer::import('Rattus_norvegicus.Rnor_6.0.104.chr.gtf')

col_target <- grep("gene.*type", names(test), value = TRUE) #字符呈现

protein_coding.gtf <- gtf |> 
  as.data.frame() |> 
  filter(.data[[col_target]] == "protein_coding")

protein_coding.gtf <- makeGRangesFromDataFrame(
  protein_coding.gtf,
  keep.extra.columns = TRUE,
  ignore.strand = FALSE
)

export(protein_coding.gtf, con = "Rattus_norvegicus_protein_coding.gtf", format = "GTF")

#获取编码蛋白

#先提取蛋白编码基因，再去重复

protein_coding.gtf <- gtf |> 
  filter(type == "gene",gene_type == "protein_coding") |> 
  distinct(gene_name,.keep_all = T) |> 
  dplyr::select(c("gene_id","gene_name","gene_type"))
  
#同理可以获取lnRNA、miRNA

lncRNA.gtf <- gtf |> 
  filter(type == "gene",gene_type == "lncRNA") |>
  distinct(gene_name,.keep_all = T) |> 
  dplyr::select(c("gene_id","gene_name","gene_type"))

miRNA.gtf <- gtf |> 
  filter(type == "gene",gene_type == "miRNA") |>
  distinct(gene_name,.keep_all = T) |> 
  dplyr::select(c("gene_id","gene_name","gene_type"))

#获取基因symtol和对应的基因类型
# 同一个 gene_name 对应多个 gene_type 时，不能只对 gene_name 去重，否则会丢失部分类型信息

gene_symtol.gtf <- gtf |> 
  dplyr::select(c("gene_name", "gene_type")) |> 
  arrange(gene_name,gene_type) |> 
  distinct(gene_name,gene_type, .keep_all = T)
