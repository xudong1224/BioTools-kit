#genecode下载gtf文件https://www.gencodegenes.org/human/
##安装包BiocManager::install("rtracklayer")
library(rtracklayer)
library(dplyr)

setwd("C:/Users/DELL/Desktop/gtf")

gtf <- rtracklayer::import('gencode.v46.annotation.gtf')
gtf<-as.data.frame(gtf)

#获取编码蛋白
protein_coding.gtf <- gtf %>%
  dplyr::select(c("gene_id", "gene_name", "gene_type")) %>%
  arrange(gene_id) %>%   
  distinct(gene_name, .keep_all = T) %>%
  filter(gene_type == "protein_coding")

#同理可以获取lnRNA、miRNA
lncRNA.gtf<-gtf %>% 
  dplyr::select(c("gene_id","gene_name","gene_type")) %>%   
  arrange(gene_id) %>%   
  distinct(gene_name,.keep_all = T) %>%   
  filter(gene_type == "lncRNA")

miRNA.gtf<-gtf %>%   
  dplyr::select(c("gene_id","gene_name","gene_type")) %>%   
  arrange(gene_id) %>%   
  distinct(gene_name,.keep_all = T) %>%
  filter(gene_type == "miRNA")

#获取基因symtol和对应的基因类型
gene_symtol.gtf <- gtf %>%
  dplyr::select(c("gene_name", "gene_type")) %>%
  arrange(gene_name) %>%   
  distinct(gene_name, .keep_all = T)




#对基因表达谱进行注释
library(readr)

gene_matrix_fpkm <- read_delim("gene_count_matrix.fpkm.xls",delim = "\t",
                                     escape_double = FALSE,trim_ws = TRUE)

colnames(gene_matrix_fpkm)[1] <- "gene_name"

fpkm <- gene_matrix_fpkm |> 
  left_join(gene_symtol.gtf,by = "gene_name") |> 
  dplyr::select(gene_name,gene_type,everything())

write.table(fpkm,"gene_matrix.fpkm.xls",row.names = F,sep = "\t")


