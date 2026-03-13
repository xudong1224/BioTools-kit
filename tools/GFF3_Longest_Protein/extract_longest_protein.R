#!/project/Research_Service/software/miniconda3/envs/R-v4.3.2/bin/Rscript

pkgs <- c("argparse", "tidyverse", "Biostrings", "readr")

invisible(lapply(pkgs, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE))
}))

# --- 1. 参数解析函数 ---
parse_args <- function() {
  parser <- ArgumentParser(description = "从 GFF3 和 FASTA 中提取每个基因最长的蛋白质序列 (以 GeneID 命名)")
  
  # 核心输入输出参数
  parser$add_argument("-g", "--gff", required = TRUE, 
                      help = "输入的 GFF3 注释文件")
  parser$add_argument("-f", "--fasta", required = TRUE, 
                      help = "输入的基因组 FASTA 文件")
  parser$add_argument("-o", "--output", default = "protein_longest_geneID.fa", 
                      help = "输出的 FASTA 文件名 [默认: %(default)s]")
  
  # 临时文件与性能参数
  parser$add_argument("-t", "--temp", default = "all_prots_temp.fa", 
                      help = "gffread 生成的临时蛋白文件 [默认: %(default)s]")
  
  parser$parse_args()
}

# --- 2. 核心逻辑函数 ---
main <- function() {
  # 解析参数
  args <- parse_args()
  
  # A. 提取 ID 映射关系
  message(">> [1/4] 正在解析 GFF3 映射关系...")
  transcript_to_gene <- read_tsv(args$gff, comment = "#", col_names = FALSE, show_col_types = FALSE) |>
    filter(X3 == "mRNA") |> 
    mutate(
      transcript_id = str_extract(X9, "(?<=ID=)[^;]+"),
      gene_id = str_extract(X9, "(?<=Parent=)[^;]+")
    ) |> 
    select(transcript_id, gene_id)
  
  # B. 调用 gffread 提取所有序列
  message(">> [2/4] 正在运行 gffread 提取原始序列...")
  system(paste("gffread", args$gff, "-g", args$fasta, "-y", args$temp))
  
  # C. 筛选最长序列 (dplyr 链式操作)
  message(">> [3/4] 正在筛选每个基因的最长异构体...")
  sequences <- readAAStringSet(args$temp)
  
  final_df <- data.frame(
    transcript_id = names(sequences) |> word(1),
    width = width(sequences),
    seq = as.character(sequences)
  ) |> 
    inner_join(transcript_to_gene, by = "transcript_id") |> 
    group_by(gene_id) |> 
    slice_max(order_by = width, n = 1, with_ties = FALSE) |> 
    ungroup()
  
  # D. 格式转换与写出
  message(">> [4/4] 正在写出最终文件: ", args$output)
  final_sequences_longest <- AAStringSet(final_df$seq) |> 
    setNames(final_df$gene_id)
  
  writeXStringSet(final_sequences_longest, args$output)
  
  # 清理
  if (file.exists(args$temp)) file.remove(args$temp)
  message(">> 处理完成！最终提取基因数: ", length(final_sequences_longest))
  message("是否存在重复 ID：", any(duplicated(names(final_sequences_longest))))
}

# --- 3. 执行脚本 ---
main()