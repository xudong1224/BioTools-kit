#!/project/Research_Service/software/miniconda3/envs/R-v4.3.2/bin/Rscript

library(janitor)
library(dplyr)
library(data.table)

metadata <- fread("meta.xls", sep = "\t", header = T)
print(metadata) |> head()
metadata_clean <- janitor::clean_names(metadata)

generate_merge_cmd <- function(metadata, path,orgin, result) {
  # 1. 内部不需要给 metadata 加 {{ }}，因为它是一个整体的数据框对象
  metadata |>
    # 2. 在 mutate 和 group_by 中，凡是涉及列名的都要用 {{ }}
    mutate(
      r1 = paste0({{ path }}, "/" ,{{ orgin }}, "_R1.fq.gz"),
      r2 = paste0({{ path }}, "/" ,{{ orgin }}, "_R2.fq.gz")
    ) |>
    group_by({{ result }}) |>
    summarise(
      r1_files = paste(r1, collapse = " "), # 这里的 r1 是刚生成的中间列，直接写
      r2_files = paste(r2, collapse = " "),
      cmd_r1 = paste0("cat ", r1_files, " > ", unique({{ result }}), "_R1.fq.gz"),
      cmd_r2 = paste0("cat ", r2_files, " > ", unique({{ result }}), "_R2.fq.gz")
    ) |>
    # 3. pull 只能提取一列！如果要提取两列命令，建议用 c() 结合 pull
    # 或者直接把两列命令合并
    mutate(all_cmds = paste(cmd_r1, cmd_r2, sep = "\n")) |>
    pull(all_cmds)
}

# 调用方式（不需要引号）：
sh_commands <- generate_merge_cmd(metadata_clean, shu_ju_lu_jing ,xia_ji_ming_cheng, jiao_fu_ming_cheng)

writeLines(sh_commands, "merge_samples.sh")

