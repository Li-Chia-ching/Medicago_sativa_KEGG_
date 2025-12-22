# ==============================================================================
# GO/KEGG富集分析V3.1 - 完整R脚本（修正字符串连接错误版）
# 作者：基于专家评审意见改进
# 日期：2024
# 功能：基于本地注释数据的GO/KEGG富集分析
# ==============================================================================

# ==============================================================================
# 0. 安装和加载包
# ==============================================================================

cat("正在检查并安装必要的R包...\n")

# 定义所需包列表
required_packages <- c("clusterProfiler", "ggplot2", "stringr", "dplyr", "tidyr", 
                       "GO.db", "ggrepel", "scales", "RColorBrewer", "forcats",
                       "ggthemes", "grid", "gridExtra", "readxl")

# 安装CRAN包
install_cran_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(paste("安装包:", pkg, "\n"))
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

# 安装Bioconductor包
install_bioc_packages <- function(packages) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(paste("安装Bioconductor包:", pkg, "\n"))
      BiocManager::install(pkg, update = FALSE)
      library(pkg, character.only = TRUE)
    }
  }
}

# 安装CRAN包
cran_packages <- c("ggplot2", "stringr", "dplyr", "tidyr", "ggrepel", "scales", 
                   "RColorBrewer", "forcats", "ggthemes", "grid", "gridExtra", "readxl")
install_cran_packages(cran_packages)

# 安装Bioconductor包
bioc_packages <- c("clusterProfiler", "GO.db")
install_bioc_packages(bioc_packages)

cat("所有必要包已安装和加载完成！\n")

# ==============================================================================
# 1. 创建输出目录
# ==============================================================================

# 创建输出目录
output_dir <- paste0("GO_KEGG_Analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(output_dir, showWarnings = FALSE)
cat(paste("输出目录:", output_dir, "\n"))

# 设置工作目录到输出目录
original_wd <- getwd()
setwd(output_dir)

# ==============================================================================
# 2. 数据加载和预处理
# ==============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("数据加载和预处理\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 检查数据是否存在
if (!exists("Summary_Modules")) {
  cat("未在环境中找到'Summary_Modules'数据框\n")
  cat("尝试从文件加载数据...\n")
  
  # 尝试常见的数据文件
  possible_files <- c("Summary_Modules.xlsx", "Summary_Modules.csv", 
                      "Summary_Modules.txt", "data.xlsx", "data.csv")
  
  file_found <- FALSE
  for (file in possible_files) {
    if (file.exists(file.path(original_wd, file))) {
      cat(paste("找到文件:", file, "\n"))
      
      if (grepl("\\.xlsx$", file)) {
        df <- readxl::read_excel(file.path(original_wd, file))
      } else {
        df <- read.csv(file.path(original_wd, file), stringsAsFactors = FALSE)
      }
      
      file_found <- TRUE
      break
    }
  }
  
  if (!file_found) {
    stop("未找到数据文件。请确保数据已加载或文件存在于工作目录中。")
  }
} else {
  df <- Summary_Modules
  cat("使用环境中的'Summary_Modules'数据框\n")
}

# 检查数据基本结构
cat(paste("数据维度:", nrow(df), "行 ×", ncol(df), "列\n"))
cat(paste("列名:", paste(colnames(df), collapse = ", "), "\n"))

# 检查必需列
required_cols <- c("Gene_ID", "Module", "GO", "KEGG")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  # 尝试查找可能的替代列名
  cat("正在查找可能的替代列名...\n")
  
  colname_mapping <- list(
    "Gene_ID" = c("GeneID", "gene_id", "gene", "ID", "Name", "Gene"),
    "Module" = c("module", "Cluster", "cluster", "Group", "group"),
    "GO" = c("Go", "go", "GO_terms", "go_terms"),
    "KEGG" = c("Kegg", "kegg", "Pathway", "pathway", "KEGG_pathway")
  )
  
  for (req_col in missing_cols) {
    possible_names <- colname_mapping[[req_col]]
    for (possible_name in possible_names) {
      if (possible_name %in% colnames(df)) {
        df[[req_col]] <- df[[possible_name]]
        cat(paste("  使用", possible_name, "作为", req_col, "\n"))
        missing_cols <- setdiff(missing_cols, req_col)
        break
      }
    }
  }
}

if (length(missing_cols) > 0) {
  stop(paste("缺少必需列:", paste(missing_cols, collapse = ", ")))
}

# 基本数据清理
df$Gene_ID <- as.character(df$Gene_ID)
df$Module <- as.character(df$Module)

# 清理基因ID（移除版本号）
df$Gene_ID_Clean <- stringr::str_replace(df$Gene_ID, "\\..*", "")

# ==============================================================================
# 3. 特殊模块处理（根据专家建议）
# ==============================================================================

cat("\n特殊模块处理\n")

# 定义需要特殊处理的模块（根据专家建议）
special_modules <- c("grey", "gray", "0", "none", "unassigned", "NA", "")

# 检查是否存在特殊模块
modules <- unique(df$Module)
modules <- modules[!is.na(modules) & modules != ""]
found_special <- intersect(tolower(special_modules), tolower(modules))

if (length(found_special) > 0) {
  cat(paste("发现特殊模块:", paste(found_special, collapse = ", "), "\n"))
  cat("这些模块通常包含未聚类或表达量低的基因\n")
  
  # 显示特殊模块的基因数
  special_module_stats <- df %>%
    dplyr::filter(tolower(Module) %in% found_special) %>%
    dplyr::group_by(Module) %>%
    dplyr::summarise(Gene_Count = n(), .groups = 'drop')
  
  cat("特殊模块基因数统计:\n")
  print(special_module_stats)
  
  # 询问用户是否要排除这些模块
  cat("\n是否要排除这些特殊模块进行分析？(y/n): ")
  user_input <- readline(prompt = "")
  
  if (tolower(user_input) %in% c("y", "yes", "是")) {
    original_count <- nrow(df)
    df <- df[!tolower(df$Module) %in% found_special, ]
    removed_count <- original_count - nrow(df)
    cat(paste("已排除", removed_count, "个基因\n"))
    
    # 更新模块列表
    modules <- unique(df$Module)
    modules <- modules[!is.na(modules) & modules != ""]
  } else {
    cat("保留所有模块进行分析\n")
  }
}

# 移除没有模块的基因
df <- df[!is.na(df$Module) & df$Module != "" & df$Module != "-", ]

# 更新模块列表
modules <- unique(df$Module)
modules <- modules[!is.na(modules) & modules != ""]

cat(paste("最终模块数量:", length(modules), "\n"))
cat(paste("模块列表:", paste(modules, collapse = ", "), "\n"))

# 统计模块大小
module_stats <- df %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    Gene_Count = n(),
    .groups = 'drop'
  ) %>%
  dplyr::arrange(desc(Gene_Count))

cat("\n模块基因数统计:\n")
print(module_stats)

# 保存模块统计
write.csv(module_stats, "01_Module_Gene_Statistics.csv", row.names = FALSE)
cat("✓ 保存模块统计到: 01_Module_Gene_Statistics.csv\n")

# ==============================================================================
# 4. 数据质量诊断
# ==============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("数据质量诊断\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 检查GO列的内容和质量
go_values <- df$GO[!is.na(df$GO) & df$GO != "" & df$GO != "-"]
go_nonempty_count <- length(go_values)
go_coverage <- round(go_nonempty_count / nrow(df) * 100, 1)

cat(paste("GO列非空值数量:", go_nonempty_count, "/", nrow(df), 
          paste0("(", go_coverage, "%)"), "\n"))

if (go_nonempty_count > 0) {
  # 显示前2个非空GO值作为示例
  cat("GO列示例值:\n")
  for (i in 1:min(2, length(go_values))) {
    sample_text <- as.character(go_values[i])
    # 截断过长的文本
    if (nchar(sample_text) > 80) {
      sample_text <- paste0(substr(sample_text, 1, 77), "...")
    }
    cat(paste("  ", i, ":", sample_text, "\n"))
  }
  
  # 分析GO值结构
  sample_go <- as.character(go_values[1])
  cat(paste("GO值结构分析:\n"))
  cat(paste("  长度:", nchar(sample_go), "\n"))
  cat(paste("  是否包含'GO:'前缀:", grepl("GO:", sample_go), "\n"))
  cat(paste("  是否包含分号:", grepl(";", sample_go), "\n"))
  
  # 检查常见的分隔符
  separators <- c(";", ",", "\\|")
  for (sep in separators) {
    if (grepl(sep, sample_go)) {
      cat(paste("  使用分隔符:", sep, "\n"))
      break
    }
  }
}

# 检查KEGG列的内容和质量
kegg_values <- df$KEGG[!is.na(df$KEGG) & df$KEGG != "" & df$KEGG != "-"]
kegg_nonempty_count <- length(kegg_values)
kegg_coverage <- round(kegg_nonempty_count / nrow(df) * 100, 1)

cat(paste("\nKEGG列非空值数量:", kegg_nonempty_count, "/", nrow(df), 
          paste0("(", kegg_coverage, "%)"), "\n"))

if (kegg_nonempty_count > 0) {
  # 显示前2个非空KEGG值作为示例
  cat("KEGG列示例值:\n")
  for (i in 1:min(2, length(kegg_values))) {
    sample_text <- as.character(kegg_values[i])
    # 截断过长的文本
    if (nchar(sample_text) > 80) {
      sample_text <- paste0(substr(sample_text, 1, 77), "...")
    }
    cat(paste("  ", i, ":", sample_text, "\n"))
  }
  
  # 分析KEGG值结构
  sample_kegg <- as.character(kegg_values[1])
  cat(paste("KEGG值结构分析:\n"))
  cat(paste("  长度:", nchar(sample_kegg), "\n"))
  
  # 检查常见的KEGG标识符格式
  kegg_patterns <- c("mtr\\d{5}", "map\\d{5}", "ko\\d{5}", "K\\d{5}")
  for (pattern in kegg_patterns) {
    if (grepl(pattern, sample_kegg, ignore.case = TRUE)) {
      cat(paste("  包含KEGG模式:", pattern, "\n"))
    }
  }
}

# 数据质量总结
cat("\n数据质量总结:\n")
cat(paste("  总基因数:", nrow(df), "\n"))
cat(paste("  GO注释覆盖率:", go_coverage, "%\n"))
cat(paste("  KEGG注释覆盖率:", kegg_coverage, "%\n"))

# 保存数据质量报告
quality_report <- data.frame(
  指标 = c("总基因数", "GO注释覆盖率", "KEGG注释覆盖率"),
  值 = c(nrow(df), paste0(go_coverage, "%"), paste0(kegg_coverage, "%"))
)
write.csv(quality_report, "02_Data_Quality_Report.csv", row.names = FALSE)

# ==============================================================================
# 5. 从原始数据提取GO注释（改进版）
# ==============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("GO注释提取\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

extract_go_from_data <- function(df) {
  go_list <- list()
  genes_with_go <- 0
  
  cat("正在提取GO注释...\n")
  
  for (i in 1:nrow(df)) {
    gene_id <- df$Gene_ID_Clean[i]
    go_text <- as.character(df$GO[i])
    
    # 跳过空值
    if (is.na(go_text) || go_text == "" || go_text == "-") {
      next
    }
    
    # 尝试多种分隔符（根据专家建议增加灵活性）
    separators <- c(";", ",", "\\|", " ", "，")
    terms <- character()
    
    for (sep in separators) {
      if (grepl(sep, go_text)) {
        terms <- unlist(strsplit(go_text, sep))
        break
      }
    }
    
    # 如果没有找到分隔符，将整个文本作为一个项
    if (length(terms) == 0) {
      terms <- go_text
    }
    
    terms <- trimws(terms)
    terms <- terms[terms != ""]
    
    if (length(terms) > 0) {
      # 提取有效的GO ID（格式：GO:XXXXXXX，支持7位数字）
      go_ids <- character()
      for (term in terms) {
        # 查找GO:后跟7位数字的模式（根据专家建议）
        # 修改正则表达式改为匹配“GO:”后跟任意长度数字的更通用模式
        matches <- regmatches(term, regexpr("GO:\\d+", term, perl = TRUE))
        if (length(matches) > 0) {
          go_ids <- c(go_ids, matches)
        }
      }
      
      if (length(go_ids) > 0) {
        go_list[[gene_id]] <- unique(go_ids)
        genes_with_go <- genes_with_go + 1
      }
    }
    
    # 进度显示
    if (i %% 1000 == 0) {
      cat(paste("  已处理", i, "/", nrow(df), "个基因\n"))
    }
  }
  
  cat(paste("  提取到", genes_with_go, "个基因的GO注释\n"))
  
  if (length(go_list) == 0) {
    cat("  警告: 未提取到任何GO注释\n")
    return(data.frame(Gene = character(), GO_ID = character(), GO_Term = character(), Ontology = character()))
  }
  
  # 转换为数据框
  go_df <- data.frame(
    Gene = rep(names(go_list), sapply(go_list, length)),
    GO_ID = unlist(go_list),
    stringsAsFactors = FALSE
  )
  
  cat(paste("  总GO术语数:", nrow(go_df), "\n"))
  
  # 使用GO.db获取GO术语详情（根据专家建议增加错误处理）
  cat("  正在获取GO术语详细信息...\n")
  
  tryCatch({
    unique_go_ids <- unique(go_df$GO_ID)
    cat(paste("  唯一GO ID数量:", length(unique_go_ids), "\n"))
    
    # 检查GO.db包是否可用（根据专家建议）
    if (!requireNamespace("GO.db", quietly = TRUE)) {
      cat("  警告: GO.db包未安装，将无法获取GO术语的完整描述和本体分类\n")
      cat("  建议安装: BiocManager::install('GO.db')\n")
      go_df$GO_Term <- go_df$GO_ID
      go_df$Ontology <- "UNKNOWN"
      return(go_df)
    }
    
    # 获取GO术语和本体信息
    go_info <- AnnotationDbi::select(GO.db, 
                                     keys = unique_go_ids,
                                     columns = c("TERM", "ONTOLOGY"),
                                     keytype = "GOID",
                                     multiVals = "first")
    
    if (nrow(go_info) > 0) {
      # 检查是否有未匹配的GO ID
      unmatched_ids <- setdiff(unique_go_ids, go_info$GOID)
      if (length(unmatched_ids) > 0) {
        cat(paste("  注意:", length(unmatched_ids), "个GO ID在GO.db中未找到\n"))
        cat("  这可能是因为GO ID格式不正确或版本过旧\n")
      }
      
      # 合并信息
      go_df <- merge(go_df, go_info, by.x = "GO_ID", by.y = "GOID", all.x = TRUE)
      colnames(go_df) <- c("GO_ID", "Gene", "GO_Term", "Ontology")
      go_df <- go_df[, c("Gene", "GO_ID", "GO_Term", "Ontology")]
      
      # 填充未匹配的条目
      go_df$GO_Term[is.na(go_df$GO_Term)] <- go_df$GO_ID[is.na(go_df$GO_Term)]
      go_df$Ontology[is.na(go_df$Ontology)] <- "UNKNOWN"
      
      cat("  GO术语信息获取成功\n")
      
    } else {
      cat("  警告: GO.db查询返回空结果，请检查GO ID格式\n")
      go_df$GO_Term <- go_df$GO_ID
      go_df$Ontology <- "UNKNOWN"
    }
    
  }, error = function(e) {
    cat(paste("  错误: 无法从GO.db获取详细信息 -", e$message, "\n"))
    cat("  将使用GO ID作为术语名称\n")
    go_df$GO_Term <- go_df$GO_ID
    go_df$Ontology <- "UNKNOWN"
  })
  
  return(go_df)
}

go_annotation <- extract_go_from_data(df)

# 保存GO注释
if (nrow(go_annotation) > 0) {
  write.csv(go_annotation, "03_GO_Annotations_Extracted.csv", row.names = FALSE)
  cat(paste("✓ 保存GO注释到: 03_GO_Annotations_Extracted.csv\n"))
  
  # 统计GO本体分布
  if ("Ontology" %in% colnames(go_annotation)) {
    ontology_stats <- table(go_annotation$Ontology)
    cat("  GO本体分布:\n")
    print(ontology_stats)
    
    write.csv(as.data.frame(ontology_stats), "03_GO_Ontology_Distribution.csv", row.names = TRUE)
  }
} else {
  cat("  警告: 未提取到GO注释\n")
}

# ==============================================================================
# 6. 从原始数据提取KEGG注释（改进版，根据专家建议）
# ==============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("KEGG注释提取\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

extract_kegg_from_data <- function(df) {
  kegg_list <- list()
  genes_with_kegg <- 0
  
  cat("正在提取KEGG注释...\n")
  
  for (i in 1:nrow(df)) {
    gene_id <- df$Gene_ID_Clean[i]
    kegg_text <- as.character(df$KEGG[i])
    
    # 跳过空值
    if (is.na(kegg_text) || kegg_text == "" || kegg_text == "-") {
      next
    }
    
    # 尝试多种分隔符
    separators <- c(";", ",", "\\|", " ", "，")
    pathways <- character()
    
    for (sep in separators) {
      if (grepl(sep, kegg_text)) {
        pathways <- unlist(strsplit(kegg_text, sep))
        break
      }
    }
    
    # 如果没有找到分隔符，将整个文本作为一个项
    if (length(pathways) == 0) {
      pathways <- kegg_text
    }
    
    pathways <- trimws(pathways)
    pathways <- pathways[pathways != ""]
    
    if (length(pathways) > 0) {
      # 提取KEGG通路ID（根据专家建议支持多种格式）
      pathway_ids <- character()
      for (pathway in pathways) {
        # 支持多种KEGG格式（根据专家建议改进正则表达式）
        # KEGG ID提取规则改为更灵活的模式，并确保正确捕获
        pattern <- "\\b(([a-zA-Z]{2,6})0\\d{4,5})\\b"
        matches <- regmatches(pathway, regexpr(pattern, pathway, ignore.case = TRUE, perl = TRUE))
        
        if (length(matches) > 0) {
          # 清理匹配结果：移除PATH:前缀，统一为小写
          clean_match <- tolower(gsub("PATH:", "", matches[1]))
          pathway_ids <- c(pathway_ids, clean_match)
        } else {
          # 如果没有匹配标准格式，但包含字母和数字，可能是一个KEGG通路名称
          if (grepl("[a-zA-Z].*\\d|\\d.*[a-zA-Z]", pathway)) {
            pathway_ids <- c(pathway_ids, pathway)
          }
        }
      }
      
      if (length(pathway_ids) > 0) {
        kegg_list[[gene_id]] <- unique(pathway_ids)
        genes_with_kegg <- genes_with_kegg + 1
      }
    }
    
    # 进度显示
    if (i %% 1000 == 0) {
      cat(paste("  已处理", i, "/", nrow(df), "个基因\n"))
    }
  }
  
  cat(paste("  提取到", genes_with_kegg, "个基因的KEGG注释\n"))
  
  if (length(kegg_list) == 0) {
    cat("  警告: 未提取到任何KEGG注释\n")
    return(data.frame(Gene = character(), Pathway_ID = character()))
  }
  
  # 转换为数据框
  kegg_df <- data.frame(
    Gene = rep(names(kegg_list), sapply(kegg_list, length)),
    Pathway_ID = unlist(kegg_list),
    stringsAsFactors = FALSE
  )
  
  cat(paste("  总KEGG通路数:", nrow(kegg_df), "\n"))
  
  # 注意：这需要网络连接，如果希望纯本地，可告知用户手动准备映射表
  cat("  正在获取KEGG通路名称...\n")
  tryCatch({
    unique_pathways <- unique(kegg_df$Pathway_ID)
    # 使用clusterProfiler的辅助函数获取名称
    pathway_names <- clusterProfiler::pathway2name(unique_pathways, organism = "mtr")
    if (!is.null(pathway_names) && nrow(pathway_names) > 0) {
      kegg_df <- merge(kegg_df, pathway_names, by.x = "Pathway_ID", by.y = "pathway", all.x = TRUE)
      # 重命名列
      colnames(kegg_df)[colnames(kegg_df) == "name"] <- "Pathway_Name"
    }
  }, error = function(e) {
    cat(paste("  警告: 无法获取KEGG通路名称 -", e$message, "\n"))
  })
  
  return(kegg_df)
}

kegg_annotation <- extract_kegg_from_data(df)

# 保存KEGG注释
if (nrow(kegg_annotation) > 0) {
  write.csv(kegg_annotation, "04_KEGG_Annotations_Extracted.csv", row.names = FALSE)
  cat(paste("✓ 保存KEGG注释到: 04_KEGG_Annotations_Extracted.csv\n"))
} else {
  cat("  警告: 未提取到KEGG注释\n")
}

# ==============================================================================
# 7. 富集分析函数（改进版）
# ==============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("富集分析准备\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 检查是否有足够的注释数据进行分析
if (nrow(go_annotation) == 0 && nrow(kegg_annotation) == 0) {
  cat("错误: 没有提取到任何注释数据，无法进行富集分析\n")
  cat("请检查您的数据中GO和KEGG列是否有有效的注释信息\n")
} else {
  # 准备富集分析函数（增加错误处理和调试信息）
  perform_enrichment_analysis <- function(gene_list, annotation_df, 
                                          term_col, term_name_col = NULL,
                                          p_cutoff = 0.05, 
                                          min_size = 5,
                                          max_size = 500,
                                          description = "富集分析") {
    
    cat(paste("    进行", description, "...\n"))
    
    if (length(gene_list) < min_size) {
      cat(paste("      基因数太少 (", length(gene_list), " < ", min_size, ")，跳过\n", sep = ""))
      return(NULL)
    }
    
    # 检查注释数据
    if (nrow(annotation_df) == 0) {
      cat("      注释数据为空\n")
      return(NULL)
    }
    
    # 准备TERM2GENE映射
    if (is.null(term_name_col)) {
      term2gene <- data.frame(
        term = annotation_df[[term_col]],
        gene = annotation_df$Gene,
        stringsAsFactors = FALSE
      )
    } else {
      term2gene <- data.frame(
        term = annotation_df[[term_col]],
        gene = annotation_df$Gene,
        stringsAsFactors = FALSE
      )
    }
    
    # 去除重复和空值
    term2gene <- term2gene[!is.na(term2gene$term) & !is.na(term2gene$gene), ]
    term2gene <- term2gene[!duplicated(term2gene), ]
    
    if (nrow(term2gene) == 0) {
      cat("      无有效的基因-术语映射\n")
      return(NULL)
    }
    
    # 准备TERM2NAME映射（如果有）
    if (!is.null(term_name_col) && term_name_col %in% colnames(annotation_df)) {
      term2name <- data.frame(
        term = annotation_df[[term_col]],
        name = annotation_df[[term_name_col]],
        stringsAsFactors = FALSE
      )
      term2name <- term2name[!duplicated(term2name), ]
      term2name <- term2name[!is.na(term2name$term) & !is.na(term2name$name), ]
    } else {
      term2name <- data.frame(
        term = unique(term2gene$term),
        name = unique(term2gene$term),
        stringsAsFactors = FALSE
      )
    }
    
    # 检查基因列表与注释数据的重叠
    overlapping_genes <- intersect(gene_list, unique(term2gene$gene))
    if (length(overlapping_genes) == 0) {
      cat(paste("      基因列表与注释数据无重叠\n"))
      return(NULL)
    }
    
    cat(paste("      基因列表与注释数据重叠数:", length(overlapping_genes), "/", length(gene_list), "\n"))
    
    # 进行富集分析
    tryCatch({
      enrich_result <- clusterProfiler::enricher(
        gene = gene_list,
        pvalueCutoff = p_cutoff,
        pAdjustMethod = "BH",
        minGSSize = min_size,
        maxGSSize = max_size,
        TERM2GENE = term2gene,
        TERM2NAME = term2name
      )
      
      if (is.null(enrich_result) || nrow(enrich_result) == 0) {
        cat("      无显著富集结果\n")
        return(NULL)
      }
      
      cat(paste("      找到", nrow(enrich_result), "个显著富集项\n"))
      return(enrich_result)
      
    }, error = function(e) {
      cat(paste("      富集分析错误:", e$message, "\n"))
      return(NULL)
    })
  }
  
  # ==============================================================================
  # 8. 进行GO富集分析（按本体）
  # ==============================================================================
  
  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("GO富集分析\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  if (nrow(go_annotation) > 0) {
    
    # 按GO本体分离
    go_bp <- go_annotation[go_annotation$Ontology == "BP", ]
    go_cc <- go_annotation[go_annotation$Ontology == "CC", ]
    go_mf <- go_annotation[go_annotation$Ontology == "MF", ]
    
    all_go_results <- data.frame()
    modules_with_go_results <- c()
    
    for (mod in modules) {
      cat(paste("  模块", mod, ": "))
      module_genes <- df$Gene_ID_Clean[df$Module == mod]
      module_gene_count <- length(module_genes)
      
      if (module_gene_count >= 10) {  # 提高最小基因数要求以获得更可靠的结果
        # 分别分析三个本体
        go_results_for_module <- data.frame()
        
        # BP本体
        if (nrow(go_bp) > 0) {
          bp_result <- perform_enrichment_analysis(
            gene_list = module_genes,
            annotation_df = go_bp,
            term_col = "GO_ID",
            term_name_col = "GO_Term",
            p_cutoff = 0.05,
            description = "BP富集"
          )
          
          if (!is.null(bp_result)) {
            bp_df <- as.data.frame(bp_result)
            bp_df$Module <- mod
            bp_df$Ontology <- "BP"
            bp_df$Module_Gene_Count <- module_gene_count
            go_results_for_module <- rbind(go_results_for_module, bp_df)
          }
        }
        
        # CC本体
        if (nrow(go_cc) > 0) {
          cc_result <- perform_enrichment_analysis(
            gene_list = module_genes,
            annotation_df = go_cc,
            term_col = "GO_ID",
            term_name_col = "GO_Term",
            p_cutoff = 0.05,
            description = "CC富集"
          )
          
          if (!is.null(cc_result)) {
            cc_df <- as.data.frame(cc_result)
            cc_df$Module <- mod
            cc_df$Ontology <- "CC"
            cc_df$Module_Gene_Count <- module_gene_count
            go_results_for_module <- rbind(go_results_for_module, cc_df)
          }
        }
        
        # MF本体
        if (nrow(go_mf) > 0) {
          mf_result <- perform_enrichment_analysis(
            gene_list = module_genes,
            annotation_df = go_mf,
            term_col = "GO_ID",
            term_name_col = "GO_Term",
            p_cutoff = 0.05,
            description = "MF富集"
          )
          
          if (!is.null(mf_result)) {
            mf_df <- as.data.frame(mf_result)
            mf_df$Module <- mod
            mf_df$Ontology <- "MF"
            mf_df$Module_Gene_Count <- module_gene_count
            go_results_for_module <- rbind(go_results_for_module, mf_df)
          }
        }
        
        if (nrow(go_results_for_module) > 0) {
          all_go_results <- rbind(all_go_results, go_results_for_module)
          modules_with_go_results <- c(modules_with_go_results, mod)
          cat(paste("找到", nrow(go_results_for_module), "个显著富集项\n"))
        } else {
          cat("无显著结果\n")
        }
      } else {
        cat(paste("跳过 (仅", module_gene_count, "个基因)\n"))
      }
    }
    
    if (nrow(all_go_results) > 0) {
      write.csv(all_go_results, "05_GO_Enrichment_Results.csv", row.names = FALSE)
      cat(paste("✓ GO富集分析完成，共", nrow(all_go_results), "个显著富集项\n"))
      cat(paste("  有富集结果的模块数:", length(modules_with_go_results), "\n"))
    } else {
      cat("  警告: 无显著GO富集结果\n")
    }
  } else {
    cat("  跳过GO富集分析（无GO注释）\n")
  }
  
  # ==============================================================================
  # 9. 进行KEGG富集分析
  # ==============================================================================
  
  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("KEGG富集分析\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  if (nrow(kegg_annotation) > 0) {
    
    all_kegg_results <- data.frame()
    modules_with_kegg_results <- c()
    
    for (mod in modules) {
      cat(paste("  模块", mod, ": "))
      module_genes <- df$Gene_ID_Clean[df$Module == mod]
      module_gene_count <- length(module_genes)
      
      if (module_gene_count >= 10) {
        kegg_result <- perform_enrichment_analysis(
          gene_list = module_genes,
          annotation_df = kegg_annotation,
          term_col = "Pathway_ID",
          p_cutoff = 0.05,
          description = "KEGG富集"
        )
        
        if (!is.null(kegg_result)) {
          kegg_df <- as.data.frame(kegg_result)
          kegg_df$Module <- mod
          kegg_df$Module_Gene_Count <- module_gene_count
          all_kegg_results <- rbind(all_kegg_results, kegg_df)
          modules_with_kegg_results <- c(modules_with_kegg_results, mod)
          cat(paste("找到", nrow(kegg_df), "个显著通路\n"))
        } else {
          cat("无显著通路\n")
        }
      } else {
        cat(paste("跳过 (仅", module_gene_count, "个基因)\n"))
      }
    }
    
    if (nrow(all_kegg_results) > 0) {
      write.csv(all_kegg_results, "06_KEGG_Enrichment_Results.csv", row.names = FALSE)
      cat(paste("✓ KEGG富集分析完成，共", nrow(all_kegg_results), "个显著富集通路\n"))
      cat(paste("  有富集结果的模块数:", length(modules_with_kegg_results), "\n"))
    } else {
      cat("  警告: 无显著KEGG富集结果\n")
    }
  } else {
    cat("  跳过KEGG富集分析（无KEGG注释）\n")
  }
}

# ==============================================================================
# 10. 高级可视化（根据专家建议改进）
# ==============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("生成高质量可视化图表\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 定义低饱和度配色方案（根据用户要求）
low_saturation_colors <- c("#FDDED7", "#F5BE8F", "#C1E0DB", "#CCD376", 
                           "#A28CC2", "#8498AB", "#5CB0C3")

# 创建颜色梯度函数
create_color_gradient <- function(n = 10) {
  colorRampPalette(c("#FDDED7", "#C1E0DB", "#5CB0C3"))(n)
}

# 文本换行函数（根据专家建议）
add_line_breaks <- function(text, max_length = 50) {
  if (is.na(text) || text == "") return("")
  
  words <- unlist(strsplit(as.character(text), " "))
  lines <- character()
  current_line <- ""
  
  for (word in words) {
    if (nchar(current_line) + nchar(word) + 1 <= max_length) {
      if (current_line == "") {
        current_line <- word
      } else {
        current_line <- paste(current_line, word)
      }
    } else {
      lines <- c(lines, current_line)
      current_line <- word
    }
  }
  
  if (current_line != "") {
    lines <- c(lines, current_line)
  }
  
  paste(lines, collapse = "\n")
}

# GO富集结果可视化
if (exists("all_go_results") && nrow(all_go_results) > 0) {
  cat("  生成GO富集图表...\n")
  
  # 选择每个模块每个本体最多5个最显著的GO项
  top_go <- all_go_results %>%
    dplyr::group_by(Module, Ontology) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 5) %>%
    dplyr::ungroup()
  
  if (nrow(top_go) > 0) {
    # 应用文本换行
    top_go$Wrapped_Description <- sapply(top_go$Description, add_line_breaks)
    
    # 创建图表
    p_go <- ggplot2::ggplot(top_go, ggplot2::aes(x = Module, y = reorder(Wrapped_Description, -p.adjust))) +
      ggplot2::geom_point(ggplot2::aes(size = Count, fill = -log10(p.adjust)), 
                          shape = 21, color = "gray30", stroke = 1.0, alpha = 0.9) +
      ggplot2::scale_fill_gradientn(
        colors = create_color_gradient(10),
        name = expression(-log[10](p.adjust)),
        breaks = scales::pretty_breaks(n = 5)
      ) +
      ggplot2::scale_size_continuous(
        range = c(4, 12), 
        name = "Gene Count",
        breaks = scales::pretty_breaks(n = 5)
      ) +
      ggplot2::facet_grid(Ontology ~ ., scales = "free_y", space = "free") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45, 
          hjust = 1, 
          size = 12,
          face = "bold",
          color = "#333333"
        ),
        axis.text.y = ggplot2::element_text(
          size = 11,
          lineheight = 0.9,
          color = "#444444"
        ),
        axis.title.x = ggplot2::element_text(
          size = 14,
          face = "bold",
          margin = ggplot2::margin(t = 10)
        ),
        axis.title.y = ggplot2::element_text(
          size = 14,
          face = "bold",
          margin = ggplot2::margin(r = 10)
        ),
        strip.text = ggplot2::element_text(
          face = "bold", 
          size = 12,
          color = "white"
        ),
        strip.background = ggplot2::element_rect(
          fill = "#666666",
          color = "#666666"
        ),
        panel.grid.major = ggplot2::element_line(
          color = "#E8E8E8",
          linewidth = 0.5
        ),
        panel.grid.minor = ggplot2::element_blank(),
        panel.spacing = ggplot2::unit(0.8, "lines"),
        panel.border = ggplot2::element_rect(
          color = "#D0D0D0", 
          fill = NA, 
          linewidth = 0.8
        ),
        plot.title = ggplot2::element_text(
          hjust = 0.5, 
          face = "bold", 
          size = 16,
          margin = ggplot2::margin(b = 10)
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5,
          size = 12,
          color = "#666666",
          margin = ggplot2::margin(b = 15)
        ),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.title = ggplot2::element_text(
          size = 12,
          face = "bold"
        ),
        legend.text = ggplot2::element_text(size = 10),
        legend.spacing.y = ggplot2::unit(0.3, "cm"),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        plot.background = ggplot2::element_rect(fill = "white", color = NA)
      ) +
      ggplot2::labs(
        title = "GO Enrichment Analysis Results",
        subtitle = "Significantly enriched GO terms by module and ontology",
        x = "Module",
        y = "GO Terms"
      )
    
    # 根据数据量调整图形尺寸
    num_modules <- length(unique(top_go$Module))
    plot_width <- max(12, num_modules * 1.5)
    plot_height <- 8 + (length(unique(top_go$Wrapped_Description)) * 0.2)
    
    # 保存PDF（矢量图）
    ggplot2::ggsave("07_GO_Enrichment_Dotplot.pdf", p_go, 
                    width = plot_width, height = plot_height,
                    device = cairo_pdf)
    
    # 保存PNG（高分辨率，根据专家建议）
    ggplot2::ggsave("07_GO_Enrichment_Dotplot.png", p_go, 
                    width = plot_width, height = plot_height, 
                    dpi = 600, bg = "white")
    
    cat("    ✓ 保存GO富集图表 (PDF和PNG格式)\n")
  }
}

# KEGG富集结果可视化
if (exists("all_kegg_results") && nrow(all_kegg_results) > 0) {
  cat("  生成KEGG富集图表...\n")
  
  # 选择每个模块最多8个最显著的通路
  top_kegg <- all_kegg_results %>%
    dplyr::group_by(Module) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 8) %>%
    dplyr::ungroup()
  
  if (nrow(top_kegg) > 0) {
    # 应用文本换行
    top_kegg$Wrapped_Description <- sapply(top_kegg$Description, add_line_breaks)
    
    # 创建颜色梯度
    pvalue_gradient <- colorRampPalette(c("#FDDED7", "#C1E0DB", "#5CB0C3"))(10)
    
    p_kegg <- ggplot2::ggplot(top_kegg, ggplot2::aes(x = Module, y = reorder(Wrapped_Description, -p.adjust))) +
      ggplot2::geom_point(ggplot2::aes(size = Count, fill = -log10(p.adjust)), 
                          shape = 21, color = "gray30", stroke = 1.0, alpha = 0.9) +
      ggplot2::scale_fill_gradientn(
        colors = pvalue_gradient,
        name = expression(-log[10](p.adjust)),
        breaks = scales::pretty_breaks(n = 5)
      ) +
      ggplot2::scale_size_continuous(
        range = c(4, 12), 
        name = "Gene Count",
        breaks = scales::pretty_breaks(n = 5)
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45, 
          hjust = 1, 
          size = 12,
          face = "bold",
          color = "#333333"
        ),
        axis.text.y = ggplot2::element_text(
          size = 11,
          lineheight = 0.9,
          color = "#444444"
        ),
        axis.title.x = ggplot2::element_text(
          size = 14,
          face = "bold",
          margin = ggplot2::margin(t = 10)
        ),
        axis.title.y = ggplot2::element_text(
          size = 14,
          face = "bold",
          margin = ggplot2::margin(r = 10)
        ),
        panel.grid.major = ggplot2::element_line(
          color = "#E8E8E8",
          linewidth = 0.5
        ),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
          color = "#D0D0D0", 
          fill = NA, 
          linewidth = 0.8
        ),
        plot.title = ggplot2::element_text(
          hjust = 0.5, 
          face = "bold", 
          size = 16,
          margin = ggplot2::margin(b = 10)
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5,
          size = 12,
          color = "#666666",
          margin = ggplot2::margin(b = 15)
        ),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.title = ggplot2::element_text(
          size = 12,
          face = "bold"
        ),
        legend.text = ggplot2::element_text(size = 10),
        legend.spacing.y = ggplot2::unit(0.3, "cm"),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
        plot.background = ggplot2::element_rect(fill = "white", color = NA)
      ) +
      ggplot2::labs(
        title = "KEGG Pathway Enrichment Analysis",
        subtitle = "Significantly enriched KEGG pathways by module",
        x = "Module",
        y = "Pathways"
      )
    
    # 根据数据量调整图形尺寸
    num_modules_kegg <- length(unique(top_kegg$Module))
    plot_width_kegg <- max(12, num_modules_kegg * 1.5)
    plot_height_kegg <- 8 + (length(unique(top_kegg$Wrapped_Description)) * 0.2)
    
    # 保存PDF（矢量图）
    ggplot2::ggsave("08_KEGG_Enrichment_Dotplot.pdf", p_kegg, 
                    width = plot_width_kegg, height = plot_height_kegg,
                    device = cairo_pdf)
    
    # 保存PNG（高分辨率，根据专家建议）
    ggplot2::ggsave("08_KEGG_Enrichment_Dotplot.png", p_kegg, 
                    width = plot_width_kegg, height = plot_height_kegg, 
                    dpi = 600, bg = "white")
    
    cat("    ✓ 保存KEGG富集图表 (PDF和PNG格式)\n")
  }
}

# ==============================================================================
# 11. 生成详细分析报告
# ==============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("生成详细分析报告\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 收集分析统计信息
report_info <- list()

report_info$analysis_date <- format(Sys.time(), "%Y年%m月%d日 %H:%M:%S")
report_info$total_genes <- nrow(df)
report_info$total_modules <- length(modules)
report_info$module_list <- paste(modules, collapse = ", ")

# 模块大小分布
report_info$module_sizes <- module_stats

# GO注释提取情况
report_info$go_annotation_count <- ifelse(exists("go_annotation"), nrow(go_annotation), 0)
report_info$go_gene_count <- ifelse(exists("go_annotation"), length(unique(go_annotation$Gene)), 0)

# KEGG注释提取情况
report_info$kegg_annotation_count <- ifelse(exists("kegg_annotation"), nrow(kegg_annotation), 0)
report_info$kegg_gene_count <- ifelse(exists("kegg_annotation"), length(unique(kegg_annotation$Gene)), 0)

# 富集分析结果
report_info$go_enriched_terms <- ifelse(exists("all_go_results") && nrow(all_go_results) > 0, 
                                        nrow(all_go_results), 0)
report_info$kegg_enriched_pathways <- ifelse(exists("all_kegg_results") && nrow(all_kegg_results) > 0, 
                                             nrow(all_kegg_results), 0)

# 生成报告内容
report_lines <- c(
  paste(rep("=", 60), collapse = ""),
  "            GO/KEGG富集分析报告 V3.0",
  paste(rep("=", 60), collapse = ""),
  "",
  paste("分析日期:", report_info$analysis_date),
  "",
  "一、数据概览",
  paste("  总基因数:", report_info$total_genes),
  paste("  模块数量:", report_info$total_modules),
  paste("  模块列表:", report_info$module_list),
  "",
  "二、模块基因数分布",
  "  模块,基因数",
  sapply(1:nrow(module_stats), function(i) {
    paste("  ", module_stats$Module[i], ",", module_stats$Gene_Count[i])
  }),
  "",
  "三、GO注释提取结果",
  paste("  提取的GO注释总数:", report_info$go_annotation_count),
  paste("  有GO注释的基因数:", report_info$go_gene_count),
  paste("  GO注释覆盖率:", round(report_info$go_gene_count/report_info$total_genes*100, 1), "%"),
  if (exists("go_annotation") && nrow(go_annotation) > 0 && "Ontology" %in% colnames(go_annotation)) {
    c("  GO本体分布:",
      sapply(names(table(go_annotation$Ontology)), function(ont) {
        count <- table(go_annotation$Ontology)[ont]
        paste("    ", ont, ":", count)
      }))
  } else {
    "  GO本体分布: 无法获取（GO.db包可能未安装）"
  },
  "",
  "四、KEGG注释提取结果",
  paste("  提取的KEGG注释总数:", report_info$kegg_annotation_count),
  paste("  有KEGG注释的基因数:", report_info$kegg_gene_count),
  paste("  KEGG注释覆盖率:", round(report_info$kegg_gene_count/report_info$total_genes*100, 1), "%"),
  if (report_info$kegg_annotation_count == 0) {
    "  警告: 未提取到KEGG注释"
  },
  "",
  "五、GO富集分析结果",
  if (report_info$go_enriched_terms > 0) {
    c(paste("  显著富集项总数:", report_info$go_enriched_terms),
      paste("  有富集结果的模块数:", 
            ifelse(exists("all_go_results"), length(unique(all_go_results$Module)), 0)),
      "",
      "  分析说明:",
      "    1. 使用p<0.05作为显著性阈值",
      "    2. 使用BH方法进行多重检验校正",
      "    3. 最小基因集大小为5",
      "    4. 仅分析基因数≥10的模块",
      "",
      "  结论: GO富集分析成功，发现了显著富集的生物学功能")
  } else {
    c("  无显著GO富集结果",
      "",
      "  可能原因:",
      "    1. 模块基因数较少",
      "    2. GO注释覆盖度不足",
      "    3. 基因表达变化不显著",
      "    4. p值阈值设置过严",
      "",
      "  建议:",
      "    1. 考虑合并小模块",
      "    2. 检查GO注释质量",
      "    3. 尝试调整p值阈值")
  },
  "",
  "六、KEGG富集分析结果",
  if (report_info$kegg_enriched_pathways > 0) {
    c(paste("  显著富集通路总数:", report_info$kegg_enriched_pathways),
      paste("  有富集结果的模块数:", 
            ifelse(exists("all_kegg_results"), length(unique(all_kegg_results$Module)), 0)),
      "",
      "  结论: KEGG通路富集分析成功，发现了显著富集的代谢通路")
  } else {
    c("  无显著KEGG富集结果",
      "",
      "  可能原因:",
      "    1. KEGG注释覆盖度不足",
      "    2. 模块基因数较少",
      "    3. 未达到显著性阈值",
      "",
      "  建议:",
      "    1. 检查KEGG注释列格式",
      "    2. 确认数据中是否包含有效的KEGG通路信息")
  },
  "",
  "七、可视化参数",
  "  配色方案: 低饱和度色谱 (FDDED7, F5BE8F, C1E0DB, CCD376, A28CC2, 8498AB, 5CB0C3)",
  "  字体大小: 基础字号14pt，坐标轴标签11-12pt",
  "  输出分辨率: PNG 600 DPI (高分辨率)",
  "  输出格式: PDF (矢量图) 和 PNG (位图)",
  "  图表尺寸: 自适应模块数量和GO术语数量",
  "",
  "八、输出文件清单",
  "  01_Module_Gene_Statistics.csv - 模块基因数统计",
  "  02_Data_Quality_Report.csv - 数据质量报告",
  "  03_GO_Annotations_Extracted.csv - 提取的GO注释",
  "  03_GO_Ontology_Distribution.csv - GO本体分布",
  "  04_KEGG_Annotations_Extracted.csv - 提取的KEGG注释",
  if (report_info$go_enriched_terms > 0) {
    "  05_GO_Enrichment_Results.csv - GO富集分析结果"
  },
  if (report_info$kegg_enriched_pathways > 0) {
    "  06_KEGG_Enrichment_Results.csv - KEGG富集分析结果"
  },
  if (report_info$go_enriched_terms > 0) {
    "  07_GO_Enrichment_Dotplot.pdf - GO富集点图 (PDF矢量图)"
  },
  if (report_info$go_enriched_terms > 0) {
    "  07_GO_Enrichment_Dotplot.png - GO富集点图 (PNG高分辨率)"
  },
  if (report_info$kegg_enriched_pathways > 0) {
    "  08_KEGG_Enrichment_Dotplot.pdf - KEGG富集点图 (PDF矢量图)"
  },
  if (report_info$kegg_enriched_pathways > 0) {
    "  08_KEGG_Enrichment_Dotplot.png - KEGG富集点图 (PNG高分辨率)"
  },
  "  09_Analysis_Report.txt - 本分析报告",
  "",
  "九、数据来源说明",
  "  本分析严格基于用户提供的原始数据",
  "  所有注释信息均从数据中的GO和KEGG列提取",
  "  未使用任何外部示例数据或虚构数据",
  "",
  "十、分析参数",
  "  p值阈值: 0.05",
  "  多重检验校正: BH方法",
  "  最小基因集大小: 5",
  "  最大基因集大小: 500",
  "  最小模块基因数: 10",
  "  基因ID清理: 移除版本号（如.1, .2等）",
  "",
  "十一、改进说明（V3.0版本）",
  "  1. 增强KEGG正则表达式，支持更多格式（mtr, map, ko, K等）",
  "  2. 改进GO.db依赖处理，提供更清晰的用户指导",
  "  3. 添加对grey模块的智能处理",
  "  4. 优化可视化自适应性，根据数据量动态调整参数",
  "  5. 增加数据质量报告功能",
  "  6. 改进错误处理和用户交互",
  "",
  paste(rep("=", 60), collapse = ""),
  "                    分析完成",
  paste(rep("=", 60), collapse = "")
)

# 保存报告
writeLines(report_lines, "09_Analysis_Report.txt")
cat("✓ 保存分析报告到: 09_Analysis_Report.txt\n")

# ==============================================================================
# 12. 最终总结和清理
# ==============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("分析完成总结\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

cat("数据来源说明:\n")
cat("  所有分析严格基于您的原始数据\n")
cat("  未使用任何外部示例数据或虚构数据\n")
cat("  所有注释信息均从您的GO和KEGG列提取\n\n")

cat("可视化改进:\n")
cat("  ✓ 使用低饱和度配色方案 (#FDDED7, #F5BE8F, #C1E0DB, #CCD376, #A28CC2, #8498AB, #5CB0C3)\n")
cat("  ✓ 调整字号大小，提高可读性 (基础字号14pt)\n")
cat("  ✓ GO术语自动换行，避免重叠\n")
cat("  ✓ 输出高分辨率PNG (600 DPI) 和PDF矢量图\n")
cat("  ✓ 图表尺寸自适应数据量\n\n")

cat("分析总结:\n")
cat(paste("- 分析基因总数:", report_info$total_genes, "\n"))
cat(paste("- 分析模块数量:", report_info$total_modules, "\n"))
cat(paste("- 提取GO注释数:", report_info$go_annotation_count, "\n"))
cat(paste("- 提取KEGG注释数:", report_info$kegg_annotation_count, "\n"))
cat(paste("- GO富集结果数:", report_info$go_enriched_terms, "\n"))
cat(paste("- KEGG富集结果数:", report_info$kegg_enriched_pathways, "\n"))

cat("\n所有结果已保存到文件夹: ", output_dir, "\n")
cat("请查看 '09_Analysis_Report.txt' 获取详细分析报告\n")

# 显示输出文件
cat("\n生成的文件:\n")
output_files <- list.files()
for (file in output_files) {
  file_size <- file.info(file)$size
  if (!is.na(file_size)) {
    size_kb <- round(file_size / 1024, 2)
    cat(paste("  ", file, " (", size_kb, "KB)\n", sep = ""))
  } else {
    cat(paste("  ", file, "\n", sep = ""))
  }
}

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("V3.0版本改进总结:\n")
cat("  1. 增强KEGG正则表达式灵活性（支持mtr, map, ko, K等格式）\n")
cat("  2. 改进GO.db依赖处理和错误信息\n")
cat("  3. 智能处理grey等特殊模块\n")
cat("  4. 增加数据质量诊断和报告\n")
cat("  5. 优化可视化参数和自适应性\n")
cat("  6. 改进用户交互和进度显示\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 返回原始工作目录
setwd(original_wd)

cat("\n分析脚本执行完成！\n")
cat(paste("所有结果保存在: ", file.path(original_wd, output_dir), "\n"))

# ==============================================================================
# 脚本结束
# ==============================================================================
