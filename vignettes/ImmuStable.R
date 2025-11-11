## ----eval=FALSE---------------------------------------------------------------
# library(ImmuStable)
# 
# # 加载参考数据库
# data(WT_database_data)
# 
# # 计算 Z-score
# zmat <- compute_zscore(seurat_obj, celltype_col = "celltype", wt_ref = WT_database_data)
# 
# # 富集分析
# res <- enrich_by_celltype(zmat, out_prefix = "example")
# 
# # 绘制结果
# plots <- enrich_results_plot(res, celltype = "B cells", direction = "up", which = c("KEGG","BP"))
# plots$KEGG$dot
# plots$BP$bar

