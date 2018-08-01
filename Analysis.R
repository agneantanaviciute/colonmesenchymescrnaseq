# Murine 10X genomics experiment - Processing pipeline
rmarkdown::render("cell_detection.Rmd", output_file = "reports/cell_detection_hc.html", 
                  params = list(source_dir = 'data/HC_3L/raw_gene_bc_matrices_mex/mm10_pL63_mCit', 
                                cells_expected = 4500,
                                min_UMI = 2500, 
                                proj_name = '10X_Ms_Healthy',
                                ctrl_cell_clust = 7,
                                batch_prefix = 'HC'))

rmarkdown::render("cell_detection.Rmd", output_file = "reports/cell_detection_dss.html", 
                  params = list(source_dir = 'data/DSS_3L/raw_gene_bc_matrices_mex/mm10_pL63_mCit', 
                                cells_expected = 4500,
                                min_UMI = 3500, 
                                proj_name = '10X_Ms_DSS',
                                ctrl_cell_clust = 9,
                                batch_prefix = 'DSS'))

rmarkdown::render("norm_clust.Rmd", output_file = "reports/norm_clust_hc.html", 
                  params = list(source_dir = 'data/HC_3L/raw_gene_bc_matrices_mex/mm10_pL63_mCit', 
                                source_metadata = 'output/celldetection/10X_Ms_Healthy.RData',
                                proj_name = '10X_Ms_Healthy',
                                calc_jackstraw = 'no',
                                pcs_use = 28,
                                cluster_res = 1))
rmarkdown::render("clust_analysis.Rmd", output_file = "reports/clust_analysis_hc.html", 
                  params = list(proj_name = '10X_Ms_Healthy'))

rmarkdown::render("norm_clust.Rmd", output_file = "reports/norm_clust_dss.html", 
                  params = list(script_name = 'norm_clust',
                                source_dir = 'data/DSS_3L/raw_gene_bc_matrices_mex/mm10_pL63_mCit', 
                                source_metadata = 'output/celldetection/10X_Ms_DSS.RData',
                                proj_name = '10X_Ms_DSS',
                                calc_jackstraw = 'no',
                                pcs_use = 28,
                                cluster_res = 0.6,
                                tsne_pcs = 28))
rmarkdown::render("clust_analysis.Rmd", output_file = "reports/clust_analysis_dss.html", 
                  params = list(proj_name = '10X_Ms_DSS'))

rmarkdown::render("cluster_biology.Rmd", output_file = "reports/cluster_biology.html")

rmarkdown::render("ont_DE.Rmd", output_file = "reports/ont_DE.html")

rmarkdown::render("DPT.Rmd", output_file = "reports/DPT.html")

rmarkdown::render("ctrl_cells.Rmd", output_file = "reports/ctrl_cells.html")
