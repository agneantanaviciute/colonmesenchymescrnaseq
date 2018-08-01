### Helper functions for Seurat
options(stringsAsFactors = FALSE)

# Modified version of DotPlot returning a dataframe for use with ggplot2
calc_dots <- function(object,genes.plot,cex.use=2,cols.use=NULL,thresh.col=2.5,dot.min=0.05,group.by=NULL,...) {
  if (!(is.null(group.by))) object=SetAllIdent(object,id = group.by)
  object@data=data.frame(t(FetchData(object,genes.plot)))
  
  #this line is in case there is a '-' in the cell name
  colnames(object@data)=object@cell.names
  avg.exp=AverageExpression(object)
  avg.alpha=ClusterAlpha(object)
  cols.use=Seurat:::set.ifnull(cols.use,myPalette(low = "red",high="green"))
  exp.scale=t(scale(t(avg.exp)))
  exp.scale=minmax(exp.scale,max=thresh.col,min=(-1)*thresh.col)
  n.col=length(cols.use)
  data.y=rep(1:ncol(avg.exp),nrow(avg.exp))
  data.x=unlist(lapply(1:nrow(avg.exp),rep,ncol(avg.exp)))
  data.avg=unlist(lapply(1:length(data.y),function(x) exp.scale[data.x[x],data.y[x]]))
  exp.col=(data.avg-min(data.avg))/max(data.avg-min(data.avg))
  data.cex=unlist(lapply(1:length(data.y),function(x) avg.alpha[data.x[x],data.y[x]])) #*cex.use+dot.min

  data.x <- factor(data.x)
  levels(data.x) <- genes.plot
  data.y <- factor(data.y, levels = unique(data.y))
  levels(data.y) <- colnames(avg.alpha)
  ret_df <- data.frame(x = data.x, y = data.y, col = exp.col, cex = data.cex)
  ret_df
}

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  return(genesV2)
}

# Save out gene lists to excel

saveXLS <- function(
  obj,
  filename,
  grp_var = "cluster",
  sort_var = "power",
  sort_thr = 0.4,
  sort_desc = TRUE,
  min_n = 100,
  max_n = Inf
) {
  wb<-createWorkbook(type="xlsx")
  TITLE_STYLE <- CellStyle(wb)+ Font(wb, isBold=TRUE)
  
  if(sort_desc) {
    obj <- obj %>% group_by(get(grp_var)) %>% dplyr::arrange(desc(get(sort_var)), .by_group = TRUE) %>% 
      dplyr::filter(row_number() <= min_n | get(sort_var) > sort_thr & row_number() < max_n) %>% as.data.frame()
  } else {
    obj <- obj %>% group_by(get(grp_var)) %>% dplyr::arrange(get(sort_var), .by_group = TRUE) %>% 
      dplyr::filter(row_number() <= min_n | get(sort_var) < sort_thr & row_number() < max_n) %>% as.data.frame()
  }
  
  grps <- levels(as.factor(obj[,grp_var]))
  
  for(i in 1:length(grps)) {
    df_out <- obj[obj[,grp_var] == grps[i],]
    num_cols <- unlist(lapply(obj,is.numeric))
    df_out[,num_cols] <- signif(df_out[,num_cols], digits = 3)
    sheet <- createSheet(wb, sheetName = gsub("/","_",grps[i]))
    addDataFrame(df_out, sheet, colnamesStyle = TITLE_STYLE, row.names = FALSE)
    autoSizeColumn(sheet, 1:ncol(df_out))
  }
  
  saveWorkbook(wb, filename)
}

# Violin plots

plotGeneViolins <- function(
  geneNames,
  subsets,
  sc,
  plot_type = "violin",
  cex_col = 16,
  cex_row = 16,
  cex_y = 16,
  cex_ylab = 12,
  cex_xlab = 16,
  angle_col = 0,
  angle_row = 0,
  inc_ylab = TRUE,
  inc_legend = TRUE,
  leg_below = FALSE
) {
  # Get data
  subsetNames <- c(paste0(subsets,".H"),paste0(subsets,".D"))
  cell_type <- factor(as.character(sc@ident[sc@ident %in% subsetNames]), levels = subsetNames)
  exp_mat <- as.matrix(sc@data[geneNames,colnames(sc@data)[sc@ident %in% subsetNames]])
  
  cell_class <- cell_type
  levels(cell_class) <- rep(subsets,2)
  
  plot_df <- data.frame(
    cell_type = cell_type,
    cell_class = cell_class,
    condition = factor(ifelse(as.numeric(cell_type) <= length(subsets),"Healthy","DSS"), levels = c("Healthy","DSS")),
    exprs = t(exp_mat), 
    stringsAsFactors = FALSE)
  colnames(plot_df)[4:ncol(plot_df)] <- geneNames
  
  fill_stat <- sapply(as.character(unique(cell_type)), function(x) colMeans(as.matrix(plot_df[cell_type == x,4:ncol(plot_df)])))
  if(!is.matrix(fill_stat)) fill_stat <- t(as.matrix(fill_stat))
  fill_max <- max(t(plot_df[,4:ncol(plot_df)]))
  fill_stat <- (fill_stat/rowMaxs(fill_stat)) * fill_max
  row.names(fill_stat) <- geneNames
  
  plot_df <- reshape2::melt(plot_df, id.vars = c("cell_type","cell_class","condition"))
  plot_df$fill_stat <- sapply(1:nrow(plot_df), function(x) fill_stat[as.character(plot_df$variable[x]),as.character(plot_df$cell_type[x])])
  plot_df$fill <- plot_df$fill_stat
  
  vio_plot <- ggplot(data = plot_df, aes(x=cell_class,y=value))
  
  
  if(plot_type == "violin") vio_plot <- vio_plot + geom_violin(aes(fill=fill), alpha = 0.7, scale = "width", show.legend = inc_legend)
  if(plot_type == "beeswarm") vio_plot <- vio_plot + ggbeeswarm::geom_quasirandom(aes(fill=fill), stroke = 0, size = 0.8, 
                                                                                  alpha = 0.7, pch = 21, bandwidth = 1, show.legend = inc_legend)
  
  vio_plot <- vio_plot +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.5) +
    ylim(c(0,NA)) +
    ylab("Expression level") +
    facet_grid(variable ~ condition, scales = "free_y")
  
  if(inc_ylab) {
    vio_plot <- vio_plot + theme_bw() + 
      theme(axis.text.y=element_text(size = cex_ylab))
  } else {
    vio_plot <- vio_plot + theme_classic() + 
      theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())
  }
  
  vio_plot <- vio_plot +
    theme(
      strip.text.x = element_text(size = cex_col, angle = angle_col),
      strip.text.y = element_text(size = cex_row, angle = angle_row),
      axis.text.x=element_text(size = cex_xlab),
      strip.background = element_blank(),
      axis.line.y=element_blank(),
      axis.title.x=element_blank(), 
      axis.title.y=element_text(size = cex_y, margin = margin(0,30,0,0)),
      legend.position = ifelse(leg_below, "bottom", "right")
    )
  
  vio_plot + 
    scale_fill_viridis(
      "", begin = 0.2, end = 0.8, option = "inferno", 
      limits = c(0,fill_max), breaks = c(0,fill_max),
      labels = c(0,"Maximum")
    )
}