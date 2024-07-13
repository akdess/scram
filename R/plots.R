

plotUMAPCellTypePerCluster <- function(object, project, cluster_id) {
    seuratObj <- object@seurat_obj
    seuratObj$celltypes <- object@cellTypeLevelAnnotationDetailed
    seuratObj$celltypes[seuratObj$celltypes == ""] <- "unknown"
    clusters <- sort(as.numeric(as.character(unique(seuratObj$seurat_clusters))))


    for (i in 1:length(clusters)) {
        cls <- clusters[i]
        srat <- subset(seuratObj, idents = (cls))
        if (dim(srat)[2] > 5) {
            threshold <- (length(colnames(srat)) * 3 / 100)
            if (threshold < 30) threshold <- 30
            srat$celltypes[srat$celltypes %in% names(which(table(srat$celltypes) < threshold))] <- "others"

            p1 <- DimPlot(srat, group.by = "celltypes") + theme(legend.position = "none")
            data <- data.frame(table(srat$celltypes) / sum(table(srat$celltypes)))
            colnames(data) <- c("group", "value")

            # Basic piechart
            p2 <- ggplot(data, aes(x = "", y = value, fill = group)) +
                geom_bar(stat = "identity", width = 1, color = "white") +
                coord_polar("y", start = 0) +
                theme_void() +
                guides(fill = guide_legend(title = NULL))
            ggsave(paste0(project, "_cluster", cls, "_plotUMAPCellType.pdf"), p1 + p2, width = 10, height = 20)
            ggsave(paste0(project, "_cluster", cls, "_plotUMAPCellType.png"), p1 + p2, width = 10, height = 20)
        }
    }
}

plotCasperBulkSignalHeatmap <- function(object, project) {

    casper_object <- object@casper_bulk

    load(paste0("./cnv_casper/", project, "_BULK_finalChrMat_thr_3.rda"))

    data <- casper_object@control.normalized[[3]]
    data <- na.omit(data)
    x.center <- mean(data)
    quantiles = quantile(data[data != x.center], c(0.01, 0.99))
    delta = max(abs(c(x.center - quantiles[1], quantiles[2] - 
        x.center)))
    low_threshold = x.center - delta
    high_threshold = x.center + delta
    x.range = c(low_threshold, high_threshold)
    data[data < low_threshold] <- low_threshold
    data[data > high_threshold] <- high_threshold

    breaks <- seq(-1, 1, length = 16)
    color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks))
    idx <- cumsum(table(casper_object@annotation.filt$Chr)[as.character(1:22)])
    xlabel <- rep("", length(rownames(casper_object@data)))
    half <- round(table(casper_object@annotation.filt$Chr)[as.character(1:22)]/2)[-1]
    xpos <- c(half[1], (idx[-22] + half))
    xlabel[xpos] <- 1:22
    xlabel[which(is.na(xlabel))] <- ""

    thr <- 1

    load(paste0("./cnv_casper/", project, "_BULK_finalChrMat_thr_", thr, ".rda"))

    p<- pheatmap(t(finalChrMat),cluster_cols = F, cluster_rows = T, fontsize_row=6,
                            color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                    show_rownames = T)


    p2<- pheatmap(t(data[ ,p$tree_row$order]), cluster_cols = F, cluster_rows = T, gaps_col = idx[1:21], 
            labels_col = xlabel,color = color, breaks = breaks, fontsize_row=6,
                    show_rownames = T)

    plot_list <- list()
    plot_list[[1]] = p2[[4]]     ##to save each plot into a list. note the [[4]]
    plot_list[[2]] = p[[4]]     ##to save each plot into a list. note the [[4]]

    g<-do.call(grid.arrange,plot_list)

    ggsave(paste0(project, ".heatmap.png"), g, width=20, height=10)
    ggsave(paste0(project, "cnv.signal.heatmap.png"), p2, width=20, height=5)

}

plotModels <- function(seuratObj,results_list, outdir="original_UMAP_model_plots" )
{

  temp <- seuratObj@meta.data 
  dir.create(outdir)
  for (i in 1:length(results_list))
  {
      dat <- results_list[[i]]
      all(colnames(seuratObj)== dat$cells)
      feats <- unique(dat[dat[,grep("_clusterScore", colnames(dat))]>0.9,grep("_cluster$", colnames(dat))])
       seuratObj@meta.data <- temp
      seuratObj@meta.data <- cbind(seuratObj@meta.data, dat)
 
      feats <- gsub(" ", ".", feats)
      feats <- gsub("-", ".", feats)
      feats <- gsub("/", ".", feats)
      feats <- gsub("\\+", ".", feats)
      feats <- gsub("\\&", ".", feats)
      feats <- gsub("\\(", ".", feats)
      feats <- gsub("\\)", ".", feats)
      feats <- na.omit(feats)

      feats_toplt <- unlist(sapply(feats, function(x) grep(x, colnames(seuratObj@meta.data ))))
      print(which(unlist(lapply(feats_toplt, length)==0)))
      p1 <- FeaturePlot(seuratObj, features = colnames(seuratObj@meta.data )[feats_toplt], min.cutoff=0.9, order=T, label=T)
      ggsave(paste0(outdir, "/", names(results_list)[i], ".png"), p1, width=20, height=25)
      p1 <- FeaturePlot(seuratObj, features = colnames(seuratObj@meta.data )[feats_toplt], min.cutoff=0.9, order=F, label=T)
      ggsave(paste0(outdir, "/", names(results_list)[i], "_orderFalse.png"), p1, width=20, height=25)
  }
}

plotPercAnnotatedCells <- function(object, project) {
    seuratObj <- object@seurat_obj
    seuratObj$cellclass2 <- rep("annotated", length(scram_obj@cellTypeLevelAnnotationDetailed))
    seuratObj$cellclass2[object@cellTypeLevelAnnotationDetailed == ""] <- "unknown"

    srat <- seuratObj
    data <- data.frame(table(srat$cellclass2) / sum(table(srat$cellclass2)))
    colnames(data) <- c("group", "value")

    # Basic piechart
    p2 <- ggplot(data, aes(x = "", y = value, fill = group)) +
        geom_bar(stat = "identity", width = 1, color = "white") +
        coord_polar("y", start = 0) +
        theme_void()
    ggsave(paste0(project, "_percAnnotatedCells.pdf"), p2, width = 5, height = 5)
    ggsave(paste0(project, "_percAnnotatedCells.png"), p2, width = 5, height = 5)
}
