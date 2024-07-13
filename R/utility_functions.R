
generate_HighExpressionFeaturesFromModel <- function(seuratObj,
 tumor_markers, k)
{
  data("IGFBP5_mclust_model")
  data("IGFBP2_mclust_model")
  data("PDGFRA_mclust_model")
  data("EGFR_mclust_model")
  data("CDK4_mclust_model")
  thresholds <- rep(Inf, length(tumor_markers))
  names(thresholds) <- tumor_markers
    data <- seuratObj@assays$RNA@data[match(tumor_markers, rownames(seuratObj@assays$RNA@data)), ]
  for (i in 1:length(tumor_markers))
  {
    model <- NULL
    human_inhouse_data <- seuratObj@assays$RNA@data[match(tumor_markers[i], rownames(seuratObj@assays$RNA@data)), ]
    c("" ,""  , "CDK4"  , "" ,"SOX2"   , "")
    if(sum(human_inhouse_data)>0){
        d <- as.numeric(human_inhouse_data)
        mod <- Mclust(na.omit(d), verbose = F, G=k)
        thresholds[i] <- min(d[which(mod$z[,k]>=0.99)])
        if(tumor_markers[i]=="EGFR")  model <- EGFR_mclust_model
        if(tumor_markers[i]=="IGFBP2")   model <- IGFBP2_mclust_model 
        if(tumor_markers[i]=="PDGFRA")   model <- PDGFRA_mclust_model 
        if(tumor_markers[i]=="IGFBP5")   model <- IGFBP5_mclust_model 
        if(tumor_markers[i]=="CDK4")   model <- CDK4_mclust_model 
    
        if(!is.null(model))
        {
            predTrainNew <- predict(model,human_inhouse_data)

            if(thresholds[i] < min(human_inhouse_data[predTrainNew$z[,2]==1]))
            {   
                thresholds[i] <- min(human_inhouse_data[predTrainNew$z[,2]==1])
            }
        } 
      }
  }

  if(!is.null(thresholds)){
    if(!all(names(thresholds)==tumor_markers))
    {
       stop("gene_thr names are not the same as genes...") 
    }
  }

  tumor_markers_mat <- matrix(0, nrow=dim(data)[1], ncol=dim(data)[2])
  rownames(tumor_markers_mat) <- rownames(data)
  colnames(tumor_markers_mat) <- colnames(data)

  for (i in 1:length(thresholds))
  {
    tumor_markers_mat[i,data[i,]>thresholds[[i]]] <- 1
  } 

  return(list(tumor_markers=tumor_markers_mat, gene_thr=thresholds))
}

writeSCRAMresults <- function(object,  project)
{
  cell_type_level=unique(bind_rows(object@setsAnnotatedCellTypeLevel, .id = "cluster_id"))
  results_maximal <- list(cell_type_level=unique(cell_type_level[cell_type_level$is.maximal, ])) 
  openxlsx::write.xlsx(results_maximal,file = paste0( project, ".SCRAM.MAXIMAL.xlsx"), overwrite = T)
   cell_type_level=unique(bind_rows(object@setsAnnotatedCellTypeLevel, .id = "cluster_id"))
  results_maximal <- list(cell_type_level=unique(cell_type_level[, ])) 
  openxlsx::write.xlsx(results_maximal,file = paste0( project, ".SCRAM.ALL.xlsx"), overwrite = T)

}



readCASPER<- function(finalChrMat, seuratObj, amp_chr, del_chr, min_threshold, plot=T, project)
{

  colnames(finalChrMat) <- paste0("chr", colnames(finalChrMat))
  if(!is.null(seuratObj)) 
  {
      if(!all(colnames(seuratObj)==rownames(finalChrMat))) {
          stop ("Casper object colnames and seurat_obj are not the same..")
      }
  }

  if(is.null(min_threshold))   stop ("min_threshold: Minimum number of cells of CNV event is missing")
  amp <- intersect(names(which(apply(finalChrMat, 2, function(x) length(which(x>0)))>min_threshold)), 
         amp_chr)
  
  del <- intersect(names(which(apply(finalChrMat, 2, function(x) length(which(x<0)))>min_threshold)), 
         del_chr)

  cnv <- NULL
  
  if(length(del)>0) 
   {
      meta <- finalChrMat[,del, drop=F]
      meta[meta>0] <- 0
      meta[meta<0] <- 1
      if(!is.null(seuratObj) & plot==T)
      { 
          seuratObj@meta.data<- cbind(seuratObj@meta.data, meta)
          p1 <- DimPlot(seuratObj, group.by=del, cols =  c("grey", "blue"), order=T)
          ggsave(paste0(project, "_CNV_del.pdf"), p1, width=length(del)*5, height=5)
          ggsave(paste0(project, "_CNV_del.png"), p1, width=length(del)*5, height=5)
       }
      
       cnv <- data.frame(meta)
       colnames(cnv) <- c(paste0(del, "_del"))     
   }
   
   if(length(amp)>0) 
   {
      meta <- finalChrMat[,amp, drop=F]
      meta[meta<0] <- 0
      meta[meta>0] <- 1
      if(!is.null(seuratObj) & plot==T)
      { 
        seuratObj@meta.data<- cbind(seuratObj@meta.data, meta)
        p1 <- DimPlot(seuratObj, group.by=amp, cols =  c( "grey", "red"), order=T)
        ggsave(paste0(project ,"_CNV_amp.pdf"), p1, width=length(amp)*5, height=5)
        ggsave(paste0(project ,"_CNV_amp.png"), p1, width=length(amp)*5, height=5)
    }
      cnv2 <-  data.frame(meta)
      colnames(cnv2) <- c(paste0(amp, "_amp"))
      if(!is.null(cnv)) cnv <- cbind(cnv, cnv2)
      if(is.null(cnv)) cnv <- cnv2
     
   }
   if(!is.null(cnv)) cnv <- t(cnv)
   return(cnv)
}


readCASPER_patientSpec<- function(finalChrMat,finalChrMatbulk, seuratObj, plot=T, project)
{

  colnames(finalChrMat_bulk) <- paste0("chr", colnames(finalChrMat_bulk))
  del <- apply(finalChrMat_bulk, 1, function(x) names(which(x<0)))
  amp <- apply(finalChrMat_bulk, 1, function(x) names(which(x>0)))

  if(!all(names(amp)==names(del))) {
          stop ("amp and del names are not same")
  }

  colnames(finalChrMat) <- paste0("chr", colnames(finalChrMat))
  if(!is.null(seuratObj)) 
  {
      if(!all(colnames(seuratObj)==rownames(finalChrMat))) {
          stop ("Casper object colnames and seuratObj are not the same..")
      }
  }

  cnv_amp  <- matrix(0, nrow=nrow(finalChrMat), ncol=ncol(finalChrMat))
  cnv_del  <- matrix(0, nrow=nrow(finalChrMat), ncol=ncol(finalChrMat))
  colnames(cnv_amp) <- colnames(finalChrMat)
  colnames(cnv_del) <- colnames(finalChrMat)

  rownames(cnv_amp) <- rownames(finalChrMat)
  rownames(cnv_del) <- rownames(finalChrMat)
  finalChrMat <- as.matrix(finalChrMat)
  itr <- max(length(amp), length(del))

  for(i in 1:itr)
  {   
     if(length(amp)>0){
      if(length(amp[[i]])>0)
      {
        cells <- colnames(seuratObj)[seuratObj$orig.ident==names(amp)[i]]

        if(length(cells)>0) cnv_amp[match(cells, rownames(finalChrMat)),amp[[i]]] <- finalChrMat[match(cells, rownames(finalChrMat)),amp[[i]]]
      }
    }
    if(length(del)>0){
      if(length(del[[i]])>0)
      {
        cells <- colnames(seuratObj)[seuratObj$orig.ident==names(del)[i]]

         if(length(cells)>0) cnv_del[match(cells, rownames(finalChrMat)),del[[i]]] <- finalChrMat[match(cells, rownames(finalChrMat)),del[[i]]]
      }
    }


  }


  cnv_del[cnv_del>0] <- 0
  cnv_del[cnv_del<0] <- 1

  cnv_amp[cnv_amp<0] <- 0
  cnv_amp[cnv_amp>0] <- 1

 
  colnames(cnv_del) <- c(paste0( colnames(cnv_del), "_del"))
  colnames(cnv_amp) <- c(paste0( colnames(cnv_amp), "_amp"))
  
  del_chr <- which(apply(cnv_del, 2, sum)>0)
  amp_chr <- which(apply(cnv_amp, 2, sum)>0)

  if(length(del_chr)>0) 
  {
    cnv_del <- cnv_del[, del_chr, drop=F]
    seuratObj@meta.data<- cbind(seuratObj@meta.data, cnv_del)
    p1 <- DimPlot(seuratObj, group.by=names(del_chr), cols =  c("grey", "blue"), order=T, raster=F)
    ggsave(paste0(project, "_CNV_del.pdf"), p1, width=length(del_chr)/2*5, height=length(del_chr)/2*5)
    ggsave(paste0(project, "_CNV_del.png"), p1, width=length(del_chr)/2*5, height=length(del_chr)/2*5)
  }

  if(length(amp_chr)>0) {
    cnv_amp <- cnv_amp[, amp_chr, drop=F]
    seuratObj@meta.data<- cbind(seuratObj@meta.data, cnv_amp)
    p1 <- DimPlot(seuratObj, group.by=names(amp_chr), cols =  c("grey", "red"), order=T, raster=F)
    ggsave(paste0(project, "_CNV_amp.pdf"), p1, width=length(amp_chr)/2*5, height=length(amp_chr)/2*5)
    ggsave(paste0(project, "_CNV_amp.png"), p1, width=length(amp_chr)/2*5, height=length(amp_chr)/2*5)
  }
 
 
  cnv <- NULL
  if(!is.null(del_chr) & !is.null(amp_chr)) cnv <- cbind(cnv_del, cnv_amp)
  if(is.null(del_chr)) cnv <- cnv_amp
  if(is.null(amp_chr)) cnv <- cnv_del  
  if(!is.null(cnv)) cnv <- t(cnv)
   return(cnv)
  
}


generate_HighExpressionFeatures <- function(seuratObj, genes, k=3, prob=0.99, gene_thr=NULL)
{

  data <- seuratObj@assays$RNA[match(genes, rownames(seuratObj@assays$RNA)), ]
  
  if(is.null(gene_thr)){
    gene_thr <- apply(data, 1, function(x){
        d <- as.numeric(x)
       # plot(density(d))
        mod <- Mclust(na.omit(d), verbose = F, G=k)
        thr <- min(d[which(mod$z[,k]>=prob)])
        return(thr)
    })
  }

  if(!is.null(gene_thr)){
    if(!all(names(gene_thr)==genes))
    {
       stop("gene_thr names are not the same as genes...") 
    }
  }

  tumor_markers <- matrix(0, nrow=dim(data)[1], ncol=dim(data)[2])
  rownames(tumor_markers) <- rownames(data)
  colnames(tumor_markers) <- colnames(data)

  for (i in 1:length(gene_thr))
  {
      tumor_markers[i,data[i,]>gene_thr[[i]]] <- 1
  } 

  return(list(tumor_markers=tumor_markers, gene_thr=gene_thr))

}


readXCVATROutput <- function(XCVATR_folderPath, name_mapping)
{
  samples <- name_mapping$filenames
  samples_seuratObj <- name_mapping$objectNames
  seuratObj@meta.data$AF <- 0

  all_cosmic_variants <- c()
  for (i in 1:length(samples))
  {

          counts <- read.delim(paste0(XCVATR_folderPath, "/", samples[i], "/ALL_IMPACTFUL_RARE_ALLELE_COUNTS.txt"))
          print(samples[i])
          print(dim(counts)[2])
      
          all_genes <- sapply(counts$REF_ALT, function(x) strsplit(x, split=" ")[[1]][4])
                
          cosmic_snvs <- read.delim(paste0(XCVATR_folderPath, "/", samples[i],"/cosmic_selected_snvs.txt"), header=F)
          cosmic_indels <- read.delim(paste0(XCVATR_folderPath, "/", samples[i], "/cosmic_selected_indels.txt"), header=F)
          
          cov <- read.delim(paste0(XCVATR_folderPath, "/", samples[i], "/pileup_snvs.op"), header=F)
          cov_all <- paste(cov[ ,1],cov[ ,2]-1,cov[ ,2], sep=":")
          cov_ratio <- cov[ ,5]/(cov[ ,5]+cov[ ,6])
          cov_all <- cov_all[which(cov_ratio>=0.1)]
      
          snvs_all1 <- paste(counts[ ,1],counts[ ,2],counts[ ,3],counts[ ,4], sep=":")
        
          snvs_all <- paste(counts[ ,1],counts[ ,2],counts[ ,3], sep=":")
          cosmic1 <- paste(cosmic_snvs[ ,1],cosmic_snvs[ ,2],cosmic_snvs[ ,3], sep=":")
          cosmic2 <- paste(cosmic_indels[ ,1],cosmic_indels[ ,2],cosmic_indels[ ,3], sep=":")
          f <- intersect(cov_all, unique(c(cosmic1, cosmic2)))
          cosmic_variants <- paste0(all_genes,sep=":", snvs_all1)[snvs_all %in% f]
          all_cosmic_variants <- c(all_cosmic_variants, cosmic_variants)
  }

  all_cosmic_variants <- unique(all_cosmic_variants)

  var_matrix <- matrix(0, nrow= length(all_cosmic_variants), ncol=length(colnames(seuratObj)))
  rownames(var_matrix) <- all_cosmic_variants

  for (i in 1:length(samples))
  {
    counts <- read.delim(paste0(XCVATR_folderPath, "/", samples[i], "/ALL_IMPACTFUL_RARE_ALLELE_COUNTS.txt"))
    all_genes <- sapply(counts$REF_ALT, function(x) strsplit(x, split=" ")[[1]][4])
    cosmic_snvs <- read.delim(paste0(XCVATR_folderPath, "/", samples[i],"/cosmic_selected_snvs.txt"), header=F)
    cosmic_indels <- read.delim(paste0(XCVATR_folderPath, "/", samples[i], "/cosmic_selected_indels.txt"), header=F)
    snvs_all <- paste(counts[ ,1],counts[ ,2],counts[ ,3], sep=":")
    cosmic1 <- paste(cosmic_snvs[ ,1],cosmic_snvs[ ,2],cosmic_snvs[ ,3], sep=":")
    cosmic2 <- paste(cosmic_indels[ ,1],cosmic_indels[ ,2],cosmic_indels[ ,3], sep=":")
    snvs_all <- paste0(all_genes,sep=":", snvs_all)
    counts2 <- counts[snvs_all %in% all_cosmic_variants,]
    genes <- sapply(counts2$REF_ALT, function(x) strsplit(x, split=" ")[[1]][4])
    vars <- snvs_all[snvs_all%in% all_cosmic_variants]

    counts3 <- counts2
    rownames(counts3) <- vars

    if(dim(counts3)[1]>0){
      AF_l <- apply(counts3 [,5:dim(counts3)[2]], 1,function(x) { toks<-strsplit(x, split=' ');names(which(unlist(lapply(toks, function(x) as.numeric(x[2])>0))))}) 
          if(length(AF_l)>0)
          {
            for (j in 1:length(AF_l))
              {
                  AF <- AF_l[[j]]
      
                  if(length(AF)>0){        
                      cellIds <- unique(gsub("\\.", "-", unlist(AF)))
                      ext <- unique(unlist(lapply(strsplit(rownames(seuratObj@meta.data)[samples_seuratObj[i]==seuratObj$orig.ident],split="_"), function(x) x[[2]])))
                      common <- intersect(paste0(cellIds, "_", ext), rownames(seuratObj@meta.data))
                      print(paste0(samples[i], ":", length(common), ":ext:", ext))
                
                      if(length(common)>0) {
                        var_matrix[match(names(AF_l)[j],rownames(var_matrix)), match(common, rownames(seuratObj@meta.data))] <- 1
                      }
                  }
                  
              }
              
          }

      }
}
  return(var_matrix)
}

# generate_pseudobulkPatientForCasper <- function(seuratObj, controlCellID="control")
# {
#   control_cells <- which(seuratObj$tumorType==controlID)
#   cells <- colnames(seuratObj)[unique(c(control_cells, grep("HIGH", seuratObj$celltype)))]
#   sub <- subset(seuratObj, cells=cells)
#   data <- (sub@assays$RNA@counts)
#   samples  <- names(which(table(sub$orig.ident)>10))
#   data2 <- lapply(1:length(samples), function(x) rowSums(data[, which(sub$orig.ident==samples[x]), na.rm = T]))
#   raw.data <- do.call(cbind, data2)
#   colnames(raw.data) <- samples
#   raw.data<- t(apply(raw.data, 1, function(x) as.numeric(x)))
#   colnames(raw.data) <- samples
#   return(raw.data)

# }