
CreateSCRAMObject <- function(seurat_obj,
                              organism,
                              min_support, max_set_size) {
    object <- new(
        Class = "scramv2", seurat_obj = seurat_obj, organism = organism, min_support = min_support, max_set_size = max_set_size
    )

    return(object)
}


readModelsPretrained <- function(path, prob_thr, refs, run)
{

    results_list <- list()

    for (i in 1:length(refs)) {
    # Try to read the CSV file and process it
    tryCatch({
        dat <- read.csv(paste0(path, "/", run, "_", refs[i], "_neuralNetwork_scores.csv"))

        # Check for the column with "_clusterScore"
        cluster_score_col <- grep("_clusterScore", colnames(dat))
        if (length(cluster_score_col) == 0) {
        stop(paste("No '_clusterScore' column found in", refs[i], "_neuralNetwork_scores.csv"))
        }

        unknown <- which(dat[, cluster_score_col] < prob_thr)

        # Check for the column with "_cluster$"
        cluster_col <- grep("_cluster$", colnames(dat))
        if (length(cluster_col) == 0) {
        stop(paste("No '_cluster' column found in", refs[i], "_neuralNetwork_scores.csv"))
        }

        dat[unknown, cluster_col] <- ""
        results_list[[i]] <- dat
        names(results_list)[i] <- refs[i]

    }, error = function(e) {
        # If an error occurs, print a message and continue with the next iteration
        message(paste("Error processing", refs[i], ": ", e$message))
        results_list[[i]] <- NULL  # Optionally store NULL or some other placeholder in case of error
        names(results_list)[i] <- refs[i]
    })
    }
  return (results_list)
}

readMODELS <- function(path, prob_thr)
{
  files <- files <- list.files(path)
  results_list <- list()

  for (i in 1:length(files))
  {
      dat <- read.csv(paste0(path, files[i], "/",files[i], "_scores.csv"))
      unknown <- which(dat[,grep("_clusterScore", colnames(dat))]<prob_thr)
      dat[unknown,grep("_cluster$", colnames(dat))] <- ""
      results_list[[i]] <- dat
      names(results_list)[i] <- files[i]
  }

  return (results_list)
}

f3 <- function(vec) {
  U <- sort(unique(vec))
  M <- matrix(0, nrow = length(vec),
              ncol = length(U),
              dimnames = list(NULL, U))
  M[cbind(seq_len(length(vec)), match(vec, U))] <- 1L
  M
}

createModelMatrix <- function(results_list,seuratObj){
  r <- list()
  for (i in 1:length(results_list))
  {
    r [[i]] <- f3(as.vector(results_list[[i]][,grep("_cluster$", colnames(results_list[[i]]))]))
  }
  models <- do.call(cbind, r)
  rownames(models) <- colnames(seuratObj)
  models
}

createCellTypeMatrix <- function(object, nn_path, prob_thr, refs, run, pretrained=T)
{
    data("celltype_mapping")
    if (pretrained) results_list <- readModelsPretrained(path=nn_path, prob_thr=prob_thr, refs=refs, run=run)
    else {  results_list <- readMODELS(path=nn_path, prob_thr=prob_thr) }
    seuratObj <- object@seurat_obj
    models <- createModelMatrix(results_list,seuratObj)

    models <- t(models)
    models <- models[rownames(models)!="", ]
    models <- models[rownames(models)!="NA_suva", ]
    #models <- models[rownames(models)!="immune_suva", ]
    models <- models[rownames(models)!="Oligo", ]
    object@models <- models
    object@results_list <- results_list

    rownames(models) <- ct[match(rownames(models),ct$celltype),"celltype_simplified"]
    simple_cell_type <- unique(rownames(models))
    simplified_models <- sapply(simple_cell_type, function(x) apply(models[rownames(models) %in% x,, drop=F ], 2, sum))
    simplified_models <- t(simplified_models)

    simplified_models["tumor",] <- simplified_models["tumor",]
    simplified_models["immune",simplified_models["immune",]<3] <- 0
    simplified_models["immune",simplified_models["immune",]>2] <- 1
 #   simplified_models <- simplified_models[is.na(rownames(simplified_models)),  ]

    simplified_models["immune",] <-  simplified_models["immune",] + apply(object@models[rownames(object@models)[grep("primary_gbm", rownames(object@models))], ], 2, sum)
    simplified_models["immune",simplified_models["immune",]>1] <- 1

    object@simplified_models <- simplified_models
    object
}



runCASPER_Bulk_Scell<- function (object, sampleCol="orig.ident", project){

    simplified_models <- object@simplified_models
    seuratObj <- object@seurat_obj
    simplified_models["tumor",] <- simplified_models["tumor",] +simplified_models["cycling progenitor cell",]

    possible_tumors  <- which(simplified_models["tumor",] >0)
    control_cells <- which(simplified_models["immune",] >0)

    seuratObj@meta.data[control_cells, sampleCol] <- "control"
    cells <- colnames(seuratObj)[unique(c(control_cells, possible_tumors))]
    sub <- subset(seuratObj, cells=cells)
    data <- GetAssayData(sub, assay='RNA', layer='counts')
    samples  <- names(which(table(sub@meta.data[, sampleCol])>10))
    data2 <- lapply(1:length(samples), function(x) apply(data[, which(sub@meta.data[, sampleCol]==samples[x]), drop = FALSE], 1, sum))
    raw.data <- do.call(cbind, data2)
    colnames(raw.data) <- samples
    raw.data<- t(apply(raw.data, 1, function(x) as.numeric(x)))
    colnames(raw.data) <- samples

    annotation <- generateAnnotation(id_type="hgnc_symbol", genes=rownames(raw.data), ishg19=F)
    raw.data <- raw.data[match( annotation$Gene,rownames(raw.data)), ]

    cytoband <- cytoband_hg38
    casper_object <- CreateCasperObject(raw.data=raw.data,loh.name.mapping=NULL,
        sequencing.type="bulk",
        cnv.scale=3, loh.scale=3, window.length=50, length.iterations=50,
        expr.cutoff = 1, filter="mean", matrix.type="raw",
        annotation=annotation, method="iterative", loh=NULL,
        control.sample.ids="control", cytoband=cytoband)

    dir.create("cnv_casper")
    casper_object<- runCaSpERWithoutLOH(casper_object, project=paste0("./cnv_casper/",project, "_BULK"))
    object@casper_bulk <- casper_object

    DefaultAssay(seuratObj) <- "RNA"
    data <-  GetAssayData(seuratObj, assay='RNA', layer='counts')
    ## cells with seuratObj$tumorType=="Normal" are used as control
    annotation <- generateAnnotation(id_type="hgnc_symbol", genes=rownames(data), ishg19=F)
    controls <- colnames(seuratObj)[control_cells]

    raw.data <- data[match( annotation$Gene,rownames(data)), ]
    cytoband <- cytoband_hg38

    print("running single cell CNV calling step, might take same time...")
    casper_object <- CreateCasperObject(raw.data=as.matrix(raw.data),loh.name.mapping=NULL,
    sequencing.type="single-cell",
    cnv.scale=3, loh.scale=3, window.length=50, length.iterations=50,
    expr.cutoff = 0.1, filter="mean", matrix.type="raw",
    annotation=annotation, method="iterative", loh=NULL,
    control.sample.ids=controls, cytoband=cytoband)

    casper_object<- runCaSpERWithoutLOH(casper_object, project=paste0("./cnv_casper/",project, "_SCELL"))
    object@casper_scell <- casper_object
    object

}

runFinalTumorAnnotation <- function(object,  loadCasper=T, loadNumbat=T, finalChrMat_bulk, loadXCVATR=T, sampleCol="orig.ident", project, bulk_thr=1, scell_thr=1, model_genes=c("PDGFRA" ,"EGFR"  ,"SOX2" ))
{
    data("celltype_mapping")
    seuratObj <- object@seurat_obj
    simplified_models <- object@simplified_models
    object@cnv  <- NULL
    object@snv <- NULL
    object@tumor <- NULL
    seuratObj$snv_num <- 0
    seuratObj$cnv_num <- 0

    if(loadCasper) {
        if(is.null(finalChrMat_bulk)){
            bulk_thr <- 1
            load(paste0("./cnv_casper/", project, "_BULK_finalChrMat_thr_", bulk_thr, ".rda"))
            finalChrMat_bulk <- t(finalChrMat)
        }
        scell_thr <- 3
        load(paste0("./cnv_casper/", project, "_SCELL_finalChrMat_thr_", scell_thr, ".rda"))
        finalChrMat <- t(finalChrMat )

        common <- intersect(rownames(finalChrMat), colnames(seuratObj))
        finalChrMat <- finalChrMat[ match(common, rownames(finalChrMat)) , ,drop=F]
        finalChrMat <- finalChrMat[ match(colnames(seuratObj), rownames(finalChrMat)) , ,drop=F]
        if(!all((rownames(finalChrMat)==colnames(seuratObj)))) stop("CNV sample names and SeuratObject sample names are different")

        if(all((rownames(finalChrMat)==colnames(seuratObj))))
        {
             cnv <- readCASPER_patientSpec(finalChrMat,finalChrMatbulk, seuratObj, plot=T, project)
        }

        object@cnv <- cnv

        amp <- apply(cnv, 2, function(x) paste0(names(which(x==1)), collapse="," ))
        del <- apply(cnv, 2, function(x) paste0(names(which(x==(-1))), collapse="," ))
        cnv_num <- apply(cnv, 2, function(x) length(which(x!=0)))

        seuratObj$cnv_num <- cnv_num

        cnv_casper <- paste0(amp, del, sep=",")
        cnv_casper_ann <- rep(1, length(cnv_casper))
        cnv_casper_ann[which(cnv_casper==",")] <- 0
        seuratObj$cnv_casper_ann <- cnv_casper_ann

    }
    if(loadXCVATR){
        load(paste0("./XCVATR/var_matrix_", project, ".rda"))
         if(!all((colnames(var_matrix)==colnames(seuratObj)))) stop("XCVATR sample names and SeuratObject sample names are different")

        var_matrix[var_matrix>1] <- 0
        object@snv <- var_matrix
        seuratObj$snv_num <- apply(var_matrix, 2, sum)

    }
    if(loadNumbat){
        if(file.exists(paste0("./cnv_casper/", project,"_cnv_numbat.rda"))) load(paste0("./cnv_casper/", project,"_cnv_numbat.rda"))
        if(!file.exists(paste0("./cnv_casper/", project,"_cnv_numbat.rda"))) stop("Numbat file do not exist")
        if(!all((colnames(cnv_numbat)==colnames(seuratObj)))) stop("NUMBAT sample names and SeuratObject sample names are different")

      #  cnv_num <- apply(cnv_numbat, 2, function(x) length(which(x!=0)))
        seuratObj$cnv_numbat <- as.vector(cnv_numbat)
        object@cnv <- rbind(object@cnv, cnv_numbat)
        seuratObj$cnv_num <- apply( object@cnv, 2, sum)


    }

    expr_features <- generate_HighExpressionFeaturesFromModel(seuratObj, tumor_markers=model_genes , k=3)
    tumor_m <- expr_features$tumor_markers
    rownames(tumor_m) <- paste0("HIGH_",rownames(tumor_m))
    object@tumor <- tumor_m
    tumor_m <- apply(tumor_m, 2, sum)
    tumor_m[tumor_m<2]<-0
    tumor_m[tumor_m==2]<-1
    #object@tumor <- tumor_m
    seuratObj$tumor_m <- as.vector(tumor_m)

    tumors <-   tumor_m+ simplified_models["tumor",] +simplified_models["Developmental_GSC",] + simplified_models["Injury_GSC",] +seuratObj$cnv_num + seuratObj$snv_num
    seuratObj$isTumor <- "non_tumor"
    seuratObj$isTumor[tumors>1] <- "tumor"
    seuratObj$num_feat <- as.vector(tumors)

    simplified_models["tumor",] <- simplified_models["tumor",] +simplified_models["Developmental_GSC",] +simplified_models["Injury_GSC",]

    seuratObj$isImmune <- "no"
    seuratObj$isImmune[simplified_models["immune", ]==1] <- "yes"
    models <- object@models
    prev <- ct[match(rownames(models),ct$celltype),"celltype_simplified"]
    models[which(prev %in% "immune"), seuratObj$isImmune=="no"] <- 0
    models[which(prev %in% "tumor"), seuratObj$isTumor=="non_tumor"] <- 0

    filter <- c(which(prev %in% "other"))

    models <- models[-filter, ]
    models <- rbind(models, as.numeric(as.factor(seuratObj$isTumor))-1)
    rownames(models)[length(rownames(models))] <- "final_tumor_annotation"
    models <- rbind(models, simplified_models["immune",])
    rownames(models)[length(rownames(models))] <- "final_immune_annotation"

    models["NPClike_suva", seuratObj$isTumor=="non_tumor"] <- 0
    models["OPClike_suva", seuratObj$isTumor=="non_tumor"] <- 0
    models["MESlike_suva", seuratObj$isTumor=="non_tumor"] <- 0
    models["AClike_suva", seuratObj$isTumor=="non_tumor"] <- 0


    simplified_models["NPClike_suva", seuratObj$isTumor=="non_tumor"] <- 0
    simplified_models["OPClike_suva", seuratObj$isTumor=="non_tumor"] <- 0
    simplified_models["MESlike_suva", seuratObj$isTumor=="non_tumor"] <- 0
    simplified_models["AClike_suva", seuratObj$isTumor=="non_tumor"] <- 0


    simplified_models["tumor", seuratObj$isTumor=="non_tumor"] <- 0
    simplified_models["immune", seuratObj$isImmune=="no"] <- 0

    seuratObj$cellclass_annotation <- apply(simplified_models, 2, function(x) paste0(rownames(simplified_models)[which(x>0)], collapse=","))
    seuratObj$celltype_annotation <- apply(models, 2, function(x) paste0(rownames(models)[which(x>0)], collapse=","))
    object@simplified_models <- simplified_models
    object@models <- models
    object@seurat_obj <- seuratObj
    object
}

runCellClassLevelSCRAM <- function(object) {
    seurat_obj <- object@seurat_obj
    models <- object@simplified_models
    cnv <- object@cnv
    snv <- object@snv
    tumor <- object@tumor
    freq_l2 <- list()
    ann_l <- list()


    if (is.null(seurat_obj$seurat_clusters)) {
            stop("seurat object clusters are missing")
    }

    clusters <- sort(as.numeric(as.character(unique(seurat_obj$seurat_clusters))))
    for (i in 1:length(clusters)) {
        cls <- clusters[i]
        index <- 1
        seuratObj_sub <- subset(seurat_obj, idents = (cls))
        print (paste0("running cluster ", cls,"..."))


        ann2 <- models[,  seurat_obj$seurat_clusters == (cls), drop = F]
        ann2 <- ann2[rownames(ann2) %in% names(which(apply(ann2, 1, sum) > 0)), , drop = F]
        freq_l2[[i]] <- data.frame()
        ann_l[[i]] <- data.frame()

        ann <- data.frame()
        snv_ann <- data.frame()
        tumor_ann <- data.frame()
        if (!is.null(cnv)) {
            ann <- cnv[, seurat_obj$seurat_clusters == (cls), drop = F]
        }

        if (!is.null(snv)) {
            snv_ann <- snv[, seurat_obj$seurat_clusters == (cls), drop = F]
        }

        if (!is.null(tumor)) {
            tumor_ann <- tumor[, seurat_obj$seurat_clusters == (cls), drop = F]

        }

        ann <- rbind(ann, snv_ann, tumor_ann)

        if (!is.null(ann2)) {
                ann <- rbind(ann2, ann)
        }

        if (!is.null(ann)) {
            if (nrow(ann) >= 0) {
                trans1 <- transactions(as.matrix(t(ann)))
                if (dim(trans1)[2] > 0) {
                    support <- object@min_support
                    #support <- object@min_support / dim(trans1)[1]
                   # if ((dim(trans1)[1] * 0.3) < object@min_support) support <- 0.3
                    freq <- eclat(trans1, parameter = list(support = support, maxlen = object@max_set_size), control = list(verbose = F))
                    if (length(freq) > 0) {
                        freq_l2[[i]] <- as(freq, "data.frame")[, -4]
                        freq_l2[[i]]$is.maximal <- is.maximal(freq)
                        freq_l2[[i]]$is.closed <- is.closed(freq)

                        freq_l2[[i]] <- freq_l2[[i]][with(freq_l2[[i]], order(-support, items)), ]
                        colnames(freq_l2[[i]]) <- c("cell_types", "cooccur_percentage", "num", "is.maximal", "is.closed")

                        ann_l[[i]] <- ann

                        names(freq_l2)[i] <- paste0("cluster", cls)
                        names(ann_l)[i] <- paste0("cluster", cls)
                    }
                }
            }
        }

    }

    object@setsAnnotatedCellClassLevel <- freq_l2
    object@cellClassLevelAnnotation <- ann_l
    object@cellClassAnnotationDetailed <- rep("", length(Idents(seurat_obj)))

    for (i in 1:length(object@cellClassLevelAnnotation)) {
        ann <- object@cellClassLevelAnnotation[[i]]
        if (dim(ann)[1] > 0) {
            cls <- as.numeric(gsub("cluster", "", names(object@cellClassLevelAnnotation)[i]))
            celltypes <- apply(ann, 2, function(x) paste0(unique(sort(names(x)[x == 1])), collapse = ","))
            object@cellClassAnnotationDetailed[Idents(seurat_obj) == cls] <- celltypes
        }
    }

    return(object)
}

runCellTypeLevelSCRAM <- function(object) {
    seurat_obj <- object@seurat_obj
    models <- object@models
    cnv <- object@cnv
    snv <- object@snv
    tumor <- object@tumor
    freq_l2 <- list()
    ann_l <- list()
    if (is.null(seurat_obj$seurat_clusters)) {
            stop("seurat object clusters are missing")
    }

    clusters <- sort(as.numeric(as.character(unique(seurat_obj$seurat_clusters))))
    for (i in 1:length(clusters)) {
        cls <- clusters[i]
        index <- 1
        seuratObj_sub <- subset(seurat_obj, idents = (cls))
        print (paste0("running cluster ", cls,"..."))


        ann2 <- models[,  seurat_obj$seurat_clusters == (cls), drop = F]
        ann2 <- ann2[rownames(ann2) %in% names(which(apply(ann2, 1, sum) > 0)), , drop = F]
        freq_l2[[i]] <- data.frame()
        ann_l[[i]] <- data.frame()

        ann <- data.frame()
        snv_ann <- data.frame()
        tumor_ann <- data.frame()
        if (!is.null(cnv)) {
            ann <- cnv[, seurat_obj$seurat_clusters == (cls), drop = F]
        }

        if (!is.null(snv)) {
            snv_ann <- snv[, seurat_obj$seurat_clusters == (cls), drop = F]
        }

        if (!is.null(tumor)) {
            tumor_ann <- tumor[, seurat_obj$seurat_clusters == (cls), drop = F]

        }

        ann <- rbind(ann, snv_ann, tumor_ann)

        if (!is.null(ann2)) {
                ann <- rbind(ann2, ann)
        }

        if (!is.null(ann)) {
            if (nrow(ann) >= 0) {
                trans1 <- transactions(as.matrix(t(ann)))
                if (dim(trans1)[2] > 0) {
                    support <- object@min_support
                    #support <- object@min_support / dim(trans1)[1]
                   # if ((dim(trans1)[1] * 0.3) < object@min_support) support <- 0.3
                    freq <- eclat(trans1, parameter = list(support = support, maxlen = object@max_set_size), control = list(verbose = F))
                    if (length(freq) > 0) {
                        freq_l2[[i]] <- as(freq, "data.frame")[, -4]
                        freq_l2[[i]]$is.maximal <- is.maximal(freq)
                        freq_l2[[i]]$is.closed <- is.closed(freq)

                        freq_l2[[i]] <- freq_l2[[i]][with(freq_l2[[i]], order(-support, items)), ]
                        colnames(freq_l2[[i]]) <- c("cell_types", "cooccur_percentage", "num", "is.maximal", "is.closed")

                        ann_l[[i]] <- ann

                        names(freq_l2)[i] <- paste0("cluster", cls)
                        names(ann_l)[i] <- paste0("cluster", cls)
                    }
                }
            }
        }

    }
    object@setsAnnotatedCellTypeLevel <- freq_l2
    object@cellTypeLevelAnnotation <- ann_l
    object@cellTypeLevelAnnotationDetailed <- rep("", length(Idents(seurat_obj)))

    for (i in 1:length(object@cellTypeLevelAnnotation)) {
        ann <- object@cellTypeLevelAnnotation[[i]]
        if (dim(ann)[1] > 0) {
            cls <- as.numeric(gsub("cluster", "", names(object@cellTypeLevelAnnotation)[i]))
            celltypes <- apply(ann, 2, function(x) paste0(unique(sort(names(x)[x == 1])), collapse = ","))
            object@cellTypeLevelAnnotationDetailed[Idents(seurat_obj) == cls] <- celltypes
        }
    }

    return(object)
}


annotateHighRibosomalCells <- function(object) {
    if (is.null(object@seurat_obj)) message("Skipping annotating high ribosomal cells. Seurat object is missing..")

    if (!is.null(object@seurat_obj)) {
        seurat_obj <- object@seurat_obj
        seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = "percent.ribo")
        object@cellTypeLevelAnnotationDetailed[object@cellTypeLevelAnnotationDetailed == "" & seurat_obj$percent.ribo > 40] <- "ribosomal"
     }
    return(object)
}


runSCRAM <- function(object) {

    print ("running cell type level SCRAM ")
    object <- runCellTypeLevelSCRAM(object)

    print ("running cell class level SCRAM ")
    object <- runCellClassLevelSCRAM(object)
    object <- annotateHighRibosomalCells(object)
    return(object)
}

