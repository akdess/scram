
# Overview

**SCRAM** is available from https://github.com/akdess/scram 
**Updated version of CaSpER:**  is here:  [https://github.com/akdess/scram](https://github.com/akdess/casper_0.2.0/)
This vignette shows the basic steps for running SCRAM.

# Installation from github

Install latest version from GitHub (requires [devtools](https://github.com/hadley/devtools) package):

```r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("akdess/scram", dependencies = TRUE, build_vignettes = FALSE)
devtools::install_github("akdess/casper_0.2.0", dependencies = TRUE, build_vignettes = FALSE)
```




# Predicting CellTypes using Pretrained Neural Network Models

We first save the R object in h5ad for predicting cell-types on our trained deep learning models. We next predict the cell types on our data using our pretrained deep learning models.  


```
cd example;
<path_to_scripts_folder>/run_conversion.sh  glioma_seuratObj.rda  ./
python3 nn_classifier_pretrained.py 041524_glioma  ./adata.h5ad --normalize_test <path to nn_models folder> ./
```

We next load the seuratobject again in R with the predicted model outcomes. 

```r
setwd("example")
load("glioma_seuratObj.rda")
project <- "example"
scram_obj <- CreateSCRAMObject(seurat_obj=seuratObj,  organism="human", min_support=0.1, max_set_size=50) 
scram_obj <-createCellTypeMatrix (object=scram_obj, nn_path="./041524_glioma_multipleNeuralNetworks/", prob_thr=0.9, refs=c('suva_idh_a_o', 'hpa_brain_simple', 'allen_class_label_main', 'allen_neurons_only', 'TissueImmune', 'aldinger', 'codex', 'suva', 'bhaduri_withAge', 'dirks_primary_gbm_combined'), run="041524_glioma", pretrained=T)

```

## Annotating Tumor Cells
Because tumor cells exhibit a wide range of transcriptional states, we employ redundant and stringent approaches to annotate tumor cells using 3 modular components: (1) marker-expression modeling, (2) genotyping of CNVs on all cells (3) RNA-inferred mutational profiling of known glioma mutations (i.e. IDH1, EGFR). 

### Large Scale CNV calls in single cell resolution
CNVs are a hallmark feature of tumor cells that can be used to classify tumor vs. non-tumor cells alongside or in the absence of expression markers. However, detection of CNVs from scRNA-seq data is inherently noisy due to a multitude of factors, including drop-outs and unmatched control sets and requires a set of cells that are known to be tumor cells. To estimate a “clean” set of CNV calls that can provide reliable CNV-based tumor scores, we used a pure tumor pseudobulk sample.
Estimation of CNV profiles using patient-specific pure tumor pseudobulk samples. We first use our expression-based marker model from Module 1 to identify tumor cells. The collection of cells that are assigned as “tumor” using Module 1 is treated as a pure tumor cell cohort. 

CNV calling of patient-specific pure pseudobulk samples: We hypothesize that the pseudobulk sample contains representative sets of CNVs with high probability and therefore should be useful to identify a clean CNV call-set. The CNV calling on the pseudobulk samples is performed using our updated CNV calling algorithm, CaSpER+, for each patient. CaSpER+ CNV calls are used as the ground truth large-scale CNV calls for each patient. 
Genotyping of CNVs on all cells. After CNVs are identified from the pseudobulk sample, we genotype the set of CNVs on all cells and generate a binary matrix that represents the existence of CNVs on the cells, i.e., CNV_(i,j).
SeuratObj should have orig.ident with sample ids, and another metadata with id "tumorType" with "Normal" annotation for control cells 

Running single cell resolution CNV: After CNVs were identified from the pseudobulk sample, we genotyped the set of CNVs in all cells and generated a binary matrix that represents the existence of CNVs in the cells 

```r
### single cell level CNV calling  takes (~5 hours for 200K cells) long for large scRNA-Seq datasets. 
### the output is provided under cnv_casper folder
scram_obj <- runCASPER_Bulk_Scell(object=scram_obj, sampleCol="orig.ident", project)
```


### SNV calls with XCAVTR:

We performed RNA-inferred rare deleterious (COSMIC-reported and dbSNP, <0.1% frequency) mutational profiling via our recently developed XCVATR tool. 
Below code shows how to generate SNV input for SCRAM using XCAVTR output: 

```r
#### read XCVATR output
### the output is saved under XCVATR folder
var_matrix <- readXCVATROutput(XCVATR_folderPath="./XCVATR/", name_mapping)
```

Final Tumor annotation using CNV, SNV, expression modelling and neural network predictions

```r
load("cnv_casper/example_BULK_finalChrMat_thr_1.rda")
scram_obj <- runFinalTumorAnnotation(object=scram_obj, loadNumbat=F, loadCasper=T, finalChrMat_bulk=finalChrMat_bulk, loadXCVATR=T, sampleCol="orig.ident", project="example", model_genes=c("PDGFRA" ,"EGFR"  ,"SOX2" ))

```


We summarized co-occurring cell types using a frequent itemset rule mining approach. CNV and SNV calls were added to provide an integrated transcriptomic and genomic summary for each cell. An example SCRAM output for a single cell is given as “glioma stem cell, mature neuron, synaptic neuron, oligodendrocyte precursor cell, chr1p_deletion, chr19q_deletion + IDH1:2:208248389 mutation”. We used the tumour and host cell assignments of the previous steps to integrate co-occurring tumour and host cell features.

```r

scram_obj <- runSCRAM (object=scram_obj) 
writeSCRAMresults(object=scram_obj,project)

```

# Session info {.unnumbered}

Here is the output of `sessionInfo()` on the system on which this document was
compiled running pandoc `r rmarkdown::pandoc_version()`:

```{r sessionInfo, echo=FALSE}
sessionInfo()
```
