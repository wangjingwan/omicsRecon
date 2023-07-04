# https://stuartlab.org/signac/articles/pbmc_vignette.html
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
args = commandArgs(T)
species = args[4]
if(species == 'human'){
    suppressMessages(library(EnsDb.Hsapiens.v75))
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
}else{
    suppressMessages(library(EnsDb.Mmusculus.v79))
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
}
# load files
# frag_path = './examples/ATAC/10X_human_brain/human_brain_3k_atac_fragments.tsv.gz'
# h5_path = './examples/ATAC/10X_human_brain/human_brain_3k_filtered_feature_bc_matrix.h5'
h5_path = args[1]
frag_path = args[2]
out_path = args[3]

counts <- Read10X_h5(h5_path)
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)
seqlevelsStyle(annotation) <- "UCSC"
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = frag_path,
  annotation = annotation
)

# check original annotation
# head(Fragments(pbmc[["ATAC"]])[[1]])
# head(Annotation(pbmc[["ATAC"]]))
# genome(Annotation(pbmc[["ATAC"]]))
# head(pbmc@assays$ATAC@fragments[[1]])
# head(pbmc@assays$ATAC@ranges@seqinfo@genome)
# 

DefaultAssay(pbmc) <- "ATAC"
# time consuming
pbmc <- NucleosomeSignal(pbmc) 
pbmc <- TSSEnrichment(pbmc)
# filtering
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
# getting RNA activity
gene.activities <- GeneActivity(pbmc)
pbmc[['peakRNA']] <- CreateAssayObject(counts = gene.activities)
rna = as.matrix(pbmc@assays$peakRNA@counts)
write.table(rna,file="atac_RNA.csv",row.names=T,col.names = T,sep = ',',quote = F)