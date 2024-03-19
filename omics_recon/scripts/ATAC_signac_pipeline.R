# https://stuartlab.org/signac/articles/pbmc_vignette.html
# remotes::install_github("Bioconductor/GenomeInfoDb",lib = '/home/wangjingwan/R/x86_64-pc-linux-gnu-library/4.2')

library(argparse)
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))

parser <- ArgumentParser(description='Check if the species is human or mouse')
parser$add_argument('-s', '--species', dest='species', required=TRUE, help='Set the species as human or mouse')
parser$add_argument('-a', '--h5_path', dest='h5_path', required=FALSE, help='Path of the h5 file, storing Seurat object')
parser$add_argument('-f', '--frag_path', dest='frag_path', required=FALSE, help='Path of the fragment file needed for h5 file')
parser$add_argument('-r', '--rds_path', dest='rds_path', required=FALSE, help='Path of the rds file, complete Seurat object')

parser$add_argument('-o', '--out_path', dest='out_path', required=FALSE, help='Set the path to the output file')

args <- parser$parse_args()

species = args$species
out_path = args$out_path

if(species == 'human'){
    # hg38
    suppressMessages(library(EnsDb.Hsapiens.v75))
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
}else{
    suppressMessages(library(EnsDb.Mmusculus.v79))
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
}


if (is.null(out_path)){
    out_path = './'
}else{
    out_path = args$out_path
}

# if not exsist, creat
if (!file.exists(out_path)){
    dir.create(out_path)
}

if(is.null(args$h5_path) & is.null(args$rds_path)) {
    stop('Please provide either h5 file or rds file')
}

if(!is.null(args$rds_path)) {
    rds_path = args$rds_path
    pbmc <- readRDS(rds_path)
    print('Object loaded')

    # obj_str <- capture.output(str(pbmc))
    # writeLines(obj_str, paste0(out_path,"strobj.txt"))
}

if(!is.null(args$h5_path)) {
    if(is.null(args$frag_path)) {
        stop('Please provide both h5 file and fragment file')
    }
    h5_path = args$h5_path
    frag_path = args$frag_path
    # frag_path = './human_brain_3k_atac_fragments.tsv.gz'
    # h5_path = './human_brain_3k_filtered_feature_bc_matrix.h5'
    counts <- Read10X_h5(h5_path)
    pbmc <- CreateSeuratObject(
      counts = counts$`Gene Expression`,
      assay = "RNA"
    )
    # having trouble swith to ucsc version, using ncbi to procede
    # seqlevelsStyle(annotation) <- "UCSC"
    pbmc[["ATAC"]] <- CreateChromatinAssay(
      counts = counts$Peaks,
      sep = c(":", "-"),
      fragments = frag_path,
      annotation = annotation
    )
    print('Object loaded')
}

# check original annotation
# head(Fragments(pbmc[["ATAC"]])[[1]])
# head(Annotation(pbmc[["ATAC"]]))
# print(genome(Annotation(pbmc[["ATAC"]])))
# print(head(pbmc@assays$ATAC@fragments[[1]]))
# head(pbmc@assays$ATAC@ranges@seqinfo@genome)
# 

DefaultAssay(pbmc) <- "ATAC"
# time consuming
pbmc <- NucleosomeSignal(pbmc) 
pbmc <- TSSEnrichment(pbmc)
# TODO del this checking line
obj_str <- capture.output(str(pbmc))
writeLines(obj_str, paste0(out_path,"strobj.txt"))

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
rna = t(as.matrix(pbmc@assays$peakRNA@counts))
write.table(rna,file=paste0(out_path,"atac_pred_RNA.tsv"),row.names=T,col.names = T,sep = '\t',quote = F)
write.table(pbmc@meta.data, file = paste0(out_path,"signac_atac_pred_meta.tsv"), sep = "\t", quote = F, row.names = T, col.names = T)
saveRDS(pbmc, file = paste0(out_path,"signac_processed_atac.rds"))