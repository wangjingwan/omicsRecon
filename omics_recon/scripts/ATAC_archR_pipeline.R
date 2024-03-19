# https://stuartlab.org/signac/articles/pbmc_vignette.html
# remotes::install_github("Bioconductor/GenomeInfoDb",lib = '/home/wangjingwan/R/x86_64-pc-linux-gnu-library/4.2')

library(argparse)
suppressMessages(library(ArchR))
set.seed(1)
addArchRThreads(threads = 10)


parser <- ArgumentParser(description='Check if the species is human or mouse')
parser$add_argument('-s', '--species', dest='species', required=TRUE, help='Set the species as human or mouse')
# parser$add_argument('-a', '--h5_path', dest='h5_path', required=FALSE, help='Path of the h5 file, storing Seurat object')
parser$add_argument('-f', '--frag_path', dest='frag_path', required=FALSE, help='Path of the fragment file needed for h5 file')
# parser$add_argument('-r', '--rds_path', dest='rds_path', required=FALSE, help='Path of the rds file, complete Seurat object')

parser$add_argument('-o', '--out_path', dest='out_path', required=FALSE, help='Set the path to the output file')

args <- parser$parse_args()

species = args$species
out_path = args$out_path

if(species == 'human'){
    # hg38
    addArchRGenome("hg38") 
}else{
    # mm10
    addArchRGenome("mm10")
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

# if(is.null(args$h5_path) & is.null(args$rds_path)) {
#     stop('Please provide either h5 file or rds file')
# }

# if(!is.null(args$rds_path)) {
#     rds_path = args$rds_path
#     pbmc <- readRDS(rds_path)
#     print('Object loaded')

#     # obj_str <- capture.output(str(pbmc))
#     # writeLines(obj_str, paste0(out_path,"strobj.txt"))
# }

inputFiles = args$frag_path
outDir = out_path
# inputFiles = './fragments_files/admae8w2_89i88tvv_atac_fragments_with_prefix.tsv.gz'
# outDir = './'

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = basename(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = outDir,
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

proj <- filterDoublets(ArchRProj = proj)
# get geneScore
geneScores = getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  threads = getArchRThreads()
)

count <- as.data.frame(as.matrix(assay(geneScores)))
row.names(count) = geneScores@elementMetadata$name
colnames(count) = sub(paste0(basename(inputFiles),'#'), "", colnames(count))
write.table(t(count),file=paste0(out_path,"archr_atac_pred_RNA.tsv"),row.names=T,col.names = T,sep = '\t',quote = F)
proj <- saveArchRProject(ArchRProj = proj)