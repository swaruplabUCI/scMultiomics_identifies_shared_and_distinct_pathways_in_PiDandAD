
# script to run cicero in parallel for all clusters:
library(optparse)

option_list = list(
  make_option(
    c('-f', '--file'), type='character', default=NULL,
    help='Cell Dataset .rds file', metavar='character'
  ),
  make_option(
    c('-o', '--outdir'), type='character', default='./',
    help='Directory to place output files', metavar='character'
  ),
  make_option(
    c('-c', '--clusters'), type='character', default=NULL,
    help='name of the cluster subset of the CDS to process cicero'
  ),
  make_option(
    c('-n', '--num'), type='numeric', default=NULL,
    help='SLURM task array number goes here, selects which cluster to process on this task'
  ),
  make_option(
    c('-g', '--genome'), type='character', default=NULL,
    help='Genome sizes file, /PATH/FILE'
  )
)

# parse arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

################################################################################

print(opt)

print('Loading libraries...')
library(Signac);
library(Seurat);
library(SeuratWrappers);
library(monocle3);
library(cicero);
print('Finished loading libraries.')

# load CDS:
print(paste('Loading CDS from file', opt$file))
NucSeq.atac_cds <- readRDS(file=opt$file)
print('Finished loading CDS')

# get a list of clusters
clusters <- unique(as.character(pData(NucSeq.atac_cds)[[opt$clusters]]))
clusters <- clusters[order(clusters)]

print(clusters)

# load genome sizes:
human.hg38.genome <- read.csv(file=opt$genome, sep='\t', header=TRUE)

# select cluster based on SLURM task array ID
cluster <- clusters[opt$num]
print(paste('Processing on group', cluster))

# subset cell dataset object by current cluster
print('Subsetting CDS...')
cur_cds <- NucSeq.atac_cds[,rownames(pData(NucSeq.atac_cds)[pData(NucSeq.atac_cds)[[opt$clusters]] == cluster,])]
print(paste('CDS dimensions:', dim(cur_cds)))

# subset by Diagnosis
print('Subsetting by Diagnosis...')
cur_cds_AD <- cur_cds[, rownames(pData(NucSeq.atac_cds)[pData(NucSeq.atac_cds)[[opt$clusters]] == cluster & pData(NucSeq.atac_cds)$Diagnosis == 'AD',])]
cur_cds_control <- cur_cds[, rownames(pData(NucSeq.atac_cds)[pData(NucSeq.atac_cds)[[opt$clusters]] == cluster & pData(NucSeq.atac_cds)$Diagnosis == 'Control',])]
print(paste('AD CDS dimensions:', dim(cur_cds_AD)))
print(paste('Control CDS dimensions:', dim(cur_cds_control)))


# make cicero cds
print('Making cicero CDS...')
# cur_cicero_cds <- make_cicero_cds(cur_cds, reduced_coordinates = reducedDims(cur_cds)$UMAP)
cur_cicero_cds_AD <- make_cicero_cds(cur_cds_AD, reduced_coordinates = reducedDims(cur_cds_AD)$UMAP)
cur_cicero_cds_control <- make_cicero_cds(cur_cds_control, reduced_coordinates = reducedDims(cur_cds_control)$UMAP)
print('Done making cicero CDS.')

# save cicero CDS:
# saveRDS(cur_cicero_cds, file=paste0(opt$outdir,'/',cluster,'_cicero_cds.rds'))
saveRDS(cur_cicero_cds_AD, file=paste0(opt$outdir,'/',cluster,'_cicero_cds_AD.rds'))
saveRDS(cur_cicero_cds_control, file=paste0(opt$outdir,'/',cluster,'_cicero_cds_control.rds'))
print('Done!')

# compute co-accessibility:
print('Computing co-accessibility')
# conns <- run_cicero(cur_cicero_cds, human.hg38.genome)
conns_AD <- run_cicero(cur_cicero_cds_AD, human.hg38.genome, window = 3e+05, sample_num = 200)
conns_control <- run_cicero(cur_cicero_cds_control, human.hg38.genome, window = 3e+05, sample_num = 200)
print('Done computing co-accessibility')

# save all the connections
print('Saving results')
# save(conns, connections, conns_AD, connections_AD, conns_control, connections_control, file=paste0(opt$outdir,'/',cluster,'_cicero_connections_qc.rda'))
# save(conns, conns_AD, conns_control, file=paste0(opt$outdir,'/',cluster,'_cicero_connections.rda'))
save(conns_AD, conns_control, file=paste0(opt$outdir,'/',cluster,'_cicero_connections.rda'))
