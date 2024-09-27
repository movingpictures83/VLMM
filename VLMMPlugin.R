#context("alpine")

source("plugins/VLMM/core.R")
source("plugins/VLMM/fit_bias.R")
source("plugins/VLMM/helper.R")
source("plugins/VLMM/vlmm.R")
source("plugins/VLMM/plots.R")
source("plugins/VLMM/estimate_abundance.R")
source("plugins/VLMM/predict.R")
library(alpineData)
  library(GenomicAlignments)
  library(rtracklayer)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.NCBI.GRCh38)


dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  pfix <<- prefix()
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
}
run <- function() {}

output <- function(outputfile) {
  gap <- ERR188088()
  #dir <- tempdir()
  #print(dir)
  bam.file <- c("ERR188088" = file.path(paste(pfix, parameters["bamfile", 2], sep="/")))
  export(gap, con=bam.file)
  load(file=paste(pfix, parameters["preproc", 2], sep="/"))
  readlength <- as.integer(parameters["readlength", 2])
  minsize <- as.integer(parameters["minsize", 2])
  maxsize <- as.integer(parameters["maxsize", 2])
  gene.names <- names(ebt.fit)[6:8]
  names(gene.names) <- gene.names
  #fragtypes <- lapply(gene.names, function(gene.name) {
  #  buildFragtypes(ebt.fit[[gene.name]],
  #                 Hsapiens, readlength,
  #                 minsize, maxsize)
  #})
  #print(fragtypes)
  #saveRDS(fragtypes, outputfile)
  fragtypes <- readRDS(paste(pfix, parameters["fragments", 2], sep="/"))

  # model missing '+ gene' gives error
  #models <- list(
  #  "GC"=list(formula="count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk)",
  #            offset=c("fraglen","vlmm"))
  #)
    #fitBiasModels(
    #  genes=ebt.fit[gene.names],bam.file=bam.file,fragtypes=fragtypes,genome=Hsapiens,
    #  models=models,readlength=readlength,minsize=minsize,maxsize=maxsize
    #)


  # works to fit only fraglen and vlmm
  models <- list("readstart"=list(formula=NULL,offset=c("fraglen","vlmm")))
  #fitpar <- fitBiasModels(
  #  genes=ebt.fit[gene.names],bam.file=bam.file,fragtypes=fragtypes,genome=Hsapiens,
  #  models=models,readlength=readlength,minsize=minsize,maxsize=maxsize
  #)

  fitpar <- list(
    ERR188088 = fitBiasModels(
      genes=ebt.fit[gene.names],bam.file=bam.file,fragtypes=fragtypes,genome=Hsapiens,
      models=models,readlength=readlength,minsize=minsize,maxsize=maxsize
    )
  )
  # specify models using a character vector
  model.names <- c("fraglen","fraglen.vlmm","GC")
  txs <- txdf.theta$tx_id[txdf.theta$gene_id == parameters["gene", 2]]
			  #]"ENSG00000198918"]
  res <- estimateAbundance(transcripts=ebt.theta[txs],
                           bam.files=bam.file,
                           fitpar=fitpar,
                           genome=Hsapiens,
                           model.names=model.names)
  write.csv(res, outputfile)
}
