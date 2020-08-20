###############################
##### FILE FOR PREDICTION #####
###############################

suppressMessages(library(optparse))
suppressMessages(library(seqinr))
suppressMessages(library(kmer))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
suppressMessages(library(parallel))
suppressMessages(library(randomForest))
suppressMessages(library(caret))

seed <- 7
set.seed(seed)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input File", metavar="character"),
  make_option(c("-t", "--thread"), type="integer", default=2, 
              help="Number of Threads", metavar="integer"),
  make_option(c("-o", "--out"), type="character", default="genes", 
              help="Output File Name [default= %default]", metavar="character"),
  make_option(c("-s", "--start"), type="integer", default=1, 
              help="Type of start codon", metavar="integer"),
  make_option(c("-int", "--intergenic"), type="integer", default=1, 
              help="Intergenics sequences", metavar="integer")
  
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Please, input a fasta file.", call.=FALSE)
}

# Number of threads
#cl <- 7
cl <- as.integer(opt$thread)
registerDoParallel(cl)
cltr <- makeCluster(cl)

# Identification 
id = opt$out

# Fasta file for predict
fileDir <- opt$input

# Type of start codon
startC <- as.numeric(opt$start)

#Output with intergenic
non <- as.numeric(opt$intergenic)

this_file = gsub("--file=", "", commandArgs()[grepl("--file", commandArgs())])
if (length(this_file) > 0){
  wd <- paste(head(strsplit(this_file, '[/|\\]')[[1]], -1), collapse = .Platform$file.sep)
}else{
  wd <- dirname(rstudioapi::getSourceEditorContext()$path)
}

dir <- wd


loadF <- paste0(dir, "/src/functions.R")
source(loadF)
readfileF <- paste0(dir, "/src/readFastaFile.R")
source(readfileF)
getFeatF <- paste0(dir, "/src/getTableFeatures.R")
source(getFeatF)
predF <- paste0(dir, "/src/predicting.R")
source(predF)
