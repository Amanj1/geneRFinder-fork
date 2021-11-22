sequencesFile <- NULL
sequencesFile <-  paste0(seqinr::read.fasta(fileDir, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE))

sequencesTable <- NULL
sequencesTable <- data.frame(matrix(unlist(sequencesFile), nrow=length(sequencesFile), byrow=T),stringsAsFactors=FALSE)

library(phylotools)
contigsNames <- as.data.frame(get.fasta.name(fileDir))
colnames(contigsNames) <- "contig"
               
closeAllConnections()

colnames(sequencesTable) <- c("sequencia")

print("Reading the fasta file")
tableSequencesPred <- data.frame(sequencia = NA, classe = NA, tamanho = NA, sentido = NA, org = NA, seq = NA)
tableSequencesPred <- cbind(sequencesTable, tableSequencesPred, contigsNames)

tableSequences <- findOrfs(tableSequencesPred, "cds", TRUE, startC)

rm(sequencesFile)
rm(sequencesTable)
rm(tableSequencesPred)
rm(contigsNames)
