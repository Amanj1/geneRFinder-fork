sequencesFile <- NULL
sequencesFile <-  paste0(read.fasta(fileDir, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE))

sequencesTable <- NULL
sequencesTable <- data.frame(matrix(unlist(sequencesFile), nrow=length(sequencesFile), byrow=T),stringsAsFactors=FALSE)

closeAllConnections()

colnames(sequencesTable) <- c("sequencia")

print("Reading the fasta file")
tableSequencesPred <- data.frame(sequencia = NA, classe = NA, tamanho = NA, sentido = NA, org = NA, seq = NA)
tableSequencesPred <- cbind(sequencesTable, tableSequencesPred)

tableSequences <- findOrfs(tableSequencesPred, "cds", TRUE, startC)

rm(sequencesFile)
rm(sequencesTable)
rm(tableSequencesPred)
