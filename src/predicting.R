#model
modelFile <- paste0(dir, "/src/model.RData")
load(modelFile)

#custom <- readRDS('./src/model.RData')

contigsId <- as.data.frame(tableFeatures$contig)
tableFeatures <- tableFeatures[,-12]

#dataset
closeAllConnections()

print("Predicting")
#removes classes for prediction
tableFeatures <- tableFeatures[,-12]

#prediction
predictClass <- predict(custom, tableFeatures)

predictClass <- as.data.frame(predictClass)

tableSequences <- cbind(tableSequences, predictClass)

tableSequencesYES <- tableSequences[which(tableSequences$predictClass == "yes") ,]
tableSequencesNO <- tableSequences[which(tableSequences$predictClass == "no") ,]

folder_and_name = paste0(id,"_genes.fasta")

tableSeq <- data.frame(tableSequencesYES)
name_seqs <- paste0(tableSeq$seq,", len=", tableSeq$tamanho,", ", tableSeq$contig)
tableSeq <- setNames(split(tableSeq[,1], seq(nrow(tableSeq))), rownames(tableSeq))
write.fasta(sequences = tableSeq, names = name_seqs, nbchar = 80, file.out = folder_and_name)

if(non == 1){
  print("Completed")
}
if(non == 2){
  folder_and_name = paste0(id,"_genes.fasta")
  tableSeq <- data.frame(tableSequencesYES)
  name_seqs <- paste0(tableSeq$seq,", len=", tableSeq$tamanho,", ", tableSeq$contig)
  tableSeq <- setNames(split(tableSeq[,1], seq(nrow(tableSeq))), rownames(tableSeq))
  write.fasta(sequences = tableSeq, names = name_seqs, nbchar = 80, file.out = folder_and_name)
  
  folder_and_name = paste0(id,"_intergenics.fasta")
  tableSeq <- data.frame(tableSequencesNO)
  name_seqs <- paste0(tableSeq$seq,", len=", tableSeq$tamanho,", ", tableSeq$contig)
  tableSeq <- setNames(split(tableSeq[,1], seq(nrow(tableSeq))), rownames(tableSeq))
  write.fasta(sequences = tableSeq, names = name_seqs, nbchar = 80, file.out = folder_and_name)
  print("Completed")
}



