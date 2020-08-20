#model
modelFile <- paste0(dir, "/src/model.RData")
load(modelFile)

#custom <- readRDS('./src/model.RData')

#dataset
closeAllConnections()

print("Predicting")
#removes classes for prediction
tableFeatures <- tableFeatures[,-12]

#prediction
predictClass <- predict(custom, tableFeatures)

predictClass <- as.data.frame(predictClass)

tableSequences <- cbind(tableSequences, predictClass)

tableSequences <- tableSequences[which(tableSequences$predictClass == "yes") ,]
tableSequencesNO <- tableSequences[which(tableSequences$predictClass == "no") ,]

folder_and_name = paste0(id,"_genes.fasta")

tableSeq <- data.frame(tableSequences)
tableSeq <- setNames(split(tableSeq[,1], seq(nrow(tableSeq))), rownames(tableSeq))
write.fasta(sequences = tableSeq, names = names(tableSeq), nbchar = 80, file.out = folder_and_name)

if(non == 1){
  print("Completed")
}
if(non == 2){
  folder_and_name = paste0(id,"_intergenics.fasta")
  tableSeq <- data.frame(tableSequencesNO)
  tableSeq <- setNames(split(tableSeq[,1], seq(nrow(tableSeq))), rownames(tableSeq))
  write.fasta(sequences = tableSeq, names = names(tableSeq), nbchar = 80, file.out = folder_and_name)
  print("Completed")
}



