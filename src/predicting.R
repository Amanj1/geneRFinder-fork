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

folder_and_name = paste0(id,".fasta")

tableSeq <- data.frame(tableSequences)
tableSeq <- setNames(split(tableSeq[,1], seq(nrow(tableSeq))), rownames(tableSeq))
write.fasta(sequences = tableSeq, names = names(tableSeq), nbchar = 80, file.out = folder_and_name)
print("Completed")



