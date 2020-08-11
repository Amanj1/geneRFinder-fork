#############################################################################################
####################### Extract Features fom CDS and intergenics orfs #######################
#############################################################################################
#print("Extract features")
print("Extracting the features")
tableFeatures<- getFeatures(tableSequences)
tableFeatures<- tableFeatures[sample(nrow(tableFeatures)),] #To shuffle the data

