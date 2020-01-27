#############################################################################################
############################### Functions used in the pipeline ##############################
#############################################################################################


#Obtains feature table
obterFeatureTableCoord <- function(nameFIle){
  t<-read.csv(nameFIle,sep="\t",header=F)
  x<-which(t$V3 == "CDS")
  featureTableCoord<-t[x,]
  featureTableCoord<-featureTableCoord[,1:3]
  
  #delete <e> character
  featureTableCoord$V1<-as.numeric(regmatches(featureTableCoord[,1],gregexpr('[0-9]+',featureTableCoord[,1])))
  featureTableCoord$V2<-as.numeric(regmatches(featureTableCoord[,2],gregexpr('[0-9]+',featureTableCoord[,2])))
  
  #Direction of sequence
  featureTableCoord<-cbind(featureTableCoord,"sentido"=NA)
  featureTableCoord[which(as.numeric(featureTableCoord[,1]) < as.numeric(featureTableCoord[,2])),4]<-"fw"
  featureTableCoord[which(as.numeric(featureTableCoord[,1]) > as.numeric(featureTableCoord[,2])),4]<-"bw"
  
  return(featureTableCoord)
}

#Obtains coding sequences
preencherComSeqsCodificadoras <- function(fasta){
  seq<-NULL
  gene<-NULL
  inicio<-1
  tableSeq<-matrix(ncol = 7)
  colnames(tableSeq)<-c("sequence","compl", "seq", "tamanho", "sentido", "org", "class")
  tableSeq<-data.frame(tableSeq)
  linhaCabecalho<- 1
  
  for(i in 1: length(fasta)){
    if(isTRUE(as.logical(substr(fasta[i],1,1) != ">"))){
      seq<-fasta[i]
      gene<-paste0(gene,seq,sep="")
      if(i == length(fasta)){
        if(as.character(regmatches(fasta[linhaCabecalho], gregexpr("complement",fasta[linhaCabecalho]))) == "complement"){
          tableSeq[inicio,"sentido"]<- "bw"
          
        }else tableSeq[inicio,"sentido"]<- "fw" 
        
        tableSeq[inicio,"sequence"]<-gene
        tableSeq[inicio,"class"]<-"yes"
        tableSeq[inicio,"tamanho"]<-length(s2c(tableSeq[inicio,1]))
        tableSeq[inicio,"seq"] <-inicio
        tableSeq[inicio,"compl"] <- 1
        inicio<-inicio+1
        gene<-NULL
        seq<-NULL
      }
    }else{
      if(i != 1){
        if(as.character(regmatches(fasta[linhaCabecalho], gregexpr("complement",fasta[linhaCabecalho]))) == "complement"){
          tableSeq[inicio,"sentido"]<- "bw"
        }else tableSeq[inicio,"sentido"]<- "fw" 
        
        tableSeq[inicio,"sequence"]<-gene
        tableSeq[inicio,"class"]<-"yes"
        tableSeq[inicio,"tamanho"]<-length(s2c(tableSeq[inicio,1]))
        tableSeq[inicio,"seq"]<-inicio
        tableSeq[inicio,"compl"] <- 1
        tableSeq[inicio,"sentido"]<-
          inicio<-inicio+1
        gene<-NULL
        seq<-NULL
        linhaCabecalho<- i
      }
    }   
  }
  return(tableSeq)
}


#Obtains genome
obterCompleteGenome <- function(nameFIle){
  openFile<-file(nameFIle, open="r")
  genome<-readLines(openFile)
  complete<-completeGenome(genome)
  
  return(complete)
}


completeGenome<-function(fastaCompleteGenome){
  
  aux<-fastaCompleteGenome[-1]
  genome<-paste0(aux,collapse="")
  return(genome)
  
}

#Extracts intergenic regions from forward sequencing
extractSeqsNegativasFita1 <- function(featureTableCoord,genome){
  tableSeqNegFita1<-matrix(ncol=6)
  colnames(tableSeqNegFita1)<-c("sequence","class","tamanho","sentido","org","seq")
  tableSeqNegFita1<-data.frame(tableSeqNegFita1)
  
  coordGenesFita1<-montarTableCoodGenesFita1(featureTableCoord)
  coordNaoGenesFita1<-montarTableCoodNaoGenesFita1(coordGenesFita1,genome)
  for(i in 1: nrow(coordNaoGenesFita1)){
    tableSeqNegFita1[i,1]<-substr(genome,coordNaoGenesFita1[i,1],coordNaoGenesFita1[i,2])
    tableSeqNegFita1[i,2]<-"no"
    tableSeqNegFita1[i,3]<-length(s2c(tableSeqNegFita1[i,1]))
    tableSeqNegFita1[i,6]<-i
  }
  
  return(tableSeqNegFita1)
}

#Assembly table of non-genes from forward sequencing
montarTableCoodNaoGenesFita1 <- function(coordGenesFita1,genome){
  coordNaoGenesFita1<-matrix(ncol=2)
  colnames(coordNaoGenesFita1)<-c("ini","fim")  
  coordNaoGenesFita1<-data.frame(coordNaoGenesFita1)
  
  vectorCompleteGenome<-s2c(genome)
  cont<-1
  for(i in 1: nrow(coordGenesFita1)){
    if(i == nrow(coordGenesFita1)){
      if(coordGenesFita1[i,2] != length(genome)){
        coordNaoGenesFita1[cont,1]<-coordGenesFita1[i,2]+1
        coordNaoGenesFita1[cont,2]<-length(vectorCompleteGenome)
      }
    }else{
      n<-i+1
      if(isTRUE(coordGenesFita1[i,2] < coordGenesFita1[n,1] & coordGenesFita1[i,2] - coordGenesFita1[n,1] != -1)){
        coordNaoGenesFita1[cont,1]<-coordGenesFita1[i,2]+1
        coordNaoGenesFita1[cont,2]<-coordGenesFita1[n,1]-1
        cont<-cont+1
      }
    }
  }
  return(coordNaoGenesFita1)
}

#Extracts intergenic regions from reverse sequencing
extractSeqsNegativasFitaComp <- function(featureTableCoord,genome){
  
  tableSeqNegFitaComp<-matrix(ncol=6)
  colnames(tableSeqNegFitaComp)<-c("sequence","class","tamanho","sentido","org","seq")
  tableSeqNegFitaComp<-data.frame(tableSeqNegFitaComp)
  
  coordGenesFitaComp<-montarTableCoodGenesFitaComp(featureTableCoord)
  coordNaoGenesFitaComp<-montarTableCoodNaoGenesFitaComp(coordGenesFitaComp,genome)
  
  seq<-NULL
  seqVetor<-NULL
  seqVetorInvertido<-NULL
  vectorNegFitaComp<-NULL
  for(i in 1: nrow(coordNaoGenesFitaComp)){
    #coord[i,]<-sort(coord[i,])
    seq<-substr(genome,coordNaoGenesFitaComp[i,1],coordNaoGenesFitaComp[i,2])
    if(seq != ""){
      seqVetor<-s2c(seq)
      seqVetorInvertido<-rev(seqVetor)
      vectorNegFitaComp<-inverteNucleotideo(seqVetorInvertido)
      seqNegFitaComp<-paste(vectorNegFitaComp,collapse="")
      tableSeqNegFitaComp[i,1]<-seqNegFitaComp
      tableSeqNegFitaComp[i,2]<-"no"
      tableSeqNegFitaComp[i,3]<-length(s2c(tableSeqNegFitaComp[i,1]))
      tableSeqNegFitaComp[i,6]<-i
    }
  }
  
  return(tableSeqNegFitaComp)
}


#genarates reverse sequencing
inverteNucleotideo <- function(sequence){
  
  sequence<-comp(sequence, forceToLower=FALSE)
  sequence<-paste(sequence, sep="", collapse="")
  return(sequence)
  
}

#Assembly table of non-genes from reverse sequencing
montarTableCoodNaoGenesFitaComp <- function(coordGenesFitaComp,genome){
  coordNaoGenesFitaComp<-matrix(ncol=2)
  colnames(coordNaoGenesFitaComp)<-c("ini","fim")  
  coordNaoGenesFitaComp<-data.frame(coordNaoGenesFitaComp)
  
  
  if(suppressWarnings(min(coordGenesFitaComp[which(coordGenesFitaComp[,"sentido"] == "bw"),2])) > 1){
    coordNaoGenesFitaComp[1,1]<-1
    coordNaoGenesFitaComp[1,2]<-coordGenesFitaComp[1,2]-1
    cont<-2
  }else cont<-1
  
  for(i in 1: nrow(coordGenesFitaComp)){
    if(cont == nrow(coordGenesFitaComp)){
      if(isTRUE(as.logical(coordGenesFitaComp[i,2] != length(genome)))){
        coordNaoGenesFitaComp[cont,1]<-coordGenesFitaComp[i,2]+1
        coordNaoGenesFitaComp[cont,2]<-length(vectorCompleteGenome)
      }
    }else{
      n<-i+1
      if((coordGenesFitaComp[i,1] > coordGenesFitaComp[n,2]) == FALSE & nrow(coordGenesFitaComp) != i){
        if(coordGenesFitaComp[i,1] - coordGenesFitaComp[n,2] != -1 & coordGenesFitaComp[i,1] - coordGenesFitaComp[n,2] != 0){
          coordNaoGenesFitaComp[cont,1]<-coordGenesFitaComp[i,1]+1
          coordNaoGenesFitaComp[cont,2]<-coordGenesFitaComp[n,2]-1
          cont<-cont+1
        }
      }else{
      }
    }    
  }
  return(coordNaoGenesFitaComp)
}

#Assembly coordinates table of genes from forward sequencing
montarTableCoodGenesFita1 <- function(featureTableCoord){
  coordGenesFita1<-matrix(ncol=2)
  colnames(coordGenesFita1)<-c("ini","fim")  
  coordGenesFita1<-data.frame(coordGenesFita1)
  
  featureTableCoord$V1<-as.integer(as.character(featureTableCoord$V1))
  featureTableCoord$V2<-as.integer(as.character(featureTableCoord$V2))
  
  cont<-1
  for(i in 1: nrow(featureTableCoord)){
    if(featureTableCoord[i,1] < featureTableCoord[i,2]){
      coordGenesFita1[cont,1]<-featureTableCoord[i,1]
      coordGenesFita1[cont,2]<-featureTableCoord[i,2]
      cont<-cont+1
    }
  }
  return(coordGenesFita1)
}

#creates coordinate table of genes from reverse sequencing
montarTableCoodGenesFitaComp <- function(featureTableCoord){
  coordGenesFitaComp<-matrix(ncol=3)
  colnames(coordGenesFitaComp)<-c("ini","fim","sentido")                                                                        
  coordGenesFitaComp<-data.frame(coordGenesFitaComp)
  
  featureTableCoord$V1<-as.integer(as.character(featureTableCoord$V1))
  featureTableCoord$V2<-as.integer(as.character(featureTableCoord$V2))
  
  cont<-0
  for(i in 1: nrow(featureTableCoord)){
    if(featureTableCoord[i,1] > featureTableCoord[i,2]){
      cont<-cont+1
      coordGenesFitaComp[cont,1]<-featureTableCoord[i,1]
      coordGenesFitaComp[cont,2]<-featureTableCoord[i,2]
      coordGenesFitaComp[cont,3]<-featureTableCoord[i,3]
    }  
  }
  return(coordGenesFitaComp)
}

#search ORFs
findOrfs<-function(sequences, origin, needSubOrfs, startC){
  minlength = 60
  if(origin == "cds"){
    orfs<- orfsFromSequences(sequences, needSubOrfs, startC, minlength)
    orfs$class<- "yes"
  }
  
  if(origin == "intergenic"){
    orfs<- orfsFromSequences(sequences, needSubOrfs, startC, minlength)
    orfs$class<- "no"
  }
  
  return(orfs)
}

#Extracts ORFs from forward sequencing
orfsFromSequences<- function(sequences, needSubOrfs, startC, minlength){
  
  sequences <- sequences[!is.na(sequences[,1]),]
  tableSeqOrfs <- data.frame(sentido = NA, org = NA, class = NA, stringsAsFactors = FALSE)
  
  #tableSeqOrfs<- tableSeqOrfs[-1,]
  sequences[,1]<- toupper(sequences[,1])
  
  if(needSubOrfs == FALSE){
    tableOrfs <- foreach(j=1:nrow(sequences), .combine= 'rbind', .packages='seqinr') %dopar% {
      seq<- sequences[j,1]
      
      if(startC == 2){
        matches <- gregexpr("(?=((ATG|GTG|TTG)(?:[ATGC|GTGC|TTGC]{3})*?)(?=TAA|TAG|TGA))",seq,perl=TRUE)
      }else{
        matches <- gregexpr("(?=(ATG(?:[ATGC]{3})*?)(?=TAA|TAG|TGA))",seq,perl=TRUE)
      }
      
      startPositions  <- as.vector(matches[[1]])
      lengths         <- c(attr(matches[[1]], "match.length"))
      stopPosition<- startPositions + lengths -1
      
      orfSequence <- as.matrix(substring(seq, first = startPositions, last  = stopPosition))
    }
    colnames(tableOrfs)<- "sequence"
    tableOrfs <- as.data.frame(tableOrfs)
    sequencesNumber <- as.data.frame(rownames(tableOrfs))
    colnames(sequencesNumber)<- "seq"
    tableOrfs$compl <- 1
    
    
  }else{
    tableOrfs <- foreach(j=1:nrow(sequences), .combine= 'rbind', .packages='seqinr') %dopar% {
      seq<- sequences[j,1]
      
      if(startC == 2){
        matches <- gregexpr("(?=((ATG|GTG|TTG)(?:[ATGC|GTGC|TTGC]{3})*?)(?=TAA|TAG|TGA))",seq,perl=TRUE)
      }else{
        matches <- gregexpr("(?=(ATG(?:[ATGC]{3})*?)(?=TAA|TAG|TGA))",seq,perl=TRUE)
      }
      
      startPositions  <- as.vector(matches[[1]])
      lengths         <- c(attr(matches[[1]], "capture.length"))
      stopPosition<- (startPositions + lengths -1) + 3
      
      
      orfSequence <- as.matrix(substring(seq, first = startPositions, last  = stopPosition))
    }
    colnames(tableOrfs)<- "sequence"
    tableOrfs <- as.data.frame(tableOrfs, stringsAsFactors = FALSE)
    sequencesNumber <- as.data.frame(rownames(tableOrfs), stringsAsFactors = FALSE)
    colnames(sequencesNumber)<- "seq"
    tableOrfs$compl <- 2
  }
  
  
  tableOrfs <- cbind(tableOrfs, sequencesNumber)
  tableOrfs <- as.data.frame(tableOrfs[which(tableOrfs$sequence != ""),], stringsAsFactors = FALSE)
  tableOrfs <- as.data.frame(tableOrfs[!duplicated(tableOrfs$sequence), ], stringsAsFactors = FALSE)
  
  # j<-1
  tamanhoSeq <- NULL
  tamanhoSeq <- foreach(j=1:nrow(tableOrfs), .combine= 'rbind', .packages='seqinr') %dopar% {
    tamanhoSeq[j] <- nchar(as.character(tableOrfs[j,1]))
  } 
  
  tamanho <- as.data.frame(tamanhoSeq, stringsAsFactors = FALSE)
  colnames(tamanho)<- "tamanho"
  tableOrfs <- cbind(tableOrfs, tamanho, tableSeqOrfs)
  tableOrfs[which.max(tableOrfs[,"tamanho"]),"compl"]<- 1
  tableOrfs<- tableOrfs[which(tableOrfs[,"tamanho"] >= minlength),]
  
  
  return(tableOrfs)
}


variancia <- function(x) {
  media <- mean(x)
  n <- length(x)
  var <- sum((x - media)^2)/n
  return(var)
}

hexamerF <- function(seq){
  matTest <- data.frame(NULL)
  aux <- kmer::kcount(seq, k = 6)
  resHexamer <- variancia(aux)
  return(resHexamer)
}

dimerF <- function(seq){
  matTest <- data.frame(NULL)
  aux <- kmer::kcount(seq, k = 2)
  resDimer <- variancia(aux)
  return(resDimer)
}

trimerF <- function(seq){
  matTest <- data.frame(NULL)
  aux <- kmer::kcount(seq, k = 3)
  resTrimer <- variancia(aux)
  return(resTrimer)
}

tetramerF <- function(seq){
  matTest <- data.frame(NULL)
  aux <- kmer::kcount(seq, k = 4)
  resTetramer <- variancia(aux)
  return(resTetramer)
}

pentamerF <- function(seq){
  matTest <- data.frame(NULL)
  aux <- kmer::kcount(seq, k = 5)
  resPentamer <- variancia(aux)
  return(resPentamer)
}

################# functions for features - begin ##########################
extractGC<-function(tableSequences){
  count<-1
  aux<- NULL
  resultado <- foreach(i=1:nrow(tableSequences), .combine= 'rbind', .packages='seqinr') %dopar% {
    seq <- s2c(as.character(tableSequences[i,1]))
    aux[i]<- GC(seq)  }
  return(resultado)
}

extractGC1<-function(tableSequences){
  count<-1
  aux<- NULL
  resultado <- foreach(i=1:nrow(tableSequences), .combine= 'rbind', .packages='seqinr') %dopar% {
    seq <- s2c(as.character(tableSequences[i,1]))
    aux[i]<- GC1(seq)  }
  return(resultado)
}

extractGC2<-function(tableSequences){
  count<-1
  aux<- NULL
  resultado <- foreach(i=1:nrow(tableSequences), .combine= 'rbind', .packages='seqinr') %dopar% {
    seq <- s2c(as.character(tableSequences[i,1]))
    aux[i]<- GC2(seq)  }
  return(resultado)
}

extractGC3<-function(tableSequences){
  count<-1
  aux<- NULL
  resultado <- foreach(i=1:nrow(tableSequences), .combine= 'rbind', .packages='seqinr') %dopar% {
    seq <- s2c(as.character(tableSequences[i,1]))
    aux[i]<- GC3(seq)  }
  return(resultado)
}

extractDimer<-function(tableSequences){
  count<-1
  aux<- NULL
  resultado <- foreach(i=1:nrow(tableSequences), .combine= 'rbind', .packages='kmer') %dopar% {
    seq <- s2c(as.character(tableSequences[i,1]))
    aux[i] <- dimerF(seq)   }
  return(resultado)
}

extractTrimer<-function(tableSequences){
  count<-1
  aux<- NULL
  resultado <- foreach(i=1:nrow(tableSequences), .combine= 'rbind', .packages='kmer') %dopar% {
    seq <- s2c(as.character(tableSequences[i,1]))
    aux[i] <- trimerF(seq)   }
  return(resultado)
}

extractTetramer<-function(tableSequences){
  count<-1
  aux<- NULL
  resultado <- foreach(i=1:nrow(tableSequences), .combine= 'rbind', .packages='kmer') %dopar% {
    seq <- s2c(as.character(tableSequences[i,1]))
    aux[i] <- tetramerF(seq)   }
  return(resultado)
}

extractPentamer<-function(tableSequences){
  count<-1
  aux<- NULL
  resultado <- foreach(i=1:nrow(tableSequences), .combine= 'rbind', .packages='kmer') %dopar% {
    seq <- s2c(as.character(tableSequences[i,1]))
    aux[i] <- pentamerF(seq)   }
  return(resultado)
}

extractHexamer<-function(tableSequences){
  count<-1
  aux<- NULL
  resultado <- foreach(i=1:nrow(tableSequences), .combine= 'rbind', .packages='kmer') %dopar% {
    seq <- s2c(as.character(tableSequences[i,1]))
    aux[i] <- hexamerF(seq)   }
  return(resultado)
}

extractC_weight<-function(tableSequences){
  count<-1
  aux<- NULL
  resultado <- foreach(i=1:nrow(tableSequences), .combine= 'rbind', .packages='seqinr') %dopar% {
    seq <- s2c(as.character(tableSequences[i,1]))
    codonusagew <- ucoweight(seq, numcode = 1)
    codonusagew <- as.table(unlist(codonusagew))
    aux[i] <- var(codonusagew) 
  }
  return(resultado)
}

extractLenght <-function(tableSequences){
  count<-1
  aux<- NULL
  resultado <- foreach(i=1:nrow(tableSequences), .combine= 'rbind', .packages='seqinr') %dopar% {
    seq <- s2c(as.character(tableSequences[i,1]))
    aux[i] <- length(seq) 
  }
  return(resultado)
}


################# functions for features - end ##########################
getFeatures <- function(tableSequences){
  
  GC <- extractGC(tableSequences)
  GC1 <- extractGC1(tableSequences)
  GC2 <- extractGC2(tableSequences)
  GC3 <- extractGC3(tableSequences)
  length <- extractLenght(tableSequences)
  dimer <- extractDimer(tableSequences)
  trimer <- extractTrimer(tableSequences)
  tetramer <- extractTetramer(tableSequences)
  pentamer <- extractPentamer(tableSequences)
  hexamer <- extractHexamer(tableSequences)
  c_weight <- extractC_weight(tableSequences)
  class <- tableSequences$class
  
  tableFeatures<- cbind(GC, GC1, GC2, GC3, length, dimer, trimer, tetramer, pentamer, hexamer, c_weight, class)
  tableFeatures <- as.data.frame(tableFeatures)
  
  colnames(tableFeatures) <- c("GC", "GC1", "GC2", "GC3", "length", "dimer", "trimer", "tetramer", "pentamer", 
                               "hexamer", "c_weight", "class")
  
  #########converte
  tableFeatures$GC <- as.numeric(as.character(tableFeatures$GC))
  tableFeatures$GC1 <- as.numeric(as.character(tableFeatures$GC1))
  tableFeatures$GC2 <- as.numeric(as.character(tableFeatures$GC2))
  tableFeatures$GC3 <- as.numeric(as.character(tableFeatures$GC3))
  tableFeatures$length <- as.numeric(as.character(tableFeatures$length))
  tableFeatures$dimer <- as.numeric(as.character(tableFeatures$dimer))
  tableFeatures$trimer <- as.numeric(as.character(tableFeatures$trimer))
  tableFeatures$tetramer <- as.numeric(as.character(tableFeatures$tetramer))
  tableFeatures$pentamer <- as.numeric(as.character(tableFeatures$pentamer))
  tableFeatures$hexamer <- as.numeric(as.character(tableFeatures$hexamer))
  tableFeatures$c_weight<- as.numeric(as.character(tableFeatures$c_weight))
  tableFeatures$class<- as.factor(tableSequences$class)
  
  return(tableFeatures)
  
}
