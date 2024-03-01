library(seqinr)
library(ggplot2)

readAli<-function(chain){
  ali=as.matrix.alignment(read.alignment(chain,format="FASTA"))
  print(dim(ali))
  return(ali)
}

getCoverage<-function(ali, minCov=0.8){
  cov=apply(ali,2,f<-function(x){1-sum(x=="-"|x=="x")/length(x)})
  return(c(sum(cov>minCov))/dim(ali)[2])
  
}

combi2<-function(x){
  return(x*(x-1)/2)
}

# the ideal alignment length is the number of amino acids in the longest sequence
# warning! this works only if the match score is 1 (I should put this match ascore as argument)
getMaxScores<-function(ali, gapChar = c("-","x")){
  nbGaps = apply(ali,1,f<-function(x){sum(x==gapChar[1])})
  if(length(gapChar)>1){
    for(k in 2:length(gapChar)){
      nbGaps = nbGaps + apply(ali,1,f<-function(x){sum(x==gapChar[k])})
    }
  }
  nbaas = max(dim(ali)[[2]] - nbGaps)
  nseqs = dim(ali)[[1]]
  maxScore = nseqs * (nseqs-1) * nbaas / 2
  print(nbaas)
  return(maxScore)
}

# sum-of-pairs score computed for a MSA, with fixed match, mismatch and gap values
# the score is normalized by the score expected for an ideal alignment 
# (i.e. a MSA with only matches and of the lenght equal to the max number of aa, computed over the sequences in the MSA)
# usage cases:
# 1) % of identity can be retrieved by setting the mismatch penalty to zero
# and disregarding the gaps and Xs:  computeSumOfPairs(readAli("6CCHB"), mismatch=0,  gapChar=c("-","x"))
# 2) global quality of the MSA can be assessed by setting the gap penalty equal to the mismatch 
# 3) quality of the gaps can be assessed and compared by setting the gap penalty to some value
computeSumOfPairs<-function(ali, match = 1, mismatch = -0.5, gap = 0, gapChar = c("-","x")){
  nseqs = dim(ali)[[1]]
  npos = dim(ali)[[2]]
  Ctot = nseqs*(nseqs-1)/2
  score = 0
  for(j in 1:npos){
    t = table(ali[,j])
    # select the match/mismatch positions
    selMatches = which(names(t)!=gapChar[1]&t>1)
    if(length(gapChar)>1){
      for(k in 2:length(gapChar)){
        selMatches = intersect(selMatches, which(names(t)!=gapChar[k]&t>1))
      }
    }
    if(length(selMatches)>0){Cmatch = sum(combi2(t[selMatches]))}
    else{Cmatch = 0} 
    tgap = 0
    for(k in 1:length(gapChar)){
      if(gapChar[k]%in%names(t)){
        tgap = tgap + t[gapChar[k]]}
    }
    if(tgap>0){Cgaps = combi2(tgap) + tgap*(nseqs-tgap)}
    else{Cgaps = 0}
    score = score + Cmatch * match + Cgaps * gap + (Ctot - Cgaps - Cmatch) * mismatch 
  }
  return(score/getMaxScores(ali,gapChar))
}



gene = "ENSG00000107643"

file_path <- paste0("/Users/louiscarrel/Documents/Alignement_Project/",gene,"/New_alignement")
file_list <- list.files(path = file_path, pattern = "\\.fasta$", full.names = TRUE)

# Initialisation du tableau
results_table <- data.frame()
alis = list

for (file_path in file_list) {
  cat("Processing file:", file_path, "\n")
  
  # Extraire le nom du fichier à partir du chemin complet
  file_name <- basename(file_path)
  file_name =  sub(".*exon_(.*)\\.fasta", "\\1", file_name)
  
  # Read alignment
  ali <- readAli(file_path)
  if (dim(ali)[1] > 13){
    color <- "canonique"
  }
  else {color <- "non canonique"}
  
  # Get coverage
   cov = getCoverage(ali)
 
  
  
  
  # Compute sum-of-pairs score
  score <- computeSumOfPairs(ali)
  
  results_table[file_name, "Mean_Coverages"] <- as.numeric(cov)
  results_table[file_name, "Sum_of_pairs_score"] <- as.numeric(score)
  results_table[file_name, "Color"] <- color
  cat("\n")
}

# Création du scatterplot
ggplot(results_table, aes(x = Sum_of_pairs_score, y = Mean_Coverages, label = rownames(results_table), col = Color )) +
  geom_point() +
  geom_text(hjust = 1.2, vjust = 0) +
  labs(title = "Scatterplot des Sum-of-pairs scores",
       x = "Sum-of-pairs score",
       y = "Nombre de cov au dessus de 80% normalisé par la longueur de la séquence") +
  theme_minimal()


