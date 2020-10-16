### Se placer dans le bon répertoire et avoir accès à df2 dans votre env
df2 
library(Biostrings)
require(MASS)
# "GGAAA" "AAAAGAAACAACAACAACAAC" rspt tête et queue de séquences (identiques)

strcmp = function(s1, s2){ # compare deux string
  if (s1 != s2) {
    return(FALSE)
  }
  return(TRUE)
}
tail_head_check = function(sequences){ # verifie juste que la queue et la tête sont présente dans tt les séquences
  verif = TRUE
  for (i in c(1:length(sequences ))) {
    seq = sequences[i]
    tete = substr(seq, start = 1, stop = 5)
    queu = substr(seq, start = 87, stop = 107)
    if (! (strcmp(tete, "GGAAA")) ) {
      verif = FALSE
    }
    if (! (strcmp(queu, "AAAAGAAACAACAACAACAAC"))) {
      verif = FALSE
    }
  }
  return(verif)
}
cut_tail_head = function( sequences){ # renvoie les séquences sans la tête et la queue
  cuted = rep(NA, length(sequences))
  for (i in c(1:length(sequences ))) {
    cuted[i] = substr(sequences[i], start = 6, stop = 86)
  }
  return(cuted)
}
align_scores = function(sequences){ # Aligne les seq 2 a 2, contruit la matrice de distance
  nseq = length(sequences)
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE, type = "RNA")
  matScores = matrix(0, nseq, nseq)
  for (i in c(1:nseq)) {
    s1 = sequences[i]
    for (j in c((i+1):nseq)) {
      if (j > nseq) {
        break
      }
      s2 = sequences[j]
      matScores[i,j] = pairwiseAlignment(s1, s2, substitutionMatrix = mat, gapOpening = 5, gapExtension = 2, scoreOnly=TRUE)
    }
  }
  return(matScores)
}
remp_zero_NA = function (matrice){ # remplace les 0 par des NA pour faciliter certaines opérations
  s = matrice
  for (i in c(1:dim(s)[1])) {
    s[c(i:dim(s)[1]),i] = NA
  }
  return(s)
}
create_myList = function(score_half,seuil){ # cree une liste en 2 dimensions avec liste[i] = vecteur de numéro de séquences qui ont un score > seuil
  # soit par exemple my_lists[[3]] =  (369, 383, 622) : 3 et 369, 3 et 383, 3 et 622 quand on les aligne ont un score > seuil
  my_lists = c()
  for (i in c(1:dim(score_half)[1])) {
    yseq = c()
    ind = 1
    for (j in c((i+1):dim(score_half)[1])) {
      if (j > dim(score_half)[1]) {
        break
      }
      if (score_half[i,j] > seuil) {
        yseq[ind] = j
        ind = ind +1
      }
    }
    my_lists = append(my_lists,list(yseq))
  }
  return(my_lists)
}
X_Y_forPlot = function(L,score_half){ # Ne fais pas ce que j'aimerai mais donne l'idee
  # permet de représenter un plot x = sequenceID, y = sequenceID ou un point traduit un score d'alignement > Score
  xseq = c()
  yseq_2 = c()
  for (i in c(1:dim(score_half)[1])) {
    if (!is.null(length(my_lists[[i]])) ) {
      xseq = append(xseq, rep(i,length(my_lists[[i]])))
      yseq_2 = append(yseq_2, my_lists[[i]])
    }
    else {
      xseq = append(xseq, i )
      yseq_2 = append(yseq_2, NA )
    }
  }
  return(list(xseq,yseq_2))
}
mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE, type = "RNA")
mat["A","G"] = -1
mat["G","A"] = -1
mat["C","U"] = -1
mat["U","C"] = -1

if (tail_head_check(df2$sequence) == TRUE) {
  cutted_sequences = cut_tail_head( df2$sequence ) # coupe les tête et queue de séquences (facultatif)
}

matScore = align_scores(cutted_sequences) # Alignement, 15 minutes environ
max(matScore) # score d alignement maximal
min(matScore) # score d alignement minimal
score_half = remp_zero_NA(matScore)
truehist(score_half, main = 'Distribution des scores d\'alignements', xlab = 'Scores')

# attention la suite c'est pas encore ça
seuil = 50
L = create_myList(score_half, seuil)
XY = X_Y_forPlot(L, score_half)
plot(XY[[1]],XY[[2]], xlab = 'sequenceID', ylab = 'sequenceID', main = 'Couples de séquences avec score d\'alignement > 50')

# Commentaires : les groupes auraient pu se faire uniquement sur cette base, en regroupant par lignes horizontales et verticales 
# comme on le voit sur le dernier plot
# j'épargne la partie hclust et PCA qui s'est avérée très fausse...
# Pour aller voir la vraie phylogénie c'est par ici : https://www.genome.jp/tools-bin/ete?id=201014012457237cdb1a2b391e1c1282e5ec0d4ebd1b9d47a89a