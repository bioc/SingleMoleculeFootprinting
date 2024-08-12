#' Design states for single TF case
#'
#' @return list of states
#' 
#' @export
#'
OneTFstates = function(){

  allPos=expand.grid(c(0,1),c(0,1),c(0,1))
  patternStrings=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
  #using only 'pure' states
  states=list(
    bound=patternStrings[6],
    accessible=patternStrings[8],
    closed=patternStrings[c(1,2,5)],
    unassigned=patternStrings[!seq_along(patternStrings)%in%c(1,2,5,6,8)]
  )

  return(states)

}

#' Design states for TF pair case
#'
#' @return list of states
#' 
#' @export
#'
TFpairStates = function(){

  #define states for each factor separately
  allPos=expand.grid(c(0,1),c(0,1))
  TF1=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
  names(TF1)=c('nucleosome','unassigned','bound','accessible')
  TF2=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
  names(TF2)=c('nucleosome','bound','unassigned','accessible')

  #create a combined state factor
  combined_states=unlist(lapply(seq_along(TF1),function(i){
    lapply(seq_along(TF2),function(j){
      paste(TF1[i],TF2[j],sep='')
    })}))
  names(combined_states)=unlist(lapply(seq_along(TF1),function(i){
    lapply(seq_along(TF2),function(j){
      paste(names(TF1)[i],names(TF2)[j],sep='_')
    })}))
  combined_statesF=as.factor(unlist(lapply(seq_along(combined_states),function(i){rep(names(combined_states[i]),length(combined_states[[i]]))}))[order(unlist(combined_states))])
  combined_statesF=factor(combined_statesF,levels=names(combined_states))
  stateM=cbind(string.split(as.character(combined_statesF),'_',1),string.split(as.character(combined_statesF),'_',2))
  grouped_states=list(
    combined_states[stateM[,1]=='bound'& stateM[,2]=='bound'],
    combined_states[stateM[,1]=='bound'& !stateM[,2]=='bound' & !stateM[,2]=='nucleosome'],
    combined_states[!stateM[,1]=='bound'& !stateM[,1]=='nucleosome' &stateM[,2]=='bound'],
    combined_states[(!stateM[,1]=='bound'& !stateM[,2]=='bound')|(stateM[,1]=='bound'& stateM[,2]=='nucleosome')|(stateM[,2]=='bound'& stateM[,1]=='nucleosome')])
  grouped_states = rev(grouped_states)
  names(grouped_states) = c("misc", "misc_bound", "bound_misc", "bound_bound")

  return(grouped_states)
}

SortReads_internal = function(SortedReads, SM_mat, isClusters){ # this orders readIDs based on arbotrary states order

  read_sort=SortedReads#[[1]]
  read_sort=read_sort[!is.na(names(read_sort))]
  if(isClusters){
    states = TFpairStates()
    stN=string.split(names(read_sort),'_',4)
  } else {
    states = OneTFstates()
    stN=names(read_sort)#unlist(lapply(names(read_sort), function(x){strsplit(x, "_")[[1]][3]}))
  }
  names(read_sort)=stN
  read_sort.e=vector("list", length(unlist(states)))
  names(read_sort.e)=unlist(states)
  read_sort.e[stN]=read_sort
  read_sort=read_sort.e
  read_sort=lapply(seq_along(read_sort),function(i){read_sort[[i]][read_sort[[i]]%in%rownames(SM_mat)]}) # I don't get it, these should already be the same
  names(read_sort)=unlist(states)
  read_sort=read_sort[(unlist(states))]
  read_sort
  return(read_sort)

}

#' Calculate colMeans after dropping zeros
#' 
#' @param MethSM one single molecule sparse matrix
#' 
#' @import Matrix
#'  
#' @return colMeans (N.b. this is +1 based)
#'  
colMeans_drop0 <- function (MethSM) {
  nnz_per_col <- diff(MethSM@p)
  # nnz_per_col[nnz_per_col == 0] <- 1  ## just avoid doing 0 / 0
  return(Matrix::colSums(MethSM) / nnz_per_col)
}

#' Calculate rowMeans after dropping zeros
#' 
#' @param MethSM one single molecule sparse matrix
#' 
#' @import Matrix
#'  
#' @return rowMeans (N.b. this is +1 based)
#'  
rowMeans_drop0 <- function (MethSM) {
  RowInd <- MethSM@i + 1
  nnz_per_row <- tabulate(RowInd)
  if (length(nnz_per_row) < MethSM@Dim[1]) {
    nnz_per_row = c(nnz_per_row, rep(0, MethSM@Dim[1] - length(nnz_per_row)))
  }
  # nnz_per_row[nnz_per_row == 0] <- 1  ## just avoid doing 0 / 0
  return(Matrix::rowSums(MethSM) / nnz_per_row)
}

#' Recalculate *_T and *_M values in MethGR object after filtering reads e.g. for conversion rate
#' 
#' @param MethGR GRanges object of methylation call
#' @param MethSM Single Molecule methylation matrix
#' @param MethSM_filtered Single Molecule methylation matrix after filtering reads
#' @param sampleIndex index for sample to treat. It serves as a correspondence between the index of the SM matrix and the order samples appear in the elementMetadata() columns
#' 
#' @import Matrix
#' @import GenomicRanges
#' 
#' @return MethGR with recalculated counts
#' 
filter_reads_from_MethGR = function(MethGR, MethSM, MethSM_filtered, sampleIndex){
  
  DiscardedReads = MethSM@Dimnames[[1]][!(MethSM@Dimnames[[1]] %in% MethSM_filtered@Dimnames[[1]])]
  DiscardedReads = ifelse(length(DiscardedReads) == 0, 0, DiscardedReads)
  AffectedCytosines = MethSM[DiscardedReads,,drop=FALSE]@Dimnames[[2]][colSums(MethSM[DiscardedReads,,drop=FALSE]) > 0]
  AffectedCytosines = ifelse(length(AffectedCytosines) == 0, 0, AffectedCytosines)
  T_counts = diff(MethSM_filtered[,AffectedCytosines,drop=FALSE]@p)
  M_counts = colMeans_drop0(MethSM_filtered[,AffectedCytosines,drop=FALSE]) - 1
  elementMetadata(MethGR)[start(MethGR) %in% AffectedCytosines, grep("_T$", colnames(elementMetadata(MethGR)), value=TRUE)[sampleIndex]] = T_counts
  elementMetadata(MethGR)[start(MethGR) %in% AffectedCytosines, grep("_M$", colnames(elementMetadata(MethGR)), value=TRUE)[sampleIndex]] = M_counts
  return(MethGR)
  
}
