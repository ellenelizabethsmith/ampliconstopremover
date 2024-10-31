has_stop <- function(seq){
  translated <- seqinr::translate(str_split(seq,"")[[1]],frame = 1)
  return("*" %in% translated)
}

# Function to check for stop codons
check_stop_codons <- function(seq, start_lhs, trim_rhs) {
  
  seq <- substring(seq,(0+start_lhs),(nchar(seq)-trim_rhs))

  if(grepl("nnnnnnnnnn",seq)){ # contains Ns
    #get LHS
    lhs <- stringr::str_split(seq,"nnnnnnnnnn")[[1]][1]
    lhs <- substring(lhs,1,(nchar(lhs) -nchar(lhs) %% 3)) # make it a multiple of 3
  
    
    #get rhs
    rhs <- str_split(seq,"nnnnnnnnnn")[[1]][2]
    rhs <- substring(rhs,(1+nchar(rhs) %% 3),nchar(rhs)) #make multiple of 3

    return (has_stop(rhs) || has_stop(lhs)) 
  }else{
    return(has_stop(seq))
  } 

}

search_sequences <- function(fasta_file,start_lhs,trim_rhs){
  
  # Process each sequence in the fasta file
  sequences_with_stop <- c()
  cleaned_fasta <- ""
  bad_fasta <- ""
  badindex <- c()
  
  for (i in 1:length(fasta_file)) {
    # Extract sequence and name
    seq <- paste(as.character(fasta_file[[i]]),collapse = "")
    name <- names(fasta_file)[i]
    
    # Check for stop codons
    has_stop_codon <- check_stop_codons(seq,start_lhs,trim_rhs)
    
    # If stop codon is found, add to list
    if (has_stop_codon) {
      sequences_with_stop <- c(sequences_with_stop, name)
      bad_fasta <- c(bad_fasta, paste0(">",name, "\n",seq))
      badindex <- c(badindex,i)
    } else {
      # If no stop codon, add to cleaned fasta
      cleaned_fasta <- c(cleaned_fasta, paste0(">",name, "\n",seq))
    }
  }

  return(list(cleaned_fasta,bad_fasta))
}

pickbest <- function(fasta_file){
  
  combos <- expand.grid(c(0,1,2),c(0,1,2))
  names <- apply(combos,1,function(combo) {paste0("Forward trim: ",combo["Var1"],". Reverse trim: ",combo["Var2"])})
  maxclean <- -1
  bestres <- c()
  bestcombo <- ""
  
  print("Checking all reading frame combinations...")
  results <- apply(combos,1,function(combo){
    start_lhs <- combo["Var1"]
    trim_rhs <- combo["Var2"]
    res <- search_sequences(fasta_file,start_lhs,trim_rhs)
  })
  print("...done.")
  names(results) <- names
  
  print("Selecting best reading frame...")
  for(name in names(results)){
    res <- results[[name]]
    nclean <- length(res[[1]])
    #check if this combination has more clean seqs than the current
    if(nclean > maxclean){
      maxclean <- length(res[[1]])
      bestres <- res
      bestcombo <- name
      }
  }
  print("...done.")
  print(bestcombo)
  print(paste0("Sequences without stop codon: ",length(bestres[[1]])))
  print(paste0("Sequences with stop codon: ",length(bestres[[2]])))
  return(bestres)
}

write_fastas <- function(result, fileoutprefix = "ASVs", write=TRUE){
  cleanout <- paste0(fileoutprefix,"_pass.fasta")
  failout <- paste0(fileoutprefix,"_fail.fasta")
  
  print("Writing new files...")
  print(paste0("Writing ", cleanout))
  writeLines(result[[1]],cleanout)
  print(paste0("Writing ", failout))
  writeLines(result[[2]],failout)
  print("...done.")
  
}

#' remove_stop_codons
#'
#' @param fasta Input file in fasta format 
#' @param fileoutprefix The base name to be used in the output files 
#'
#' @return None, writes FASTA files with and without stop codons.
#' @export
#'
remove_stop_codons <- function(fasta, fileoutprefix="amplicon",writefiles=TRUE){
  bestres <- pickbest(fasta_file = fasta)
  if(writefiles) write_fastas(bestres,fileoutprefix = fileoutprefix)
  return(read.fasta(textConnection(bestres[[1]]))) #returns a fasta object
}

