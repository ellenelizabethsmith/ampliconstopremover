# ampliconstopremover
This package provides a function to remove sequences containing stop codons from a FASTA file. It can be useful for filtering amplicon data from functional genes after ASV or OTU creation. If the sequences contain a string of 10 "N"s, such as added when concatenating unmergeable paired reads in DADA, both sides of the gap will be resolved independently.
The tool searches within all possible reading frames to find the most likely (the assumption being that this is the reading frame which produces the fewest stop codons).

The default behaviour is to write both the "passing" and "failing" reads to new fasta files and return a new fasta object containing only sequences which did not contain a stop codon.

## Example usage
```{r}
devtools::install_github("ellenelizabethsmith/ampliconstopremover")
library(ampliconstopremover)
library(seqinr)

#load fasta file
aoa_fasta <- read.fasta("AOA.fasta")
#set writefiles = FALSE to return only the cleaned FASTA object without writing the sequences to a file.
clean_fasta <- remove_stop_codons(aoa_fasta,fileoutprefix = "AOA",writefiles = FALSE)
```
