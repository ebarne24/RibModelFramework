aminoAcids()
#aminoAcids is a function

#write new amino acid function to be able to change number of aa included

Std <- c("A", "A", "C", "D", "E", "F", "G", "G", "H", "I", "I", "K", "L", "L", "M", "N", "P", "P", "Q", "R", "R", "S", "S", "T", "T", "V", "V", "W","X", "Y")

AA <- c("AR", "AY", "CY", "DY", "ER", "FY", "GR", "GY", "HY", "IR", "IY", "KR", "LR", "LY", "MR", "NY", "PR", "PY", "QR", "RR", "RY", "SR", "SY", "TR", "TY", "VR", "VY", "WR","XR", "YY")

C1 <- c("GCA", NA, NA, NA, "GAA", NA, "GGA", NA, NA, "ATA", NA, "AAA", "CTA", NA, NA, NA, "CCA", NA, "CAA", "AGA", NA, "TCA", NA, "ACA", NA, "GTA", NA, NA, "TAA", NA)

C2 <- c("GCG", NA, NA, NA, "GAG", NA, "GGG", NA, NA, NA, NA, "AAG", "CTG", NA, "ATG", NA, "CCG", NA, "CAG", "AGG", NA, "TCG", NA, "ACG", NA, "GTG", NA, "TGG", "TAG", NA)

C3 <- c(NA, "GCC", "TGC", "GAC", NA, "TTC", NA, "GGC", "CAC", NA, "ATC", NA, NA, "CTC", NA, "AAC", NA, "CCC", NA, NA, "CGC", NA, "TCC", NA, "ACC", NA, "GTC", NA,  NA, "TAC")

C4 <- c(NA, "GCT", "TGT", "GAT", NA, "TTT", NA, "GGT", "CAT", NA, "ATT", NA, NA, "CTT", NA, "AAT", NA, "CCT", NA, NA, "CGT", NA, "TCT", NA, "ACT", NA, "GTT", NA, NA, "TAT")

C5 <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "TTA", NA, NA, NA, NA, NA, NA, "CGA", NA, NA, NA, NA, NA, NA, NA, NA, "TGA", NA)

C6 <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "TTG", NA, NA, NA, NA, NA, NA, "CGG", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

CTAG <- c("AG", "CT", "CT", "CT", "AG", "CT", "AG", "CT", "CT", "AG", "CT", "AG", "AG", "CT", "AG", "CT", "AG", "CT", "AG", "AG", "CT", "AG", "CT", "AG", "CT", "AG", "CT", "AG", "AG", "CT")



AATable <- data.frame(Std, AA, C1, C2, C3, C4, C5, C6)
amino.acids.new <- AATable$AA

#AATable <- AATable %>% 
  #mutate_if(is.character, ~replace_na(.,""))

codons <- AATable["AA", c"C1", "C2", "C3", "C4", "C5", "C6"]
aa.names <- AATable$AA

#could try subsetting data by which columns in those rows contain data and then write command to drop the rest

rows.two.purines |> slice(AATable, c(1, 5, 7, 12, 17, 19, 22, 24, 26))
for(rows.two.purines)
  rows.two.purines |> select_if(~ !any(is.na(.)))


rows.two.pyrimidines |> slice(AATable, c(2, 3, 4, 6, 8, 9, 11, 14, 16, 18, 21, 23, 25, 27, 30))
for(rows.two.pyrimidines)
  rows.two.pyrimidines |> select_if(~ !any(is.na(.)))

rows.four.purines |> slice(AATable, c(13, 20))
for(rows.four.purines)
  rows.four.purines |> select_if(~ !any(is.na(.)))

AATable[AATable$AA =="AR"]
#gave me the Std column? 


#cols_codon <- AATable[c("C1", "C2", "C3", "C4", "C5", "C6")]
 #for("AA" in AATable)
   #info_aa <- AATable$AA=="AR"
 # codon <- na.omit


# na.omit(codons)
 # sort(na.omit(codons), decreasing = TRUE)
 # if(!ROC.or.FONSE) codons <- codons[-1]
 # codons <- sort(codons, increasing)
#cannot get this code to work, should omit NAs, sort codons by decreasing, pull out last alphabetically for reference, then reorder
  
#write test case to make sure plotting works for different lists of AAs, 
  #will be inserted ~line 146 of plotTraceObject code
  if( is.null(aa.names) ) {
    aa.names <- aminoAcids()
  } else
    aa.match <- (aa.names %in% aminoAcids())
  ## test to ensure there's no aa being called that don't exist in trace
  aa.mismatch <- aa.names[!aa.match]
  if(length(aa.mismatch) > 0){
    warning("Members ", aa.mismatch, "of aa.names argument absent from trace object and will be excluded.",
            call. = TRUE, immediate. = FALSE, noBreaks. = FALSE,
            domain = NULL)
  }
  aa.names <- aa.names[aa.match]
  }

if(aa.names = (amino.acids.new) ) {
  aa.names <- AATable$AA
  warning("Alternative amino acid list selected, using amino acids with codons split by ending of purine/pyramidine.")
} else 
  aa.names <- aminoAcids()
    warning("Traditional amino acid list being used.")
}




















