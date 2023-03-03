aminoAcids()
#aminoAcids is a function

#write new amino acid function to be able to change number of aa included

Std <- c("A", "A","A", "A","C","C", "D", "D", "E", "E", "F", "F", "G", "G", "G", "G", "H", "H", "I", "I", "I", "K", "K", "L", "L", "L", "L", "L", "L", "M", "N", "N", "P", "P", "P", "P", "Q", "Q", "R", "R", "R", "R", "R", "R", "S", "S", "S", "S", "T", "T", "T", "T","W", "V", "V", "V", "V", "Y", "Y")

AA <- c("AR", "AR", "AY", "AY", "CY", "CY", "DY", "DY", "ER", "ER", "FY", "FY", "GR", "GR", "GY", "GY", "HY", "HY", "IR", "IY", "IY", "KR", "KR", "LR", "LR", "LR2", "LR2", "LY", "LY", "MR", "NY", "NY", "PR", "PR", "PY", "PY", "QR", "QR", "RR", "RR", "RR2", "RR2", "RY", "RY", "SR", "SR", "SY", "SY", "TR", "TR", "TY", "TY","WR","VR", "VR", "VY", "VY", "YY", "YY")

Codons <- c("GCA", "GCG","GCC", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG", "TTC", "TTT", "GGA", "GGG", "GGC", "GGT", "CAC", "CAT", "ATA", "ATC", "ATT", "AAA", "AAG", "TTA", "TTG", "CTA", "CTG", "CTC", "CTT", "ATG", "AAC", "AAT", "CCA", "CCG", "CCC", "CCT", "CAA", "CAG", "CGA", "CGG", "AGA", "AGG", "CGC", "CGT", "TCA", "TCG", "TCC", "TCT", "ACA", "ACG", "ACC", "ACT","TGG","GTA", "GTG", "GTC", "GTT", "TAC", "TAT")


length(Std)
length(AA)
length(Codons)

AATable <- data.frame(Std, AA, Codons)
