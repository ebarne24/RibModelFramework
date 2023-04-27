Std <- c("A", "A", "C", "D", "E", "F", "G", "G", "H", "I", "I", "K", "L", "L", "M", "N", "P", "P", "Q", "R", "R", "S", "S", "T", "T", "V", "V", "W","X", "Y")

AA <- c("AR", "AY", "CY", "DY", "ER", "FY", "GR", "GY", "HY", "IR", "IY", "KR", "LR", "LY", "MR", "NY", "PR", "PY", "QR", "RR", "RY", "SR", "SY", "TR", "TY", "VR", "VY", "WR","XR", "YY")

C1 <- c("GCA", "GCC", "TGC", "GAC", "GAA", "TTC", "GGA", "GGC", "CAC", "ATA", "ATC", "AAA", "TTA", "CTC", "ATG", "AAC", "CCA", "CCC", "CAA", "AGA", "CGC", "TCA", "TCC", "ACA", "ACC", "GTA", "GTC", "TGG", "TAC")

C2 <- c("GCG", "GCT", "TGT", "GAT", "GAG", "TTT", "GGG", "GGT", "CAT", NA, "ATT", "AAG", "TTG", "CTT", NA, "AAT", "CCG", "CCT", "CAG", "AGG", "CGT", "TCG", "TCT", "ACG", "ACT", "GTG", "GTT", NA, "TAT")

C3 <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "CTA", NA, NA, NA, NA, NA, NA, "CGA", NA, NA, NA, NA, NA, NA, NA, NA, NA)

C4 <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "CTG", NA, NA, NA, NA, NA, NA, "CGG", NA, NA, NA, NA, NA, NA, NA, NA, NA)


#make data frame from this data
AATablealt <- data.frame(Std, AA, C1, C2, C3, C4)

#tidyAATablealt <- do.call("cbind", lapply(AATablealt, function(x) ts(na.omit(x))))[TRUE, ]
