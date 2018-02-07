library(tidyverse)

eln_risk_caller <- function(full.sheet) {
  updated.sheet <- full.sheet %>% 
    mutate(risk = as.factor(
      #Assigns ELN prognostic risk based on decision tree
      ifelse(PML_RARA == 1, "Favorable",
             ifelse(TP53, "Adverse",
                    ifelse(RUNX1_RUNX1T1==1 | CBFB_MYH11==1, "Favorable",
                    ifelse(DEK_NUP214 == 1 | MLL_rearranged == 1 | BCR_ABL == 1 |
                                    RPN_EVI1_Inv3 == 1 | Monosomy_Deletion_5 == 1 | Monosomy_7 == 1 |
                                    Abnormal_17 == 1 | complex_karyotype == 1 | monosomal_karyotype == 1 |
                                    double_minutes == 1, "Adverse",
                    ifelse(MLLT3_KMT2A == 1, "Intermediate",
                    ifelse(other_abnormalities == 1 & !(RUNX1|ASXL1), "Intermediate",
                    ifelse(other_abnormalities == 1 & (RUNX1|ASXL1), "Adverse",
                    ifelse(FLT3_ITD == 1 & !(RUNX1|ASXL1), "Intermediate",
                    ifelse(FLT3_ITD == 1 & (RUNX1|ASXL1), "Adverse",      
                    ifelse(NPM1, "Favorable",
                    ifelse(CEBPA, "Favorable",                
                    ifelse(RUNX1 | ASXL1, "Adverse",
                      "Intermediate"))))))))))))))      
  
  risk_col <- as.factor(updated.sheet$risk)
  
  return(risk_col)
  
}
