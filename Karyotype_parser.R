karyotype_parse <- function(karyotype_col) {
  library(tidyverse)
  library(stringr)
  
  
  karyo.df <- karyotype_col %>% 
    data.frame("Karyotype" = ., stringsAsFactors=FALSE)
  
  karyotype_abnormalities <- function(full.sheet) {
    updated.sheet <- full.sheet %>% 
      mutate(abnormalities =
               (str_count(Karyotype, ",") +
               str_count(Karyotype, "(?=[^,]\\+[123456789XY])")) - #Add one for cases where a "+" is not preceded by a ","; catching miswritten karyotypes
               str_count(Karyotype, "/") - 1)
    return(updated.sheet)
  }
  
  search.cyto <- function(full.sheet, search.term) {
    matches <- ifelse(str_detect(full.sheet$Karyotype, search.term), 1, 0)
    return(matches)
  }
  
  search.cyto2 <- function(full.sheet, search1, search2) {
    matches <- ifelse(str_detect(full.sheet$Karyotype, search1), 1,
                      ifelse(str_detect(full.sheet$Karyotype, search2), 1, 0))
    return(matches)
  }
  
  ## Specific abnormality searches
  
  find.PML <- function(full.sheet) {
    search.result <- search.cyto(full.sheet, "t\\(15;17\\).*\\(q22.*q21")
    updated.sheet <- mutate(full.sheet, PML_RARA = search.result)
    return(updated.sheet)
  }
  
  find.RUNX1 <- function(full.sheet) {
    search.result <- search.cyto2(full.sheet, "t\\(8;21\\).*\\(q22.*q22", "t\\(21;8\\).*\\(q22.*q22")
    updated.sheet <- mutate(full.sheet, RUNX1_RUNX1T1 = search.result)
    return(updated.sheet)
  } 
  
  find.CBFB <- function(full.sheet) {
    search.result <- search.cyto2(full.sheet, "t\\(16;16\\).*\\(p13.*q22", "inv\\(16\\).*\\(p13.*q22")
    updated.sheet <- mutate(full.sheet, CBFB_MYH11 = search.result)
    return(updated.sheet)
  }
  
  find.KMT2A.MLLT3 <- function(full.sheet) {
    search.result <- search.cyto2(full.sheet, "t\\(9;11\\)\\(p22.*;q23", "t\\(11;9\\)\\(q23.*;p22")
    updated.sheet <- mutate(full.sheet, MLLT3_KMT2A = search.result)
    return(updated.sheet)
  }
  
  find.other.abnormalities <- function(full.sheet) {
    full.sheet %>% 
      mutate(other_abnormalities =
               ifelse(
                 RUNX1_RUNX1T1 + CBFB_MYH11 + MLLT3_KMT2A + 
                   Monosomy_Deletion_5 + Monosomy_7 + Abnormal_17 + complex_karyotype + monosomal_karyotype +
                   double_minutes == 0 & abnormalities > 0, 1, 0))
    
  }
  
  
  
  
  find.DEK <- function(full.sheet) {
    search.result <- search.cyto(full.sheet, "t\\(6;9\\).*\\(p23.*q34")
    updated.sheet <- mutate(full.sheet, DEK_NUP214 = search.result)
    return(updated.sheet)
  }
  
  find.KMT2A.other <- function(full.sheet) {
    search.result <- ifelse(
      #Check for 11q23 translocations when listed first
      (
        search.cyto(
          full.sheet, "t\\(11;.*\\(q23;") |
          #Check for 11q23 translocations when listed second
          (
            search.cyto(
              full.sheet, "t\\([[:digit:]]+;11\\).*\\([[:alpha:]]{1}[[:digit:]]+;q23"))) &
        #Exclude t(9;11) translocations
        (
          full.sheet$MLLT3_KMT2A == 0),
      1, 0)
    updated.sheet <- mutate(full.sheet, MLL_rearranged = search.result)
    return(updated.sheet)
  }
  
  find.bcrabl <- function(full.sheet) {
    search.result <- search.cyto(full.sheet, "t\\(9;22\\).*\\(q34.*q11")
    updated.sheet <- mutate(full.sheet, BCR_ABL = search.result)
    return(updated.sheet)
  }
  
  find.GATA2 <- function(full.sheet) {
    search.result <- search.cyto2(full.sheet, "inv\\(3.*\\(q21.*q26", "t\\(3;3\\).*\\(q21.*q26")
    updated.sheet <- mutate(full.sheet, RPN_EVI1_Inv3 = search.result)
    return(updated.sheet)
  }
  
  find.mono.del5 <- function(full.sheet) {
    search.result <- search.cyto2(full.sheet, "-5", "del\\(5\\)\\(q")
    updated.sheet <- mutate(full.sheet, Monosomy_Deletion_5 = search.result)
    return(updated.sheet)
  }
  
  
  find.mono7 <- function(full.sheet) {
    search.result <- search.cyto(full.sheet, "-7")
    updated.sheet <- mutate(full.sheet, Monosomy_7 = search.result)
    return(updated.sheet)
  }
  
  find.abn17 <- function(full.sheet) {
    search.result <- search.cyto2(full.sheet, "-17", "del\\(17\\)\\(p|t\\(17;[[:digit:]]+\\)\\(p|\\+17|i\\(17\\)\\(q|t\\([[:digit:]]+;17\\)\\([[:alpha:]]{1}[[:digit:]]+;p")
    updated.sheet <- mutate(full.sheet, Abnormal_17 = search.result)
    return(updated.sheet)

  }

  find.complex <- function(full.sheet) {
    updated.sheet <- full.sheet %>% 
    mutate(complex_karyotype =
             ifelse(
               RUNX1_RUNX1T1 + CBFB_MYH11 + MLLT3_KMT2A + 
                 MLL_rearranged + DEK_NUP214 + RPN_EVI1_Inv3 + BCR_ABL >= 1, 0,
              ifelse(
                abnormalities >= 3, 1, 0
              )))
    return(updated.sheet)
  }

  find.monosomy <- function(full.sheet) {
    updated.sheet <- full.sheet %>% 
      mutate(monosomal_karyotype = 
              ifelse(abnormalities - RUNX1_RUNX1T1 - CBFB_MYH11 <= 1, 0, #Monosomal must have at least 2 abnormalities (monosomy + one other)
                     # Core binding factor abnormalities don't count for second abn.
                      ifelse(str_detect(full.sheet$Karyotype, "-[[:digit:]]"), 1, 0)))
    return(updated.sheet)
  }
  
  find.double.minutes <- function(full.sheet) {
    search.result <- search.cyto(full.sheet, "dmin")
    updated.sheet <- mutate(full.sheet, double_minutes = search.result)
    return(updated.sheet)
  }
  
  find.normal.karyotype <- function(full.sheet){
    updated.sheet <- full.sheet %>% 
      mutate(normal_karyotype = 
               ifelse(abnormalities == 0 &
                        RUNX1_RUNX1T1 + CBFB_MYH11 + MLLT3_KMT2A + 
                        MLL_rearranged + DEK_NUP214 + RPN_EVI1_Inv3 + BCR_ABL + PML_RARA +
                        Monosomy_Deletion_5 + Monosomy_7 + Abnormal_17 + complex_karyotype + monosomal_karyotype +
                        double_minutes == 0, 1, 0
                        ))
    return(updated.sheet)
  }

  parsed.karyotypes <- karyo.df %>% 
    karyotype_abnormalities() %>% 
    find.PML() %>% 
    find.RUNX1() %>% 
    find.CBFB() %>% 
    find.KMT2A.MLLT3() %>% 
    find.DEK() %>% 
    find.KMT2A.other() %>% 
    find.bcrabl() %>% 
    find.GATA2() %>% 
    find.mono.del5() %>% 
    find.mono7() %>% 
    find.abn17() %>%
    find.complex() %>% 
    find.double.minutes() %>% 
    find.monosomy() %>% 
    find.other.abnormalities() %>% 
    find.normal.karyotype() %>% 
    select(-Karyotype)
    
  return(parsed.karyotypes)
  }
  
