#' PE
#' Incorporated both spoke and matrix model.
#'
#' @title PE
#' @param datInput A dataframe with column names: idRun, idBait, idPrey, 
#' Each row represent one unique protein captured in one pull-down 
#' experiment.
#' @param rBait The value of the 'r' parameter as desribed in the publication.
#'
#' @param cntPseudo The value of the 'pseudo count' parameter.
#' @return A dataframe consists of protein-protein interactions from both Spoke 
#' and Matrix model.
#'
#' @author Qingzhou Zhang, \email{zqzneptune@hotmail.com}
#' @references Collins, Sean R., et al. "Toward a comprehensive atlas of the 
#' physical interactome of Saccharomyces cerevisiae." 
#' Molecular & Cellular Proteomics 6.3 (2007): 439-450.
#' \url{https://doi.org/10.1074/mcp.M600381-MCP200}
#' 
#' @importFrom tidyr separate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr n
#' @importFrom magrittr %>%
#' @importFrom utils unstack
#' @importFrom RcppAlgos comboGeneral
#' @useDynLib SMAD
#' @exportPattern '^[[:alpha:]]+'
#' @export
#' @examples
#' data(TestDatInput)
#' datScore <- PE(TestDatInput, 0.37, 1)
#' head(datScore)

PE <- function(datInput, rBait = 0.37, cntPseudo = 1){
    colInput <-
        c("idRun", "idBait", "idPrey")
    
    if(!is.data.frame(datInput)){
        stop("Input data should be data.frame")
    }
    
    if(!all(colInput %in% colnames(datInput))){
        missingCol <-
            setdiff(colInput, 
                    colnames(datInput)[match(colInput, colnames(datInput))])
        stop("Input data missing: ", paste(missingCol, collapse = ", "))
    }
    # .getPPI <- function(prts){
    #     if(length(prts) < 2){
    #         return(NA)
    #     }else{
    #         pps <-
    #             combn(prts, 2)
    #         s <-
    #             apply(pps, 2, sort)
    #         return(paste(s[1, ], s[2, ], sep = "~"))
    #     }
    # }
    . <- NULL
    BP <- NULL
    Freq <- NULL
    cntBP <- NULL
    cntBait <- NULL
    cntPry <- NULL
    f <- NULL
    idBait <- NULL
    idPrey <- NULL
    idRun <- NULL
    nonBP <- NULL
    obsPry <- NULL
    p <- NULL
    sSum <- NULL
    totPPS <- NULL
    datInput <-
        datInput[datInput$idBait != datInput$idPrey, ]
    # Generating Bair-idPrey space
    s <-
        apply(datInput[, c("idBait", "idPrey")], 1, sort)
    pkBP <-
        unique(paste(s[1, ], s[2, ], sep = "~"))
    
    # Generating idPrey-idPrey space
    pkPrys <-
        unstack(datInput[, c("idPrey", "idRun")])
    pkPPS <-
        lapply(pkPrys, function(prts){
            if(length(prts) < 2){
                return(NA)
            }else{
                sPrts <-
                    sort(prts)
                pps <-
                    comboGeneral(sPrts, 2)
                return(paste(pps[, 1], pps[, 2], sep = "~"))
            }
        })
        # bplapply(pkPrys, FUN = .getPPI, BPPARAM = MulticoreParam(workers = 8))
    
    pkPP <-
        unique(unlist(pkPPS))
    # Final PPI space
    ppsOut <-
        data.frame(`PPI` = unique(c(pkBP, pkPP)),
                   stringsAsFactors = FALSE) %>% 
        separate(data = ., col = "PPI",
                 into = c("InteractorA", "InteractorB"),
                 sep = "~",
                 remove = FALSE)
    ppsOut[, "BP"] <-
        paste(ppsOut$InteractorA, ppsOut$InteractorB, sep = ":")
    ppsOut[, "PB"] <-
        paste(ppsOut$InteractorB, ppsOut$InteractorA, sep = ":")
    
    totObsPry <- length(datInput$idPrey) ### Total number of idPrey observations. (Redundant)
    nDistPry <- length(unique(datInput$idPrey)) # Total number of distinct preys.
    pkBait <-
        unique(datInput[, c("idRun", "idBait")])
    nPkBait <-
        as.data.frame(table(pkBait$idBait), stringsAsFactors = FALSE)
    colnames(nPkBait) <-
        c("idBait", "cntBait") # Bait and number of times each bait purified
    nObsPry <-
        as.data.frame(table(datInput$idPrey), stringsAsFactors = FALSE)
    colnames(nObsPry) <-
        c("idPrey", "obsPry") # idPrey and numebr of times each idPrey observed
    fPrey <-  # Frequency of idPrey captured across all purifications as in (Eq.9)
        nObsPry %>% 
        mutate(`f` = (`obsPry` + cntPseudo)/(totObsPry + (nDistPry * cntPseudo)))
    
    cPkPry <- # Number preys in each purifications, replicates calculated individually
        datInput %>% 
        group_by(`idBait`, `idRun`) %>% 
        summarise(`cntPry` = n())
    pkBP <- # Number of bait-idPrey pairs observed
        datInput %>% 
        group_by(`idBait`, `idPrey`) %>% 
        summarise(`totBP` = n()) %>% 
        mutate(`BP` = paste(`idBait`, `idPrey`, sep = ":"))
    dS <- # S for each idRun-bait-idPrey 
        datInput %>% 
        left_join(., nPkBait, by = "idBait") %>% 
        left_join(., fPrey[, c("idPrey", "f")], by = "idPrey") %>% 
        left_join(., cPkPry[, c("idRun", "cntPry")], by = "idRun") %>% 
        mutate(`p` = 1 - exp((-1) * `f` * `cntPry` * `cntBait`)) %>% #(Eq.8)
        mutate(`s` = log10((rBait + (1 - rBait) * `p`)/`p`)) %>% #(Eq.6)
        mutate(`BP` = paste(`idBait`, `idPrey`, sep = ":"))
    pkS <- # Sum across bait-idPrey captured
        dS[, c("BP", "s")] %>% 
        group_by(`BP`) %>% 
        mutate(`cntBP` = n(), `sSum` = sum(`s`))
    pkSU <-
        unique(pkS[, c("BP", "cntBP", "sSum")])
    sNoPrey <-
        log10(1 - rBait) # (Eq.7)
    
    spacePk <-
        expand.grid(`idBait` = unique(datInput$idBait),
                    `idPrey` = unique(datInput$idPrey),
                    stringsAsFactors = FALSE) %>% 
        filter(`idBait` != `idPrey`) %>% 
        mutate(`BP` = paste(`idBait`, `idPrey`, sep = ":")) %>% 
        left_join(., pkSU, by = "BP") %>% 
        left_join(., nPkBait, by = "idBait")
    spacePk[is.na(spacePk)] <- 0
    finalS <-
        spacePk %>% 
        mutate(`nonBP` = `cntBait` - `cntBP`) %>%
        mutate(`finalS` = `sSum` + `nonBP` * sNoPrey)
    spokeBP <-
        finalS[, c("BP", "finalS")]
    colnames(spokeBP) <-
        c("BP", "spokeBP")
    spokePB <-
        finalS[, c("BP", "finalS")]
    colnames(spokePB) <-
        c("PB", "spokePB")
    totPkList <-
        unstack(datInput[, c("idPrey", "idRun")])
    f_1 <- fPrey[, c("idPrey", "f")]
    colnames(f_1) <- c("pry1", "f_1")
    f_2 <- fPrey[, c("idPrey", "f")]
    colnames(f_2) <- c("pry2", "f_2")
    
    totPP <-
        as.data.frame(table(unlist(pkPPS)), stringsAsFactors = FALSE)
    colnames(totPP) <-
        c("PPI", "Freq")
    
    totPP[, "totPPS"] <-
        length(unlist(pkPPS))
    pkM <-
        separate(data = totPP,
                 col = "PPI",
                 into = c("pry1", "pry2"),
                 sep = "~",
                 remove = FALSE) %>% 
        left_join(., f_1, by = "pry1") %>% 
        left_join(., f_2, by = "pry2") %>% 
        mutate(`p` = 1 - exp(-`f_1` * `f_2` * `totPPS`)) %>% 
        mutate(`finalM` = (log10((rBait + (1 - rBait) * `p`)/`p`)) * `Freq`)
    
    finalM <-
        unique(pkM[, c("PPI", "finalM")])
    colnames(finalM) <-
        c("PPI", "matrixPP")
    peOut <-
        ppsOut %>% 
        left_join(., spokeBP, by = c("BP")) %>% 
        left_join(., spokePB, by = c("PB")) %>% 
        left_join(., finalM, by = "PPI")
    peOut[, "PE"] <-
        rowSums(peOut[, c("spokeBP", "spokePB", "matrixPP")], na.rm = TRUE)
    return(peOut)
}
