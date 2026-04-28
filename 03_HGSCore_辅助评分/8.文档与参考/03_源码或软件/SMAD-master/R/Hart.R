#' Hart
#' Scoring algorithm based on a hypergeometric distribution error model (Hart et al.,2007).
#'
#' @title Hart
#' @param datInput A dataframe with column names: idRun, idPrey. 
#' Each row represent one unique protein captured in one pull-down experiment.
#'
#' @return A dataframe consists of pairwise combindation of preys identified in 
#' the input with Hart scores indicating interacting probabilities computed from 
#' negative log transformed Hypergeometric test P-values.
#'
#' @author Qingzhou Zhang, \email{zqzneptune@hotmail.com}
#' @references Hart, G. Traver, Insuk Lee, and Edward M. Marcotte. 
#' 'A high-accuracy consensus map of yeast protein complexes reveals modular 
#' nature of gene essentiality.' BMC bioinformatics 8.1 (2007): 236.
#' \url{https://doi.org/10.1186/1471-2105-8-236}

#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr n
#' @importFrom tidyr separate
#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#' @importFrom stats phyper
#' @importFrom RcppAlgos comboGeneral
#' @useDynLib SMAD
#' @exportPattern '^[[:alpha:]]+'
#' @export
#' @examples
#' data(TestDatInput)
#' datScore <- Hart(TestDatInput)
#' head(datScore)


Hart <- function(datInput) {
    colInput <-
        c("idRun", "idPrey")
    
    if(!is.data.frame(datInput)){
        stop("Input data should be data.frame")
    }
    
    if(!all(colInput %in% colnames(datInput))){
        missingCol <-
            setdiff(colInput, 
                    colnames(datInput)[match(colInput, colnames(datInput))])
        stop("Input data missing: ", paste(missingCol, collapse = ", "))
    }
    
    
    . <- NULL
    idRun <- NULL
    idPrey <- NULL
    Tn <- NULL
    UniprotID <- NULL
    tnA <- NULL
    tnB <- NULL
    NMinTn <- NULL
    ppiTN <- NULL
    Freq <- NULL
    TnA <- NULL
    totTn <- NULL
    TnB <- NULL
    
    datCnt <-
        unique(datInput[, c("idRun", "idPrey")])
    tnInteractor <- 
        datCnt %>% 
        group_by(`idPrey`) %>% 
        summarise(`Tn` = n())
    tnA <-
        tnInteractor
    colnames(tnA) <-
        c("InteractorA", "TnA")
    tnB <-
        tnInteractor
    colnames(tnB) <-
        c("InteractorB", "TnB")
    listCnt <-
        unstack(datCnt[, c("idPrey", "idRun")])
    listPPI <-
        lapply(listCnt, function(prt){
            if(length(prt) > 1){
                sPrt <- 
                    sort(prt)
                ps <-
                    comboGeneral(sPrt, 2)
                return(paste(ps[, 1], ps[, 2], sep = "~"))
            }else{
                return(NA)
            }
        })
    datPPI <-
        as.data.frame(table(unlist(listPPI)), stringsAsFactors = FALSE)
    colnames(datPPI) <-
        c("PPI", "Freq")
    scorePPI <- 
        datPPI %>% 
        separate(data = ., col = "PPI",
                 into = c("InteractorA", "InteractorB"),
                 sep = "~", remove = FALSE) %>% 
        left_join(., tnA, by = "InteractorA") %>% 
        left_join(., tnB, by = "InteractorB") %>% 
        mutate(`totTn` = sum(tnInteractor$Tn)) %>% 
        mutate(`Hart` = -phyper(`Freq`, `TnA`, `totTn` - `TnB`, 
                                `TnB`, lower.tail = FALSE, log.p = TRUE)) %>% 
        filter(!is.infinite(`Hart`))
    return(scorePPI)
}
