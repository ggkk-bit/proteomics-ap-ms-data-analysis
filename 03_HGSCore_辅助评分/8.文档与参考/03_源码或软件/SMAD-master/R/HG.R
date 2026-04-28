#' HGScore
#' Scoring algorithm based on a hypergeometric distribution error model 
#' (Hart et al.,2007) with incorporation of NSAF (Zybailov, Boris, et al., 2006)
#' . This algorithm was first introduced to predict the protein complex network 
#' of Drosophila melanogaster (Guruharsha, K. G., et al., 2011). This scoring 
#' algorithm was based on matrix model.
#'
#' @title HGScore
#' @param datInput A dataframe with column names: idRun, idPrey, countPrey, 
#' lenPrey. Each row represent one unique protein captured in one pull-down 
#' experiment.
#'
#' @return A dataframe consists of pairwise combindation of preys identified in 
#' the input with HG scores indicating interacting probabilities computed from 
#' negative log transformed Hypergeometric test P-values.
#'
#' @author Qingzhou Zhang, \email{zqzneptune@hotmail.com}
#' @references Guruharsha, K. G., et al. 'A protein complex network of 
#' Drosophila melanogaster.' Cell 147.3 (2011): 690-703.
#' \url{https://doi.org/10.1016/j.cell.2011.08.047}
#' @references Hart, G. Traver, Insuk Lee, and Edward M. Marcotte. 
#' 'A high-accuracy consensus map of yeast protein complexes reveals modular 
#' nature of gene essentiality.' BMC bioinformatics 8.1 (2007): 236.
#' \url{https://doi.org/10.1186/1471-2105-8-236}
#' @references Zybailov, Boris, et al. 'Statistical analysis of membrane 
#' proteome expression changes in Saccharomyces c erevisiae.' Journal of 
#' proteome research 5.9 (2006): 2339-2347.
#' \url{https://doi.org/10.1021/pr060161n}

#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr n
#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @importFrom stats setNames
#' @importFrom stats phyper
#' @importFrom RcppAlgos comboGeneral
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib SMAD
#' @exportPattern '^[[:alpha:]]+'
#' @export
#' @examples
#' data(TestDatInput)
#' datScore <- HG(TestDatInput)
#' head(datScore)


HG <- function(datInput) {
    colInput <-
        c("idRun", "idPrey", "countPrey", "lenPrey")
    
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
    countPrey <- NULL
    lenPrey <- NULL
    NormalSpec <- NULL
    SumNS <- NULL
    NSAF <- NULL
    NormalNSAF <- NULL
    Tn <- NULL
    UniprotID <- NULL
    tnA <- NULL
    tnB <- NULL
    NMinTn <- NULL
    HG <- NULL
    ppiTN <- NULL
    s <- NULL
    InteractorA <- NULL
    InteractorB <- NULL
    datCnt <- 
        datInput %>% 
        mutate(`NormalSpec` = `countPrey`/`lenPrey`) %>% 
        group_by(`idRun`) %>% 
        mutate(`SumNS` = sum(`NormalSpec`)) %>% 
        mutate(`NSAF` = `NormalSpec`/`SumNS`) %>% 
        group_by(`idRun`) %>% 
        mutate(`NormalNSAF` = `NSAF`/min(`NSAF`)) %>% 
        mutate(`Tn` = as.integer(sqrt(`NormalNSAF`)))
    d <- spread(datCnt[, c("idRun", "idPrey", "Tn")], 
                `idRun`, `Tn`)
    g <- as.matrix(d[, -1])
    g[is.na(g)] <- 0
    g <- t(g)
    colnames(g) <- d$idPrey
    pps <- 
        comboGeneral(colnames(g), 2)
    PPN <- .GetPPN(g)
    CppPPN <- PPN[lower.tri(PPN, diag = FALSE)]
    datPPI <- data.frame(cbind(pps[CppPPN != 0, ], 
                                CppPPN[CppPPN != 0]), stringsAsFactors = FALSE)
    colnames(datPPI) <- 
        c("InteractorA", "InteractorB", "ppiTN")
    datPPI$ppiTN <- 
        as.numeric(datPPI$ppiTN)
    tnInteractorA <- 
        datPPI[, c("InteractorA", "ppiTN")]
    colnames(tnInteractorA) <- 
        c("UniprotID", "ppiTN")
    tnInteractorB <- 
        datPPI[, c("InteractorB", "ppiTN")]
    colnames(tnInteractorB) <- 
        c("UniprotID", "ppiTN")
    tnProtein <- 
        bind_rows(tnInteractorA, tnInteractorB) %>% 
        group_by(`UniprotID`) %>% 
        summarise(minTn = sum(`ppiTN`))
    sumMinTnInteractorA <- 
        tnProtein
    colnames(sumMinTnInteractorA) <- 
        c("InteractorA", "tnA")
    sumMinTnInteractorB <- 
        tnProtein
    colnames(sumMinTnInteractorB) <- 
        c("InteractorB", "tnB")
    scorePPI <- 
        datPPI %>% 
        left_join(., sumMinTnInteractorA, by = "InteractorA") %>% 
        left_join(., sumMinTnInteractorB, by = "InteractorB") %>% 
        mutate(`PPI` = paste(`InteractorA`, `InteractorB`, sep = "~")) %>% 
        mutate(`NMinTn` = sum(tnProtein$minTn)/2) %>% 
        mutate(`HG` = -phyper(`ppiTN`, `tnA`, `NMinTn` - `tnB`, 
                                `tnB`, lower.tail = FALSE, log.p = TRUE))
    return(scorePPI)
}
