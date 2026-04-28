#' DICE
#' The Dice coefficient is used to score the interaction affinity between two proteins. 
#'
#' @title DICE
#' @param datInput A dataframe with column names: idRun, idPrey.
#' Each row represent one unique protein captured in one pull-down 
#' experiment.
#'
#' @return A dataframe consists of pairwise combindation of preys identified in 
#' the input with DICE scores.
#'
#' @author Qingzhou Zhang, \email{zqzneptune@hotmail.com}
#' @references Bing Zhang et al., From pull-down data to protein interaction 
#' networks and complexes with biological relevance, Bioinformatics, 
#' Volume 24, Issue 7, 1 April 2008, Pages 979–986, 
#' \url{https://doi.org/10.1093/bioinformatics/btn036}

#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#' @importFrom RcppAlgos comboGeneral
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib SMAD
#' @exportPattern '^[[:alpha:]]+'
#' @export
#' @examples
#' data(TestDatInput)
#' datScore <- DICE(TestDatInput)
#' head(datScore)


DICE <- function(datInput) {
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
    countPrey <- NULL
    Tn <- NULL
    s <- NULL
    datCnt <-
        unique(datInput[, c("idRun", "idPrey")])
    datCnt[, "Tn"] <- 1
    d <- spread(datCnt[, c("idRun", "idPrey", "Tn")], # idPrey sorted already
                `idRun`, `Tn`)
    g <- as.matrix(d[, -1])
    g[is.na(g)] <- 0
    rownames(g) <- d$idPrey
    datDice <-
        .GetDICE(t(g))
    pps <- 
        comboGeneral(d$idPrey, 2)
    ppDice <- 
        datDice[lower.tri(datDice, diag = FALSE)]
    datPPI <- 
        data.frame(cbind(pps, ppDice), 
                   stringsAsFactors = FALSE)
    colnames(datPPI) <-
        c("InteractorA", "InteractorB", "DICE")
    datPPI$DICE <-
        as.numeric(datPPI$DICE)
    datPPI[, "PPI"] <-
        paste(datPPI$InteractorA, datPPI$InteractorB, sep = "~")
    return(datPPI)
}
