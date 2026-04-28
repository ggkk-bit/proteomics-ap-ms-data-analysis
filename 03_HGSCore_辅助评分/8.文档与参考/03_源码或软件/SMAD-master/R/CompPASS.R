#' CompPASS
#' Comparative Proteomic Analysis Software Suite (CompPASS) is based on spoke
#' model. This algorithm was developed by Dr. Mathew Sowa for defining the
#' human deubiquitinating enzyme interaction landscape
#' (Sowa, Mathew E., et al., 2009). The implementation of this algorithm was
#' inspired by Dr. Sowa's online tutorial
#' (\url{http://besra.hms.harvard.edu/ipmsmsdbs/cgi-bin/tutorial.cgi}).
#' The output includes Z-score, S-score, D-score and WD-score.
#' This function also computes entropy and normalized 
#' WD-score. The source  code for this function was based on the source code. 
#'  \url{https://github.com/dnusinow/cRomppass}
#'
#'
#' @title CompPASS
#' @param datInput A dataframe with column names: idRun, idBait, idPrey, 
#' countPrey. Each row represent one unique protein captured in one 
#' pull-down experiment
#' @return A data frame consists of unique bait-prey pairs with 
#' Z-score, S-score,D-score and WD-score indicating interacting probabilities.

#' @author Qingzhou Zhang, \email{zqzneptune@hotmail.com}
#' 
#' @references Sowa, Mathew E., et al. "Defining the human deubiquitinating
#' enzyme interaction landscape." Cell 138.2 (2009): 389-403.
#' \url{https://doi.org/10.1016/j.cell.2009.04.042}
#' @references Huttlin, Edward L., et al. "The BioPlex network: a systematic 
#' exploration of the human interactome." Cell 162.2 (2015): 425-440. 
#' \url{https://doi.org/10.1016/j.cell.2015.06.043}
#' @references Huttlin, Edward L., et al. "Architecture of the human 
#' interactome defines protein communities and disease networks." 
#' Nature 545.7655 (2017): 505. 
#' \url{https://www.nature.com/articles/nature22366}

#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr n
#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#' @importFrom stats quantile
#' @export
#' @examples
#' data(TestDatInput)
#' datScore <- CompPASS(TestDatInput)
#' head(datScore)
CompPASS <- function(datInput){
    colInput <-
        c("idRun", "idBait", "idPrey", "countPrey")
    
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
    idPrey <- NULL
    idBait <- NULL
    idRun <- NULL
    countPrey <- NULL
    MaxTSC <- NULL
    nBait <- NULL
    f_sum <- NULL
    AvePSM <- NULL
    MeanDiff <- NULL
    Mean <- NULL
    SD <- NULL
    WD_inner <- NULL
    WD_raw <- NULL
    WD_raw.factor <- NULL
    
    # Use total number of baits instead of actual AP-MS runs as
    # total number of experiments.
    # Multiple runs for the same bait were considered as replicated
    k <-
        length(unique(datInput$idBait))
    # f_sum: number of runs capturing the same prey
    f <-
        unique(datInput[, c("idBait", "idPrey")]) %>%
        group_by(`idPrey`) %>%
        summarise(`f_sum` = n())
    # p: number of replicates the bait-prey pair captured
    p <-
        datInput[, c("idBait", "idPrey")] %>%
        group_by(`idBait`, `idPrey`) %>%
        summarise(`p` = n()) %>%
        mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
    # Using MAXIMUM spectral count to compute Entropy
    e <-
        datInput %>%
        group_by(`idBait`, `idPrey`, `idRun`) %>%
        summarise(`MaxTSC` = max(`countPrey`)) %>% 
        mutate(`p` = (`MaxTSC` + 1/length(`MaxTSC`))/(sum(`MaxTSC`) + 1)) %>% 
        mutate(`Entropy` = sum(unlist(lapply(`p`, function(x) {
            -1*x*log(x, 2) 
        })))) %>%
        mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
    # AvePSM = average speactral counts across replicates.
    stats <-
        datInput %>%
        group_by(`idBait`, `idPrey`) %>%
        summarise(`AvePSM` = mean(`countPrey`)) %>%
        mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
    dat <-
        stats %>%
        left_join(., unique(e[, c("BP", "Entropy")]), by = "BP") %>%
        left_join(., f, by = "idPrey") %>%
        left_join(., p[, c("BP", "p")], by = "BP") %>%
        mutate(`k` = k)
    
    preyWithBaitOnly <-
        stats %>%
        group_by(`idPrey`) %>%
        summarise(`nBait` = n()) %>%
        filter(`nBait` == 1) %>%
        .$idPrey
    dat <-
        dat %>%
        mutate(`f_num_no_prey` =
                ifelse((`idBait` == `idPrey`)&(`idPrey` %in% preyWithBaitOnly),
                            k, k-`f_sum`))
    df_num_no_prey <-
        unique(dat[, c("idPrey", "f_num_no_prey")])
    f_num_no_prey <-
        df_num_no_prey$f_num_no_prey
    names(f_num_no_prey) <-
        df_num_no_prey$idPrey
    statsM <-
        stats %>%
        mutate(`AvePSM` = ifelse(`idBait` != `idPrey`,
                                    `AvePSM`,
                                    ifelse(`idPrey` %in% preyWithBaitOnly, 
                                        `AvePSM`, NA)))
    statsMatrix <-
        spread(`statsM`[, c("idBait", "idPrey", "AvePSM")], `idBait`, `AvePSM`)
    m <-
        as.matrix(statsMatrix[, -1])
    rownames(m) <-
        statsMatrix$idPrey 
    prey.mean <-
        rowSums(m, na.rm = TRUE) / k
    
    prey.stats <-
        data.frame(`idPrey` = names(prey.mean),
                    `Mean` = prey.mean,
                    `MeanDiff` = rowSums((m - prey.mean)^2, na.rm = TRUE),
                    `f_num_no_prey` = f_num_no_prey[names(prey.mean)],
                    stringsAsFactors = FALSE) %>%
        mutate(`SD` =  
                sqrt((`MeanDiff` + ((`Mean`^2) * (`f_num_no_prey`)))/(k - 1)))
    avePSM <-
        left_join(`dat`, prey.stats[, c("idPrey", "Mean", "SD")], by = "idPrey")
    output <-
        avePSM %>%
        mutate(`scoreZ` = (`AvePSM` - `Mean`) / (`SD`)) %>%
        mutate(`scoreS` = sqrt((`AvePSM`) * (`k`) / (`f_sum`))) %>%
        mutate(`scoreD` = sqrt((`AvePSM`) * (((`k`) / (`f_sum`))^`p`))) %>%
        mutate(`WD_inner` = (`k` / `f_sum`) * (`SD` / `Mean`)) %>%
        mutate(`WD_raw` = sqrt(`AvePSM` * (`WD_inner`^`p`)))
    
    # Weighted WD score
    i = !(is.na(output$WD_raw) | is.nan(output$WD_raw))
    WD_raw <-
        output$WD_raw[i]
    WD_raw.factor <-
        unname(quantile(WD_raw, 0.98)[1])
    
    output[i, "scoreWD"] <-
        WD_raw / WD_raw.factor
    
    output <-
        as.data.frame(output[, c("idBait", "idPrey", "AvePSM",
                                 "scoreZ", "scoreS", "scoreD",
                                    "Entropy", "scoreWD")])
    return(output)
}
