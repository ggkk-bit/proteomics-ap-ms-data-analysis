# Statistical Modelling of AP-MS Data (SMAD)

This R package implements statistical modelling of affinity purification–mass spectrometry (AP-MS) data to compute confidence scores to identify *bona fide* protein-protein interactions (PPI).

## Installation

The development version can be installed through github:
```{r}
 devtools::install_github(repo="zqzneptune/SMAD")
 library(SMAD)
```
## Input Data
A demo data.frame was provided as a hint how the input data should strcutured in order to run the scoring functions:

```{r}
data(TestDatInput)
colnames(TestDataInput)

[1] "idRun" "idBait" "idPrey" "countPrey" "lenPrey" 
```

|idRun|idBait|idPrey|countPrey|lenPrey|
|-----|:----:|:----:|:-------:|-------|
|Unique ID of one AP-MS run|Bait ID|Prey ID|Prey peptide count|Protein sequence length of the prey|

In case of duplcates, a suffix or prefix of e.g. "A", "B" could be added to **idRun** in order to make **"idRun-idBait"** combination unique to each replicate.

## Run scoring

### 1. CompPASS

Comparative Proteomic Analysis Software Suite (CompPASS) is based on spoke model. This algorithm was developed by Dr. Mathew Sowa for defining the human deubiquitinating enzyme interaction landscape [(Sowa, Mathew E., et al., 2009)][1]. The implementation of this algorithm was inspired by Dr. Sowa's [online tutorial][2]. The output includes Z-score, S-score, D-score and WD-score. In its implementation in BioPlex 1.0 [(Huttlin, Edward L., et al., 2015)][3] and 
BioPlex 2.0 [(Huttlin, Edward L., et al., 2017)][4], a naive 
Bayes classifier that learns to distinguish true interacting proteins from 
non-specific background and false positive identifications was included in the 
compPASS pipline. This function was optimized from the [source code][5].

The input data.frame, *datInput*, should include:**idRun**, **idBait**, **idPrey** and **countPrey**.

```{r}
datScore <- CompPASS(datInput)
```

### 2. DICE

The Dice coefficient is used to score the interaction scores across prey pair-wise combinations, which was proposed by [(Bing Zhang et al., 2008)][9]

The input data.frame, *datInput*, should include:**idRun** and **idPrey**.

```{r}
datScore <- DICE(datInput)
```

### 3. Hart

Hart scoring algorithm is based on a hypergeometric distribution error model [(Hart et al., 2007)][6].

The input data.frame, *datInput*, should include:**idRun** and **idPrey**.

```{r}
datScore <- Hart(datInput)
```


### 4. HGScore

HGScore algorithm is based on a hypergeometric distribution error model [(Hart et al., 2007)][6] with incorporation of NSAF [(Zybailov, Boris, et al., 2006)][7]. This algorithm was first introduced to predict the protein complex network of Drosophila melanogaster [(Guruharsha, K. G., et al., 2011)][8]. This scoring algorithm was based on matrix model. 

The input data.frame, *datInput*, should include:**idRun**, **idPrey**, **countPrey** and **lenPrey**.

```{r}
datScore <- HG(datInput)
```

### 5. PE
PE incorporated both spoke and matrix model as repored in [(Sean R. Collins, et al., 2007)][10].

The input data.frame, *datInput*, should include:**idRun**, **idBait** and **idPrey**.

```{r}
datScore <- PE(datInput)
```



## License

MIT @ Qingzhou Zhang

[1]: https://doi.org/10.1016/j.cell.2009.04.042
[2]: http://besra.hms.harvard.edu/ipmsmsdbs/cgi-bin/tutorial.cgi
[3]: https://doi.org/10.1016/j.cell.2015.06.043
[4]: https://www.nature.com/articles/nature22366
[5]: https://github.com/dnusinow/cRomppass
[6]: https://doi.org/10.1186/1471-2105-8-236
[7]: https://doi.org/10.1021/pr060161n
[8]: https://doi.org/10.1016/j.cell.2011.08.047
[9]: https://doi.org/10.1093/bioinformatics/btn036
[10]: https://doi.org/10.1074/mcp.M600381-MCP200