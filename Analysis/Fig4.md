Figure 4 and related supplement tables and figures
================
Meng-Ching Ko (MaggieMCKO)
05/06/2023

This [R Markdown](http://rmarkdown.rstudio.com) Notebook contain codes
for reproducing Fig.4 and related supplement tables and figures of Ko et
al. (<https://www.biorxiv.org/content/10.1101/2022.06.13.495861v1>).
Data deposited on dryad: <https://doi.org/10.5061/dryad.5hqbzkh8c>

### load data

``` r
library(tidyverse) # v.2.0.0
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.2     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(RColorBrewer) # 1.1-3
library(WGCNA); # 1.71
```

    ## Loading required package: dynamicTreeCut
    ## Loading required package: fastcluster
    ## 
    ## Attaching package: 'fastcluster'
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     hclust
    ## 
    ## 
    ## 
    ## Attaching package: 'WGCNA'
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     cor

``` r
set.seed(100)

# trait 
path = paste0(getwd(), "/Data/WGCNA/Traits.tsv")
datTraits = read_tsv(path)
```

    ## Rows: 40 Columns: 25
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): sample
    ## dbl (24): T_pgml, Body Weight (g), Brain (mg), Ovi duct (mg), ImplantDuratio...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
datTraits0 = datTraits[, -1] %>% as.matrix()
row.names(datTraits0) = datTraits$sample

# exp 
path = paste0(getwd(), "/Data/WGCNA/Expression_est_perBird.tsv")
datExprori = read_tsv(path)
```

    ## New names:
    ## Rows: 40 Columns: 12361
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "\t" chr
    ## (1): sample dbl (12360): A1CF, A2M, A2ML1, A4GALT, A4GNT, AACS, AADACL4,
    ## AADACP1, AADAT,...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `PGM2` -> `PGM2...7817`
    ## • `PGM2` -> `PGM2...7818`

``` r
datExpr0 = datExprori[, -1] %>% as.matrix()
row.names(datExpr0) = datExprori$sample
```

### filter samples

``` r
gsg = goodSamplesGenes(datExpr0, verbose = 3);
```

    ##  Flagging genes and samples with too many missing values...
    ##   ..step 1

``` r
gsg$allOK
```

    ## [1] TRUE

``` r
#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the oending genes and samples from the data:
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

## smaple tree clustering 
# cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers. 
sampleTree = hclust(dist(datExpr0), method = "average");
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1, cex.axis = 1, cex.main = 1)
```

![](Fig4_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# choosing samples
MEDissThres = 60
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = MEDissThres, minSize = 1)
table(clust)
```

    ## clust
    ##  1  2 
    ## 39  1

``` r
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
deletSamples = row.names(datExpr0)[clust == 2]; deletSamples
```

    ## [1] "SDf_T3h_H_764R"

``` r
# "SDf_T3h_H_764R"

# remove traits from the outlier
ind = which(row.names(datTraits0)%in% deletSamples); ind
```

    ## [1] 9

``` r
datTraits0 = datTraits0[-ind, ]

# make a new expression data frame
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr); nGenes     # 12360
```

    ## [1] 12360

``` r
nSamples = nrow(datExpr); nSamples # 39
```

    ## [1] 39

``` r
rm(datExpr0)

## match trait data =
collectGarbage();

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits0, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(datTraits0),
                    marAll = c(1, 8, 2, 1), # addGuide = FALSE, setLayout = FALSE,
                    cex.dendroLabels = 0.7, 
                    cex.colorLabels = 0.7,
                    main = "Sample dendrogram and trait heatmap")
```

![](Fig4_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

### Network construction and module detection

``` r
## Choose a set of soft-thresholding powers
networktype = "signed hybrid"
powers = c(c(1:12), seq(from = 14, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, dataIsExpr = TRUE, powerVector = powers, networkType = networktype, verbose = 5) 
```

    ## pickSoftThreshold: will use block size 3619.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 3619 of 12360

    ## Warning: executing %dopar% sequentially: no parallel backend registered

    ##    ..working on genes 3620 through 7238 of 12360
    ##    ..working on genes 7239 through 10857 of 12360
    ##    ..working on genes 10858 through 12360 of 12360
    ##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1   0.0302 -0.149         0.0395  2510.0  2200.000   4880
    ## 2      2   0.5340 -0.560         0.7540  1340.0   784.000   3550
    ## 3      3   0.7210 -0.681         0.8150   836.0   359.000   2770
    ## 4      4   0.7720 -0.745         0.7730   570.0   186.000   2250
    ## 5      5   0.8200 -0.791         0.7830   410.0   103.000   1880
    ## 6      6   0.8390 -0.828         0.7940   306.0    59.900   1590
    ## 7      7   0.8440 -0.865         0.8060   235.0    35.900   1360
    ## 8      8   0.8380 -0.897         0.8160   184.0    22.400   1180
    ## 9      9   0.8250 -0.931         0.8220   147.0    14.300   1040
    ## 10    10   0.8160 -0.962         0.8290   119.0     9.360    914
    ## 11    11   0.7980 -0.994         0.8290    97.3     6.150    811
    ## 12    12   0.8020 -1.020         0.8450    80.4     4.140    723
    ## 13    14   0.7780 -1.080         0.8480    56.3     1.960    582
    ## 14    16   0.7700 -1.130         0.8610    40.5     0.973    475
    ## 15    18   0.7730 -1.170         0.8760    29.8     0.509    392
    ## 16    20   0.7920 -1.200         0.9020    22.3     0.275    326

``` r
# Plot the results:
plot.new()
par(mfrow = c(1,2));
cex1 = 1;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
abline(h=0.85,col="red")
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```

![](Fig4_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#  choose the power value which is the lowest power that reaches high "Scale Free Topology Model Fit,signed R^2"


## module manual calculation 
# Co-expression similarity and adjacency
softPower = 6 ; 
adjacency = adjacency(datExpr, power = softPower, type = networktype)

# export to cytoscape
# cytoscape = exportNetworkToCytoscape(adjacency,
#                          edgeFile = paste0(getwd(), "/Data/WGCNA/edge.txt"),
#                          nodeFile = paste0(getwd(), "/Data/WGCNA/node.txt"))
```

### Topological Overlap Matrix (TOM)

``` r
# Turn adjacency into topological overlap
collectGarbage()
TOM = TOMsimilarity(adjacency, TOMType = "unsigned");
```

    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.

``` r
dissTOM = 1-TOM
rm(TOM)


## Clustering using TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
```

![](Fig4_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
## Module identification using dynamic tree cut:
minModuleSize = 100;
deepSplit = 4;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid", 
                            deepSplit = deepSplit, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
```

    ##  ..cutHeight not given, setting it to 0.996  ===>  99% of the (truncated) height range in dendro.
    ##  ..done.

``` r
## Convert numeric lables into colors
# to remove color that is too light
color_to_remove = c(8, 16, 19, 27, 31, 41:44, 50, 51, 58, 87, 90, 115, 119)
colorSeq = labels2colors(setdiff(1:(length(table(dynamicMods))+length(color_to_remove)), color_to_remove))

dynamicColors = labels2colors(dynamicMods, colorSeq = colorSeq)
table(dynamicColors)
```

    ## dynamicColors
    ##        black         blue        brown         cyan        green  greenyellow 
    ##          438          891          603          190          547          372 
    ##         grey       grey60      magenta midnightblue       purple          red 
    ##           71          122          417          184          396          483 
    ##       salmon          tan    turquoise       yellow 
    ##          239          342         6485          580

``` r
# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
```

![](Fig4_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

### Merging of modules whose expression profiles are very similar

``` r
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# We choose a height cut of 0.2, corresponding to correlation of 0.8, to merge
MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
```

![](Fig4_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
plot.new()
par(cex = 1, cex.axis= 1, cex.main = 1.1, cex.lab = 1);
par(mar=c(0, 5, 2, 0)) # c(bottom, left, top, right)
par(mgp=c(2.5, 1,.5)) # axis title, axis labels and axis line.
plot(METree, main = "Clustering of module eigengenes_SDf",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
```

![](Fig4_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
```

    ##  mergeCloseModules: Merging modules whose distance is less than 0.2
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 16 module eigengenes in given set.
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 14 module eigengenes in given set.
    ##    Calculating new MEs...
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 14 module eigengenes in given set.

``` r
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
table(mergedColors)
```

    ## mergedColors
    ##        black         blue        brown         cyan        green  greenyellow 
    ##          780          891          603          190          943          372 
    ##         grey       grey60      magenta midnightblue          red       salmon 
    ##           71          122          417          184          483          239 
    ##    turquoise       yellow 
    ##         6485          580

``` r
length(unique(mergedColors))
```

    ## [1] 14

``` r
plot.new()
par(mgp=c(2.5, 1,.5)) # axis title, axis labels and axis line.
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03, marAll = c(1, 6, 3, 1), 
                    addGuide = TRUE, guideHang = 0.05)
```

![](Fig4_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
# Rename to mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50)); 
moduleLabels = match(mergedColors, colorOrder)-1;
MEs = mergedMEs

# Recalculate MEs with color labels 
nGenes = ncol(datExpr); nGenes # transcript:26061; gene: 12360
```

    ## [1] 12360

``` r
nSamples = nrow(datExpr); nSamples # 39
```

    ## [1] 39

``` r
# Calculates module eigengenes 
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)

# Put close eigenvectors next to each other
MET=orderMEs(cbind(MEs, datTraits0)) 
# Plot the relationships among the eigengenes and the trait 
plot.new()
par(cex = 1, cex.axis= 1, cex.main = 1.1, cex.lab = 1);
plotEigengeneNetworks(MET,"", marDendro=c(0,11,1,5.5), marHeatmap=c(10,12,.5,1.5), 
                      cex.lab = 1, xLabelsAngle=90, printPreservation=T)
```

![](Fig4_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

    ## Warning in mapply(textFnc, x = labPos$xMid[xTextLabInd], y = xLabYPos, labels =
    ## xLabels.show[xTextLabInd], : longer argument not a multiple of length of
    ## shorter

    ## Warning in mapply(textFnc, x = labPos$xMid[xTextLabInd], y = xLabYPos, labels =
    ## xLabels.show[xTextLabInd], : longer argument not a multiple of length of
    ## shorter

    ## Warning in mapply(textFnc, x = labPos$xMid[xTextLabInd], y = xLabYPos, labels =
    ## xLabels.show[xTextLabInd], : longer argument not a multiple of length of
    ## shorter

    ## Warning in mapply(textFnc, x = labPos$xMid[xTextLabInd], y = xLabYPos, labels =
    ## xLabels.show[xTextLabInd], : longer argument not a multiple of length of
    ## shorter

![](Fig4_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->

### Module-trait relationships

``` r
datTraits_tmp = datTraits0
datTraits_tmp1_6 = datTraits_tmp[,1:6]
datTraits_tmp7 = datTraits_tmp[,7:ncol(datTraits_tmp)]
datTraits_tmp7[is.na(datTraits_tmp7)]  = 0 
# ind_bird2 = grep(788, row.names(datTraits_tmp7))
# datTraits_tmp7[ind_bird2, 1] = 0.01529500 # songrate (not songrate_2s)
datTraits_tmp = cbind(datTraits_tmp1_6, datTraits_tmp7)


# reorder MEs
MEs0 = MEs
names(MEs0)
```

    ##  [1] "MEblue"         "MEturquoise"    "MEgreen"        "MEred"         
    ##  [5] "MEgrey60"       "MEmidnightblue" "MEsalmon"       "MEbrown"       
    ##  [9] "MEgreenyellow"  "MEblack"        "MEcyan"         "MEmagenta"     
    ## [13] "MEyellow"       "MEgrey"

``` r
sel = c("MEgrey60", "MEbrown", 
        "MEblue", "MEturquoise", "MEgreen", "MEred", 
        "MEblack", "MEcyan", "MEyellow", "MEgrey", "MEmagenta",
        "MEsalmon", "MEmidnightblue", "MEgreenyellow"
)
ind_c = sapply(sel, function(s){grep(paste0(s, "$"), names(MEs0))})
MEs = MEs0[, ind_c]

datTraits_tmp_0 = datTraits_tmp
sel = c("MedianSlope", "IqrSlope")
ind_c = sapply(sel, function(s){grep(paste0(s, "$"), colnames(datTraits_tmp_0))})
ind_c_r = c( setdiff(1:ncol(datTraits_tmp_0), ind_c), ind_c)
datTraits_tmp = datTraits_tmp_0[, ind_c_r]

moduleTraitCor0 = cor(MEs, datTraits_tmp, use = "p");
moduleTraitPvalue0 = corPvalueStudent(moduleTraitCor0, nSamples);
moduleTraitPvalue.Adj0 = apply(moduleTraitPvalue0, MARGIN = 2, FUN = function(s){p.adjust(s, "bonferroni")})
uniqModuleCol0 = substring(row.names(moduleTraitCor0), 3)
ModuleCount0 = sapply(1:length(uniqModuleCol0), FUN = function(i){
  mduleI = uniqModuleCol0[i];
  signif(length(which(mergedColors == mduleI)), 3)})
textMatrix0 = paste(signif(moduleTraitCor0, 2), "\n(", signif(moduleTraitPvalue0, 1), ")", sep = "");
dim(textMatrix0) = dim(moduleTraitCor0)
ind_non_sig0 = moduleTraitPvalue0 > 0.05
textMatrix0[ind_non_sig0] = ""
```

#### Fig. 4A plot

``` r
plot.new()
par(mar = c(7.5, 12, 3, 1), cex = 0.7);
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor0,
               xLabels = colnames(datTraits_tmp),
               yLabels = names(MEs),
               ySymbols = paste(names(MEs), ": ", ModuleCount0, sep=""),
               # ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix0,
               setStdMargins = FALSE,
               cex.text = 0.7,
               cex.lab = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
```

![](Fig4_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### geneTraitSignificance

#### Fig. 4 - Figure supplement 1A

``` r
geneTraitSignificance0 = as.data.frame(cor(datExpr, datTraits_tmp, use = "p"));
GSPvalue0 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance0), nSamples));
GSPvalue0.Adj = apply(GSPvalue0, MARGIN = 2, FUN = function(s){p.adjust(s, "fdr")})
names(geneTraitSignificance0) = paste0("GS.", colnames(datTraits_tmp));
names(GSPvalue0) = paste0("p.GS.", colnames(datTraits_tmp));
names(GSPvalue0.Adj) = paste0("padj.GS.", colnames(datTraits_tmp));
dim(geneTraitSignificance0) # 12360    24
```

    ## [1] 12360    24

``` r
dim(GSPvalue0.Adj)          # 12360    24
```

    ## [1] 12360    24

``` r
# geneTraitSignificance 
GS_dataframe = sapply(1:ncol(geneTraitSignificance0), function(i){
  c(geneTraitSignificance0[, i], GSPvalue0[, i], GSPvalue0.Adj[, i])
})
GS_dataframe = as.data.frame(t(matrix(GS_dataframe, ncol = nrow(geneTraitSignificance0), byrow = T)))
dim(GS_dataframe) #  12360    72
```

    ## [1] 12360    72

``` r
names(GS_dataframe) = paste0(c("GS.", "p.GS", "padj.GS"),  
                             rep(colnames(datTraits_tmp), each= 3))
row.names(GS_dataframe) = colnames(datExpr)
GS_dataframe = cbind(moduleColor = mergedColors, GS_dataframe)

# plot
plot.new()
par(cex = .5, cex.axis=.5, cex.main = 1, cex.lab = .5);
par(mgp=c(2.5, 1,.5)) # axis title, axis labels and axis line.
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors, numbers2colors(geneTraitSignificance0, signed = TRUE)),
                    c("Dynamic Tree Cut", "Merged dynamic", names(geneTraitSignificance0)),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = F, guideHang = 0.05, marAll = c(0, 6, 2, 0))
```

![](Fig4_files/figure-gfm/unnamed-chunk-8-1.png)<!-- --> \###
geneModuleMembership

``` r
modNames = substring(names(MEs), 3)
geneModuleMembership = cor(datExpr, MEs, use = "p");
MMPvalue = corPvalueStudent(as.matrix(geneModuleMembership), nSamples)
MMPvalue.Adj = apply(MMPvalue, MARGIN = 2, FUN = function(s){p.adjust(s, "fdr")})
names(geneModuleMembership) = paste0("MM", modNames);
names(MMPvalue) = paste0("p.MM", modNames);
names(MMPvalue.Adj) = paste0("padj.MM", modNames);

## geneModuleMembership 
dim(geneModuleMembership) # 12360    13
```

    ## [1] 12360    14

``` r
dim(MMPvalue)             # 12360    13
```

    ## [1] 12360    14

``` r
MM_dataframe = sapply(1:ncol(geneModuleMembership), function(i){
  c(geneModuleMembership[, i], MMPvalue[, i], MMPvalue.Adj[, i])
})
MM_dataframe = as.data.frame(t(matrix(MM_dataframe, ncol = nrow(geneModuleMembership), byrow = T)))
dim(MM_dataframe) #  12360    39
```

    ## [1] 12360    42

``` r
names(MM_dataframe) = paste0(c("MM", "p.MM", "padj.MM"),  rep(modNames, each = 3))
row.names(MM_dataframe) = colnames(datExpr)
MM_dataframe = cbind(moduleColor = mergedColors, MM_dataframe)

# plot
plot.new()
par(cex = .8, cex.axis=.8, cex.main = 1, cex.lab = .8);
par(mgp=c(2.5, 1,.5)) # axis title, axis labels and axis line.
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors, numbers2colors(geneModuleMembership, signed = TRUE)),
                    c("Dynamic Tree Cut", "Merged dynamic", names(geneModuleMembership)),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = F, guideHang = 0.05, marAll = c(0, 6, 2, 0))
```

![](Fig4_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

### module data

#### Supplementary Table 10

``` r
intramodularConnectivity = intramodularConnectivity(adjacency, mergedColors, scaleByMax = FALSE)
geneInfo0 = data.frame(GeneSymbol = colnames(datExpr),
                       moduleColor = mergedColors,
                       intramodularConnectivity,
                       geneTraitSignificance0)
```

#### Fig. 4 - Figure supplement 1B-E

``` r
library(grid) # v4.3.0
library(gridExtra) # v2.3
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
theme_m = theme_classic() + 
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

AveConnectivity = geneInfo0 %>%
  group_by(moduleColor) %>%
  summarise(Mean_ktotal = mean(kTotal))
ind_sort = sort(AveConnectivity$Mean_ktotal, index.return = TRUE, decreasing = TRUE)
Module_col = AveConnectivity$moduleColor[ind_sort$ix]; Module_col
```

    ##  [1] "turquoise"    "blue"         "red"          "brown"        "green"       
    ##  [6] "yellow"       "black"        "magenta"      "cyan"         "greenyellow" 
    ## [11] "midnightblue" "grey60"       "salmon"       "grey"

``` r
AveConnectivity = within(AveConnectivity, moduleColor <- factor(moduleColor, levels = Module_col))
Mean_ktotal_p = ggplot(AveConnectivity, aes(x = moduleColor, y = Mean_ktotal, col = moduleColor)) +
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.2) +
  geom_point() + 
  scale_y_continuous(name = "Mean connectivity") +
  scale_color_manual(values = Module_col) +
  ggtitle("Overall") +
  theme_m ; Mean_ktotal_p
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](Fig4_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
AveConnectivity = geneInfo0 %>%
  group_by(moduleColor) %>%
  summarise(Mean_kDiff = mean(kDiff))
ind_sort = sort(AveConnectivity$Mean_kDiff, index.return = TRUE, decreasing = TRUE)
Module_col = AveConnectivity$moduleColor[ind_sort$ix]; Module_col
```

    ##  [1] "turquoise"    "green"        "grey"         "salmon"       "black"       
    ##  [6] "yellow"       "magenta"      "greenyellow"  "midnightblue" "grey60"      
    ## [11] "red"          "brown"        "cyan"         "blue"

``` r
AveConnectivity$moduleColor <- factor(AveConnectivity$moduleColor, levels = Module_col)
Mean_kDiff_p = ggplot(AveConnectivity, aes(x = moduleColor, y = Mean_kDiff, col = moduleColor)) +
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.2) +
  geom_point() + 
  scale_y_continuous(name = "Mean connectivity") +
  scale_color_manual(values = Module_col) +
  ggtitle("Differential\n(Intra - Inter)") +
  theme_m ; Mean_kDiff_p
```

![](Fig4_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
AveConnectivity = geneInfo0 %>%
  group_by(moduleColor) %>%
  summarise(Mean_kWithin = mean(kWithin))
ind_sort = sort(AveConnectivity$Mean_kWithin, index.return = TRUE, decreasing = TRUE)
Module_col = AveConnectivity$moduleColor[ind_sort$ix]; Module_col
```

    ##  [1] "turquoise"    "blue"         "green"        "red"          "brown"       
    ##  [6] "black"        "yellow"       "magenta"      "greenyellow"  "midnightblue"
    ## [11] "grey60"       "cyan"         "salmon"       "grey"

``` r
AveConnectivity$moduleColor <- factor(AveConnectivity$moduleColor, levels = Module_col)
Mean_kWithin_p = ggplot(AveConnectivity, aes(x = moduleColor, y = Mean_kWithin, col = moduleColor)) +
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.2) +
  geom_point() + 
  scale_y_continuous(name = "Mean connectivity") +
  scale_color_manual(values = Module_col) +
  ggtitle("Intramodular") +
  theme_m ; Mean_kWithin_p
```

![](Fig4_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

``` r
AveConnectivity = geneInfo0 %>%
  group_by(moduleColor) %>%
  summarise(Mean_kOut = mean(kOut))
ind_sort = sort(AveConnectivity$Mean_kOut, index.return = TRUE, decreasing = TRUE)
Module_col = AveConnectivity$moduleColor[ind_sort$ix]; Module_col
```

    ##  [1] "blue"         "red"          "brown"        "cyan"         "turquoise"   
    ##  [6] "yellow"       "magenta"      "black"        "greenyellow"  "green"       
    ## [11] "grey60"       "midnightblue" "salmon"       "grey"

``` r
AveConnectivity$moduleColor <- factor(AveConnectivity$moduleColor, levels = Module_col)
Mean_kOut_p = ggplot(AveConnectivity, aes(x = moduleColor, y = Mean_kOut, col = moduleColor)) +
  geom_hline(yintercept = 0, color = "darkgrey", size = 0.2) +
  geom_point() + 
  scale_y_continuous(name = "Mean connectivity") +
  scale_color_manual(values = Module_col) +
  ggtitle("Intermodular") +
  theme_m ; Mean_kOut_p

p = arrangeGrob(grobs = list(Mean_ktotal_p, Mean_kDiff_p, Mean_kWithin_p, Mean_kOut_p), ncol = 4)
grid.draw(p)
```

![](Fig4_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->

### Fig. 4B: prep inputs needed for circos

``` r
library(tidyverse) # v.2.0.0

# load node file
path =  paste0(getwd(), "/Data/WGCNA/Modules.tsv")
geneInfo0 = read_tsv(path)
```

    ## Rows: 12360 Columns: 30
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (2): GeneSymbol, moduleColor
    ## dbl (28): kTotal, kWithin, kOut, kDiff, GS.T_pgml, GS.Body.Weight..g., GS.Br...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
node = geneInfo0 %>% dplyr::select(GeneSymbol, moduleColor, kWithin)

# load edge file
path = paste0(getwd(), "/Data/WGCNA/edge.txt") # 335845
edge = read_tsv(path)[, 1:3] 
```

    ## Rows: 335845 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): fromNode, toNode, direction
    ## dbl (1): weight
    ## lgl (2): fromAltName, toAltName
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
## make karyotype
node_m = node %>% arrange(moduleColor, desc(kWithin))
sort_col = sort(table(node_m$moduleColor), decreasing = TRUE)

chr = paste0("chr - ", names(sort_col))
id = 1:length(chr)
zero = rep(0, length(chr))
cord = sort_col

karyo = data.frame(chr,  zero, cord, id)
ind = sapply(c("chr", "Var1", "zero", "Freq", "id"), function(s){
  grep(s, names(karyo))
})
karyo = karyo[, ind]
# write_tsv(karyo, file = paste0(getwd(), "/Data/circos/module_karyotype.txt"), col_names = F) 


## make data files
cord_end = lapply(names(sort_col), function(s){
  indd = which(node_m$moduleColor == s);
  df = arrange(node_m[indd, ], desc(kWithin));
  df$end = 1:length(indd)
  return(df)
})
node_m2 = cord_end %>% bind_rows()
node_m2$start = node_m2$end-1

# node
names(node_m2)
```

    ## [1] "GeneSymbol"  "moduleColor" "kWithin"     "end"         "start"

``` r
ind = sapply(c("moduleColor", "start", "end", "GeneSymbol"), function(s){
  grep(s, names(node_m2))
})
cord_out = node_m2[, ind]
# write_tsv(cord_out, file = paste0(getwd(), "/Data/circos/cordinates.txt"), col_names = F) 

## kWithin
ind = sapply(c("moduleColor", "start", "end", "kWithin"), function(s){
  grep(s, names(node_m2))
})
kWithin_out = node_m2[, ind]
# write_tsv(kWithin_out, file = paste0(getwd(), "/Data/circos/kWithin.txt"), col_names = F)

## edge
ind = sapply(c("moduleColor", "start", "end"), function(s){
  grep(s, names(node_m2))
})
line_out = node_m2[, ind]
line_out$y = 1
# write_tsv(line_out, file = paste0(getwd(), "/Data/circos/line_out.txt"), col_names = F)

## make text label
target = c("SP8", "ESR1", "ESR2", "AR", "GATA2", "GATA3", "GATA4", "MYBL1", "SIN3A", 
           "BDNF", "VEGFA", "VEGFC", "KDR", "FLT4", "FLT1")

sapply(target, function(s){ind = which(node$GeneSymbol == s); node$GeneSymbol[unique(ind)]})
```

    ##     SP8    ESR1    ESR2      AR   GATA2   GATA3   GATA4   MYBL1   SIN3A    BDNF 
    ##   "SP8"  "ESR1"  "ESR2"    "AR" "GATA2" "GATA3" "GATA4" "MYBL1" "SIN3A"  "BDNF" 
    ##   VEGFA   VEGFC     KDR    FLT4    FLT1 
    ## "VEGFA" "VEGFC"   "KDR"  "FLT4"  "FLT1"

``` r
sapply(target, function(s){ind = grep(s, node$GeneSymbol); node$GeneSymbol[unique(ind)]})
```

    ## $SP8
    ## [1] "CASP8"    "CASP8AP2" "CRSP8P"   "DUSP8"    "SP8"      "USP8"    
    ## 
    ## $ESR1
    ## [1] "ESR1"
    ## 
    ## $ESR2
    ## [1] "ESR2"
    ## 
    ## $AR
    ##   [1] "AAR2"      "AARS"      "AARS2"     "AARSD1"    "ADAR"      "ADARB1"   
    ##   [7] "ADARB2"    "AP1AR"     "AR"        "ARAP1"     "ARAP2"     "ARC"      
    ##  [13] "ARCN1"     "AREG"      "AREL1"     "ARF1"      "ARF4"      "ARF6"     
    ##  [19] "ARFGAP1"   "ARFGAP2"   "ARFGAP3"   "ARFGEF1"   "ARFGEF2"   "ARFGEF3"  
    ##  [25] "ARFIP1"    "ARFRP1"    "ARG2"      "ARGLU1"    "ARHGAP1"   "ARHGAP10" 
    ##  [31] "ARHGAP11A" "ARHGAP11B" "ARHGAP12"  "ARHGAP15"  "ARHGAP17"  "ARHGAP18" 
    ##  [37] "ARHGAP19"  "ARHGAP20"  "ARHGAP21"  "ARHGAP22"  "ARHGAP24"  "ARHGAP25" 
    ##  [43] "ARHGAP26"  "ARHGAP28"  "ARHGAP29"  "ARHGAP31"  "ARHGAP32"  "ARHGAP36" 
    ##  [49] "ARHGAP39"  "ARHGAP4"   "ARHGAP40"  "ARHGAP42"  "ARHGAP5"   "ARHGAP6"  
    ##  [55] "ARHGDIB"   "ARHGEF10"  "ARHGEF12"  "ARHGEF16"  "ARHGEF17"  "ARHGEF18" 
    ##  [61] "ARHGEF26"  "ARHGEF28"  "ARHGEF3"   "ARHGEF33"  "ARHGEF37"  "ARHGEF38" 
    ##  [67] "ARHGEF39"  "ARHGEF4"   "ARHGEF5"   "ARHGEF6"   "ARHGEF7"   "ARHGEF9"  
    ##  [73] "ARID1A"    "ARID1B"    "ARID2"     "ARID3B"    "ARID3C"    "ARID4A"   
    ##  [79] "ARID4B"    "ARID5B"    "ARIH1"     "ARIH2"     "ARL1"      "ARL11"    
    ##  [85] "ARL13B"    "ARL14"     "ARL14EP"   "ARL15"     "ARL16"     "ARL2BP"   
    ##  [91] "ARL3"      "ARL4A"     "ARL4C"     "ARL5A"     "ARL5B"     "ARL6"     
    ##  [97] "ARL6IP1"   "ARL6IP4"   "ARL6IP5"   "ARL8A"     "ARL8B"     "ARL9"     
    ## [103] "ARMC1"     "ARMC10"    "ARMC2"     "ARMC3"     "ARMC4"     "ARMC4P1"  
    ## [109] "ARMC6"     "ARMC7"     "ARMC8"     "ARMC9"     "ARMT1"     "ARNT"     
    ## [115] "ARNT2"     "ARNTL2"    "ARPC1A"    "ARPC2"     "ARPC3"     "ARPC4"    
    ## [121] "ARPC5"     "ARPC5L"    "ARPP19"    "ARPP21"    "ARR3"      "ARRB1"    
    ## [127] "ARRDC1"    "ARRDC2"    "ARRDC3"    "ARRDC4"    "ARSB"      "ARSD"     
    ## [133] "ARSG"      "ARSI"      "ARSJ"      "ARSK"      "ART4"      "ARTN"     
    ## [139] "ARV1"      "ARVCF"     "ARX"       "BARD1"     "BARHL1"    "BARHL2"   
    ## [145] "BARX1"     "BARX2"     "BCAR1"     "BCAR3"     "BFAR"      "C3AR1"    
    ## [151] "CARD10"    "CARD11"    "CARD19"    "CARD9"     "CARF"      "CARHSP1"  
    ## [157] "CARKD"     "CARNMT1"   "CARS"      "CARS2"     "CARTPT"    "CBARP"    
    ## [163] "CCAR1"     "CCKAR"     "CFLAR"     "DARS"      "EARS2"     "EDAR"     
    ## [169] "EDARADD"   "ERMARD"    "FAR1"      "FAR2"      "FARP1"     "FARP2"    
    ## [175] "FARS2"     "FARSB"     "FFAR4"     "GABARAPL1" "GABARAPL2" "GAR1"     
    ## [181] "GAREM1"    "GARNL3"    "GARS"      "GART"      "HARBI1"    "HARS"     
    ## [187] "HARS2"     "IARS"      "IARS2"     "IFNAR1"    "IFNAR2"    "JARID2"   
    ## [193] "KARS"      "LARGE"     "LARP1"     "LARP1B"    "LARP4B"    "LARP6"    
    ## [199] "LARP7"     "LARS"      "LARS2"     "LPAR1"     "LPAR3"     "LPAR4"    
    ## [205] "LPAR5"     "LPAR6"     "LYAR"      "MARC2"     "MARCH1"    "MARCH11"  
    ## [211] "MARCH2"    "MARCH3"    "MARCH4"    "MARCH5"    "MARCH6"    "MARCH7"   
    ## [217] "MARCH8"    "MARCO"     "MARK1"     "MARK3"     "MARS2"     "MARVELD1" 
    ## [223] "MARVELD2"  "MARVELD3"  "NARF"      "NARFL"     "NARG2"     "NARS2"    
    ## [229] "OARD1"     "PARD3"     "PARD3B"    "PARD6A"    "PARD6B"    "PARD6G"   
    ## [235] "PARG"      "PARK2"     "PARK7"     "PARL"      "PARM1"     "PARN"     
    ## [241] "PARP1"     "PARP11"    "PARP12"    "PARP14"    "PARP15"    "PARP16"   
    ## [247] "PARP3"     "PARP6"     "PARP8"     "PARP9"     "PARPBP"    "PARS2"    
    ## [253] "PARVB"     "PARVG"     "PPARA"     "PPARD"     "PPARG"     "PPARGC1A" 
    ## [259] "PPARGC1B"  "PRKAR1A"   "PRKAR1B"   "PRKAR2A"   "PRKAR2B"   "PTAR1"    
    ## [265] "QARS"      "RARB"      "RARRES1"   "RARRES2"   "RARS"      "RARS2"    
    ## [271] "SAR1A"     "SAR1B"     "SARAF"     "SARDH"     "SARM1"     "SARS"     
    ## [277] "SART3"     "SCARB1"    "SCARB2"    "SCARF1"    "SCARF2"    "SIGMAR1"  
    ## [283] "SMARCA1"   "SMARCA2"   "SMARCA5"   "SMARCAD1"  "SMARCAL1"  "SMARCB1"  
    ## [289] "SMARCC1"   "SMARCD2"   "SMARCD3"   "SMARCE1"   "SPARC"     "SPARCL1"  
    ## [295] "STAR"      "STARD10"   "STARD13"   "STARD3"    "STARD3NL"  "STARD4"   
    ## [301] "STARD5"    "STARD8"    "STARD9"    "TAAR1"     "TARBP1"    "TARDBP"   
    ## [307] "TARS"      "TARS2"     "TARSL2"    "TIGAR"     "TSNARE1"   "TSPEAR"   
    ## [313] "WARS"      "WARS2"     "YARS"      "YARS2"     "ZAR1"      "ZAR1L"    
    ## 
    ## $GATA2
    ## [1] "GATA2"
    ## 
    ## $GATA3
    ## [1] "GATA3"
    ## 
    ## $GATA4
    ## [1] "GATA4"
    ## 
    ## $MYBL1
    ## [1] "MYBL1"
    ## 
    ## $SIN3A
    ## [1] "SIN3A"
    ## 
    ## $BDNF
    ## [1] "BDNF"
    ## 
    ## $VEGFA
    ## [1] "VEGFA"
    ## 
    ## $VEGFC
    ## [1] "VEGFC"
    ## 
    ## $KDR
    ## [1] "KDR"    "PKDREJ"
    ## 
    ## $FLT4
    ## [1] "FLT4"
    ## 
    ## $FLT1
    ## [1] "FLT1"

``` r
text_lab = subset(cord_out, su = GeneSymbol %in% target)
# write_tsv(text_lab, file = paste0(getwd(), "/Data/circos/text_labels_0726.txt"), col_names = F)


## make link file
link_from = cord_out[match(edge$fromNode, cord_out$GeneSymbol),]
link_to   = cord_out[match(edge$toNode, cord_out$GeneSymbol),]
link = cbind(link_from, link_to)
names(link)
```

    ## [1] "moduleColor" "start"       "end"         "GeneSymbol"  "moduleColor"
    ## [6] "start"       "end"         "GeneSymbol"

``` r
ind = sapply(c("GeneSymbol"), function(s){
  grep(s, names(link))
})
link = link[, -ind]
link = as.data.frame(link)

## binary color SP8 nodes (more color) vs non-SP8node
link$para = paste0("color=", "vlgrey", "_a3") # a1(16%)~a5(83%) for alpha

ind_from_sp8 = which(edge$fromNode == "SP8")
ind_to_sp8 = which(edge$toNode == "SP8")
ind_sp8 = c(ind_from_sp8, ind_to_sp8)

color_blues = c("vvlblue", "vlblue", "blue") # the color defined in circos
color_blues_2 = rep(color_blues, length.out = length(ind_sp8))
link$para[ind_sp8] = paste0("color=", color_blues_2, "_a3")

link = arrange(link, desc(para))
# write_tsv(link, file = paste0(getwd(), "/Data/circos/link_sp8_blues3_vlgrey_0726.txt"), col_names = F)


## for circos.conf
# chromosomes_color
paste0("/", names(sort_col), "/:", names(sort_col), collapse = ";")
```

    ## [1] "/turquoise/:turquoise;/green/:green;/blue/:blue;/black/:black;/brown/:brown;/yellow/:yellow;/red/:red;/magenta/:magenta;/greenyellow/:greenyellow;/salmon/:salmon;/cyan/:cyan;/midnightblue/:midnightblue;/grey60/:grey60;/grey/:grey"

``` r
# link color
paste0("/", names(sort_col), "/:", names(sort_col), "_a4", collapse = ";")
```

    ## [1] "/turquoise/:turquoise_a4;/green/:green_a4;/blue/:blue_a4;/black/:black_a4;/brown/:brown_a4;/yellow/:yellow_a4;/red/:red_a4;/magenta/:magenta_a4;/greenyellow/:greenyellow_a4;/salmon/:salmon_a4;/cyan/:cyan_a4;/midnightblue/:midnightblue_a4;/grey60/:grey60_a4;/grey/:grey_a4"

``` r
# color def
col_def = sapply(names(sort_col), function(s){
  paste0(s, " = ", paste0(col2rgb(s), collapse = ","))
})
cat(col_def, sep = "\n")
```

    ## turquoise = 64,224,208
    ## green = 0,255,0
    ## blue = 0,0,255
    ## black = 0,0,0
    ## brown = 165,42,42
    ## yellow = 255,255,0
    ## red = 255,0,0
    ## magenta = 255,0,255
    ## greenyellow = 173,255,47
    ## salmon = 250,128,114
    ## cyan = 0,255,255
    ## midnightblue = 25,25,112
    ## grey60 = 153,153,153
    ## grey = 190,190,190

``` r
# then run circos-0.69-4  ## http://circos.ca
# circos -conf ./circos.conf 
```

### Fig. 4C

``` r
library(tidyverse) # v.2.0.0
library(HiveR) # 0.3.63

# load string output
path = paste0(getwd(), "/Data/WGCNA/SP8_string_nodes_20180531.tsv") 
node = read_tsv(path)
```

    ## Rows: 12 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (12): Category, Category2, Color, description, gene.symbol, id.passing.f...
    ## dbl (10): SUID, ...2, gene.id, input.gene.id, T14d, T1h, T3d, T3h, T7d, T8h
    ## lgl  (1): selected
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# select modulecolor for axises
unique(node$regulatory.function)
```

    ## [1] "kinase" "actin"  "TF_co"  "TF"

``` r
sort(table(node$regulatory.function))
```

    ## 
    ##  actin  TF_co kinase     TF 
    ##      2      2      3      5

``` r
node$regulatory.function = gsub("TF_co", "TF", node$regulatory.function)
sel_col = c("TF", "kinase", "actin")
node_sub = subset(node, su = regulatory.function %in% sel_col, se = c("regulatory.function", "gene.symbol"))
unique(node_sub$regulatory.function)
```

    ## [1] "kinase" "actin"  "TF"

``` r
# load edge file
path = paste0(getwd(), "/Data/WGCNA/string_interactions_tur_adjwithSP8_1151_ZF.tsv")
edge = read_tsv(path)
```

    ## Rows: 3035 Columns: 15
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): node1, node2, node1_external_id, node2_external_id
    ## dbl (11): node1_string_internal_id, node2_string_internal_id, neighborhood_o...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# subset edges according to node
## only direct neighbors
edge_sub = subset(edge, su = node1 %in% c("SP8") | node2 %in% c("SP8"))
all_edges_sub_node = unique(c(edge_sub$node1, edge_sub$node2))
node_sub = subset(node_sub, su = gene.symbol %in% all_edges_sub_node)

# check how many module colors are left
all_edges_sub_node = unique(c(edge_sub$node1, edge_sub$node2))
unique(node_sub$regulatory.function[which(node_sub$gene.symbol %in% all_edges_sub_node)])
```

    ## [1] "kinase" "actin"  "TF"

``` r
## manuelly making NODE
# lab
names(node_sub) = gsub("gene.symbol", "lab", names(node_sub))
# col
df = data.frame(regulatory.function = sel_col, color = c("#FF9900", "#00cc33", "#00ccff"))
node_sub$color = df$color[match(node_sub$regulatory.function, df$regulatory.function)]
# id
node_sub$id = 1:nrow(node_sub)
# size
node_sub$size  = 3
range(node_sub$size)
```

    ## [1] 3 3

``` r
# axis
temp = data.frame(mod_col = sel_col, axis = 1:length(sel_col)); temp
```

    ##   mod_col axis
    ## 1      TF    1
    ## 2  kinase    2
    ## 3   actin    3

``` r
node_sub$axis = temp$axis[match(node_sub$regulatory.function, temp$mod)]

# load WGCNA edge file
path = paste0(getwd(), "/Data/WGCNA/edge.txt") # 335845
edge_W = read_tsv(path)[, 1:3] 
```

    ## Rows: 335845 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): fromNode, toNode, direction
    ## dbl (1): weight
    ## lgl (2): fromAltName, toAltName
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# radius
connections = lapply(node_sub$lab, function(s){
  ind_c = c(which(edge_W$fromNode == s), which(edge_W$toNode == s)); # using WGCNA edge here
  data.frame("lab" = s, "count" = length(ind_c))
}) %>% bind_rows()
node_sub = merge(node_sub, connections, by = "lab")
radius = lapply(sel_col, function(s){
  ind = which(node_sub$regulatory.function == s);
  # df = arrange(node_sub[ind, ], desc(count))
  df = arrange(node_sub[ind, ], desc(lab))
  df$radius = seq(20, 200, by = 25)[1:length(ind)]
  df
}) %>% bind_rows()

# node_out
ind = sapply(c("id", "lab", "axis", "radius", "size", "color"), function(s){
  grep(s, names(radius))
})
node_out = radius[, ind]

## manuelly making EDGE
edge_sub$id1 = node_out$id[match(edge_sub$node1, node_out$lab)]
edge_sub$id2 = node_out$id[match(edge_sub$node2, node_out$lab)]
ind = sapply(c("id1", "id2", "combined_score"), function(s){
  grep(s, names(edge_sub))
})
edge_out = edge_sub[, ind] 
names(edge_out) = gsub("combined_score", "weight", names(edge_out))
edge_out$weight = scales::rescale(edge_out$weight, c(0.5, 2))
edge_out$color = "lightgray" 
edge_out = arrange(edge_out, desc(color)) %>% as.data.frame() # first color will be at the bottom
str(edge_out)
```

    ## 'data.frame':    11 obs. of  4 variables:
    ##  $ id1   : int  2 12 3 8 10 9 11 5 7 11 ...
    ##  $ id2   : int  11 11 11 11 11 11 6 11 11 1 ...
    ##  $ weight: num  2 2 1.39 1.35 1.3 ...
    ##  $ color : chr  "lightgray" "lightgray" "lightgray" "lightgray" ...

``` r
# manuelly creating HPD
H1 = list()
H1$nodes = node_out
H1$edges = edge_out
H1$type = "2D"
H1$desc = "cytoscape"
H1$axis.cols = rep('darkgray', length(sel_col)) # make invisible
H1$axLabs = sel_col
class(H1) = "HivePlotData"
chkHPD(H1, confirm = TRUE)
```

    ## You must be awesome: This hive plot data looks dandy!

    ## [1] FALSE

``` r
# sumHPD(H1)

anNodes_path = paste0(getwd(), "/Data/WGCNA/HiveR_anNodes_temp.tsv")
df = data.frame("node.lab" = node_out$lab)
df$node.text = node_out$lab
df$angle = 0
df$radius = 0
df$offset = 0
df$hjust = 0.5
df$vjust = 0.5
# write_csv(df, file = anNodes_path)

plotHive(H1, axLabs = sel_col, ch = 5, bkgnd = "transparent",
         axLab.pos = 20, 
         anNodes = anNodes_path,
         anNode.gpar = gpar(col = "black", fontsize = 6, lwd = 0), 
         axLab.gpar = gpar(col = "black", fontsize = 6) ) 
```

![](Fig4_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->