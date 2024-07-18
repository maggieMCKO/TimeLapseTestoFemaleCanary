<!---
This README uses Markdown syntax
-->

# Reference Information

## Provenance for this README

This README.md was generated on 17-July-2024 by Meng-Ching Ko.

This depository included all datasets needed for reproducing figures, supplementary files, and supplementary figures in "Extensive, transient, and long-lasting gene regulation in a song-controlling brain area during testosterone-induced song development in adult female canaries."

The scripts for reproducing figures, supplementary files, and supplementary figures were deposited on GitHub ([https://github.com/maggieMCKO/TimeLapseTestoFemaleCanary](https://github.com/maggieMCKO/TimeLapseTestoFemaleCanary)).

* AUTHOR INFORMATION
  * Ko, Meng-Ching, Department of Behavioural Neurobiology, Max Planck Institute for Biological Intelligence
    * orcid: [https://orcid.org/0000-0002-2234-9380](https://orcid.org/0000-0002-2234-9380)
    * email: [mengching.ko@bi.mpg.de](mailto:mengching.ko@bi.mpg.de)
    * alternate email: [mengchingko@gmail.com](mailto:mengchingko@gmail.com)
  * Frankl-Vilches, Carolina, Department of Behavioural Neurobiology, Max Planck Institute for Biological Intelligence
  * Bakker, Antje, Department of Behavioural Neurobiology, Max Planck Institute for Biological Intelligence
  * Sohnius-Wilhelmi, Nina, Department of Behavioural Neurobiology, Max Planck Institute for Biological Intelligence
  * Alcami, Pepe, Department of Behavioural Neurobiology, Max Planck Institute for Biological Intelligence; Division of Neurobiology, Faculty of Biology, Ludwig-Maximilians-University Munich
  * Gahr, Manfred, Department of Behavioural Neurobiology, Max Planck Institute for Biological Intelligence

## FILE OVERVIEW

### Table of Contents

# Dataset Overview

This repository contains the following files that are critical for reproducing the data analysis and results. Each file is briefly described to help you understand its purpose and content.

## Table of Contents

### Primary Data Files

* **DiffExpression.tsv**
  * The output of the differential analysis is stored in this file.
* **HVC_volume.tsv**
  * This file includes data on HVC volume and normalized HVC volume.
* **Physiological_measurements.tsv**
  * This file contains measurements of brain weight, body weight, and oviduct weight.
* **PlasmaAndrogenLv.tsv**
  * This file contains data on plasma androgen levels.
* **RNAscope_NormalizedStainedArea.tsv**
  * This file includes the quantification of RNAscope staining.
* **Expression_est_perBird.tsv**
	* output of RMA (Robust Multichip Analysis).

### Song analysis-related files

* **DayLv_data.tsv**
  * This file contains daily singing rates.
* **SongLv_data.tsv**
  * Song parameters such as song length, syllable repetition rate, and the number of syllables per song are stored in this file.

### File/Folder Details

#### Details for: DiffExpression.tsv

* Description: This tab-delimited file contains the filtered output of the differential expression analysis performed using the R package limma. It includes only genes with human orthologous annotation meeting the significance criteria (adjusted P-value, fold change, and power thresholds) and showing consistent direction of regulation, as described in the 'Analysis of differentially expressed genes' section of the associated manuscript. The file contains normalized expression values, log fold changes, adjusted p-values, and power values for these differentially expressed genes across different time points.
* Format(s): .tsv
* Variables:
  * GeneSymbol: Differentially expressed genes
  * logFC: Log2 fold change of the gene
  * AveExpr: Average expression of the gene across all samples
  * t: t-statistic
  * P.Value: Raw p-value
  * adj.P.Val: Adjusted p-value (Benjamini-Hochberg method)
  * B: B-statistic (log-odds that the gene is differentially expressed)
  * contrast: The group in which the gene was identified as differentially expressed (relative to the control group)
  * direction: Direction of regulation (-1 for down-regulated, 1 for up-regulated)
  * logFC_ave: The average log2 fold change across all probes for each differentially expressed gene
  * power: Power value for the gene (≥ 0.8)
  * TF: Indicates whether the gene is a transcription factor (1) or not (0)
* Missing data codes: Cells with missing data are coded as NA ("NA")

#### Details for: HVC\_volume.tsv

* Description: This tab-delimited file provides the HVC volume data for the testosterone-treated group. Both raw and body weight-normalized HVC volumes are included.
* Format(s): .tsv
* Variables:
  * Group: Testosterone-treated group
  * ProcessNum: Bird (process) ID
  * HVC Volume (mm3): HVC volume
  * HVC_nor: HVC volume normalized by body weight

#### Details for: Physiological\_measurements.tsv

* Description: This tab-delimited file contains measurements of brain weight, body weight, and oviduct weight.
* Format(s): .tsv
* Variables:
  * Group: Testosterone-treated group
  * ProcessNum: Bird (process) ID
  * Brain (mg): Brain weight
  * Body Weight (g): Body weight
  * Oviduct (mg): Oviduct weight
* Missing data codes: Cells with missing data are coded as NA ("NA")

#### Details for: PlasmaAndrogenLv.tsv

* Description: This tab-delimited file contains data on plasma androgen levels.
* Format(s): .tsv
* Variables:
  * Group: Testosterone-treated group
  * PrePost: Blood sampled before or after testosterone implantation
  * ProcessNum: Bird (process) ID
  * ngml: Plasma androgen levels
  * date: Date of blood sampling (yyyy/mm/dd)

#### Details for: RNAscope\_NormalizedStainedArea.tsv

* Description: This tab-delimited file contains the quantification of RNAscope staining. Please see the methods section "RNAScope® in situ hybridization assay" in the manuscript for further details.
* Format(s): .tsv
* Variables:
  * Probe: Gene examined using RNAscope
  * Group: Testosterone-treated group
  * Individual: Bird ID
  * Slide: Slide ID
  * sum_slide: The sum of the area of accumulated chromogenic particles in the stained HVC sections
  * Area: The HVC volume of the section
  * n_particle: The number of accumulated chromogenic particles in the stained HVC sections
  * Density: The density of accumulated chromogenic particles within the HVC (n_particle/Area)
  * percent_area_slide: The area of accumulated chromogenic particles normalized by HVC volume
* Missing data codes: Cells with missing data are coded as NA ("NA")

#### Song analysis-related files

#### Details for: DayLv\_data.tsv

* Description: This tab-delimited file contains daily singing rates.
* Format(s): .tsv
* Variables:
  * Ind: Bird ID
  * Individual: Bird ID as integer
  * Group: Testosterone-treated group
  * n_Group: Testosterone-treated group as integer
  * n_Day: days after testosterone implantation
  * Daily_SongRate: daily song rate

#### Details for: SongLv\_data.tsv

* Description: This tab-delimited file contains song parameters such as song length, syllable repetition rate, the number of syllables per song, and the slope coefficient (α).
* Format(s): .tsv
* Variables:
  * Ind: Bird ID
  * Individual: Bird ID as integer
  * Group: Testosterone-treated group
  * n_Group: Testosterone-treated group as integer
  * FisrtDate: Testosterone implantation date
  * n_Day: days after testosterone implantation
  * StartTime: The recording date
  * Daily_SongRate: daily song rate
  * SongLen: Song length
  * SylperSong: number of syllables in a song
  * RepetitionRate: syllable repetition rate
  * slope: the slope coefficient (α)

#### Details for: Expression\_est\_perBird.tsv

* Description: This tab-delimited file contains the output of the Robust Multichip Analysis (RMA). Each column from the second to the 12361st represents the average expression of a specific gene.
* Format(s): .tsv
* Variables:
  * Sample: Sample (Bird) ID
  * Column 2 to 12360: Average expression of each gene (each column corresponds to a specific gene)

