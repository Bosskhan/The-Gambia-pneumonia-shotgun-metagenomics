# README

## PERCH MAG and Microbiome Analysis

This repository contains the R code used to generate the main and supplementary figures and perform the statistical analyses in the revised manuscript.
The supporting input files can be requested from Dam kHAN: dam.khan@lshtm.ac.uk


## 1. Analysis overview

The analysis workflow consists of the following main components:

1. Construction of a `phyloseq` object from OTU, taxonomy, and sample metadata tables.
2. Filtering of taxa
3. Generation of the main figures:

   * Figure 1: Sample availability
   * Figure 2: MAG antimicrobial resistance profiles
   * Figure 3: Bacterial composition by specimen type
   * Figure 4: Shannon diversity
   * Figure 5: Jaccard similarity between paired blood and respiratory specimens
   * Figure 6: NPS bacterial composition by disease severity
   
4. Generation of supplementary figures:
   * Figure S1: Individual-level bacterial composition
   * Figure S2: Dominant species by specimen type
   * Figure S3: NPS Shannon diversity by disease severity
   * Figure S4: NPS beta diversity and PCoA
   
6. Statistical analyses including:

   * Linear mixed-effects modelling
   * Estimated marginal means and Tukey-adjusted pairwise comparisons
   * Wilcoxon rank-sum tests
   * PERMANOVA
   * MaAsLin3 differential-abundance analysis
   * Permutation testing of paired versus randomly paired Jaccard similarities

---

## 2. Input files

The following files are required by the analysis script.


| File              | Purpose                                                       |
| ----------------- | ------------------------------------------------------------- |
| `otu_table.csv`   | OTU/ASV abundance table                                       |
| `tax_table.csv`   | Taxonomic assignments corresponding to the OTU/ASV table      |
| `sample_data.csv` | Sample-level metadata used to construct the `phyloseq` object |

These three files are used to construct the initial phyloseq object.

### Figure 1

| File                      | Purpose                                                   |
| ------------------------- | --------------------------------------------------------- |
| `Sample_availability.csv` | Sample availability information used to generate Figure 1 |

### Figure 2

| File                       | Purpose                                                          |
| -------------------------- | ---------------------------------------------------------------- |
| `amr_prevalence_table.csv` | Antimicrobial resistance gene prevalence among high-quality MAGs |
| `amr_data.csv`             | Per-MAG antimicrobial resistance gene presence/absence data      |

### Figure 5

| File                              | Purpose                                                                                 |
| --------------------------------- | --------------------------------------------------------------------------------------- |
| `combined_jaccard_index_data.csv` | Selected pairwise Jaccard similarity values used for the descriptive Figure 5 summaries |

The complete 71 ?? 71 Jaccard similarity matrix is generated directly from the microbiome data during the analysis and is not required as a separate input file.

### Figure S2

| File                   | Purpose                                                 |
| ---------------------- | ------------------------------------------------------- |
| `dominant_species.csv` | Dominant species information used to generate Figure S2 |

---

## 3. Output files

Running the analysis produces the following figure files:

### Main figures


Figure 1_sample availability.png
figure_2_MAGs_AMR.png
Figure 3_compositional barplot.png
figure_4_shannon_diversity.png
figure_5_jaccard.png
figure_6_NPS_severity.png


### Supplementary figures

figure_s1_individual_barplots.png
figure_s2_dominant_species_pies.png
figure_s3_NPS_diversity_by_severity.png
figure_s4_NPS_PCoA.png


The Jaccard permutation analysis also produces a console output table containing:

```text
pair
n_pairs
observed_mean
null_mean
null_sd
p_value
q_value
significant
```

---

## 4. Software requirements

The analysis was performed in R.

Recommended R version:

```text
R >= 4.3
```

The following R packages are used by the analysis:

```r
library(ggplot2)
library(phyloseq)
library(maaslin3)
library(dplyr)
library(readxl)
library(tidyr)
library(writexl)
library(patchwork)
library(microViz)
library(ggtext)
library(ggpubr)
library(stringr)
library(tibble)
library(ggsignif)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggrepel)
```

Additional packages/functions used explicitly through their namespace include:

```r
ggh4x
cowplot
grid
purrr
```

For reproducibility, these packages should also be installed.

---

## 5. Running the analysis

Place the R analysis script and all required input files in the same directory.

The core files should include:

```text
otu_table.csv
tax_table.csv
sample_data.csv
Sample_availability.csv
amr_prevalence_table.csv
amr_data.csv
combined_jaccard_index_data.csv
dominant_species.csv
```

Open the R script in RStudio and set the working directory to the directory containing these files, or use an RStudio Project.

The analysis can then be run sequentially from beginning to end.

---

## 6. Microbiome preprocessing

The three input tables are used to construct a `phyloseq` object.

Taxonomic assignments containing `Firmicutes_` prefixes are standardised to:

```text
Firmicutes
```

Relative abundance is calculated within each sample.

Taxa contributing less than 0.1% relative abundance within an individual sample are set to zero.

The following taxa are subsequently removed because they were considered potential contaminants:

```text
Cutibacterium acnes
Staphylococcus epidermidis
Staphylococcus saprophyticus
```

The resulting object is used for downstream community composition and diversity analyses.

---

## 7. Jaccard similarity and permutation analysis

Species-level Jaccard similarity is calculated after aggregating the abundance data to species level and converting the data to binary presence/absence.

The Jaccard similarity matrix contains pairwise similarities for all 71 specimens.

For each blood-respiratory specimen comparison, participants with both specimen types are identified.

The observed mean Jaccard similarity is calculated across the true participant-matched pairs.

A permutation test is then performed by keeping one specimen type fixed and randomly shuffling the other specimen type among the same participants. This breaks the original participant-level pairing while retaining the observed specimens and number of pairs.

The mean Jaccard similarity is recalculated for each permutation.

The procedure is repeated 9,999 times for each comparison.

The one-sided permutation P value is calculated as:

```text
number of permuted mean similarities >= observed mean similarity
---------------------------------------------------------------
                         9,999
```

Separate tests are performed for:

```text
NPS vs Blood
IS vs Blood
NP-OP vs Blood
```

The resulting P values are adjusted using the Benjamini???Hochberg procedure.

Statistical significance is defined as:

```text
q < 0.05
```

The permutation tests use:

```r
seed = 123
```

to ensure reproducibility.

---

## 8. Alpha diversity

Species-level Shannon diversity is calculated using `microViz`.

For the main comparison among specimen types, a linear mixed-effects model is fitted using:

```text
Shannon diversity ~ specimen type + (1 | participant ID)
```

Participant ID is included as a random intercept to account for repeated specimens from the same participant.

Pairwise specimen-type comparisons are obtained using estimated marginal means with Tukey adjustment.

---

## 9. Beta diversity

For the NPS severity analysis, Aitchison distance is calculated following compositional transformation.

Principal coordinates analysis (PCoA) is used for visualisation.

PERMANOVA is performed to assess associations between community composition and:


case-control classification
age
vaccination status


The analysis uses 999 permutations and a fixed seed of 123.

---

## 10. Differential abundance analysis

MaAsLin3 is used to assess differences in species abundance between severe and very severe pneumonia cases among NPS specimens.

The model includes:


case-control classification
age


The abundance data are CLR-normalised without an additional transformation.

The MaAsLin3 analysis uses a minimum prevalence threshold of 15%.

---

## 11. Vaccination analyses

For the NPS analyses, vaccination status is collapsed into two categories:

full
incomplete


where:

text
full       = full vaccination
incomplete = none or partial vaccination


Separate Wilcoxon rank-sum tests are used to compare:

* *Streptococcus pneumoniae* relative abundance by PCV vaccination status
* *Haemophilus influenzae/aegyptius* relative abundance by Hib vaccination status

The two resulting P values are adjusted using the Benjamini???Hochberg procedure.
