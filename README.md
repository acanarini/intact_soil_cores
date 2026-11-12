# Deuterium isotope probing of microbial growth in intact soil cores

This repository contains the **analysis-ready data and R scripts** associated with the accepted manuscript *Deuterium isotope probing of microbial growth in intact soil cores* by Alberto Canarini, Gaëlle Marmasse, Lucia Fuchslueger, Andreas Richter, and Hannes Schmidt.

The study compares microbial growth and related processes in **intact soil cores versus sieved soil**, using deuterium (²H) labelling of soil water by water-vapour equilibration and phospholipid fatty acid (PLFA) analysis. Samples were obtained from three ecosystems: **forest, cropland, and meadow**. The study includes five paired intact-core/sieved-soil replicates per ecosystem (30 observations in the main analysis tables).

## Repository contents

### R scripts

| File | Description |
| --- | --- |
| `Code.R` | Main analysis script: PLFA growth and abundance ordinations; PERMANOVAs; distances to group centroids and dispersion analyses; microbial-group growth, respiration, and carbon use efficiency (CUE) plots; mixed-effects analyses and associated diagnostics. Generates figure files with names such as `Fig. 2.png`, `Fig. S2.png`, `Fig. S3.png`, and `Fig. S4.png`. |
| `2Hincorporation.R` | Fits water-vapour ²H equilibration curves for intact and sieved soil in all three ecosystems, plots observed and predicted enrichment, and calculates fitted-curve integrals. Saves individual ecosystem plots and a combined SVG plot. |

### Input data

| File | Contents | Observations |
| --- | --- | ---: |
| `MultivariateMS_GROWTH.csv` | Mass-specific growth values for individual PLFA biomarkers, plus ecosystem, treatment and core identifiers. | 30 |
| `MultivariateABUNDANCE.csv` | PLFA abundance values and corresponding sample identifiers. **Semicolon-delimited**. | 30 |
| `Taxon-specific-bmcgrowth.csv` | Aggregated growth measurements for microbial biomarker groups (e.g., Gram-positive and Gram-negative bacteria, actinobacteria, fungi and AMF). | 30 |
| `Resp_Growth_CUE.csv` | Community-level respiration, mass-specific growth and carbon use efficiency (CUE). | 30 |
| `2Hincorporation.csv` | Measured ²H enrichment (atom %) of labelled water over equilibration time (hours) for intact and sieved treatments. | 66 |

For the four 30-row analysis tables, `Ecosystem` identifies the soil type; `Treatment` distinguishes intact cores (`IC`) from `Sieved`; `Core` identifies the paired samples within each ecosystem; and `Combined` combines ecosystem and treatment labels. The equilibration dataset uses `type` (Forest, Meadow, Field) and `Treatment` (`IC`, `S`); **Field** corresponds to **Cropland** in the main analysis tables.

These CSVs are inputs to the supplied analyses, not raw instrument export files.

## Running the analyses

1. Download or clone this repository, then open it in an interactive **R / RStudio** session.
2. Install the necessary packages if they are not already available:

   ```r
   install.packages(c(
     "ggplot2", "ggpubr", "reshape2", "forcats", "tidyverse",
     "dplyr", "ggthemes", "RColorBrewer", "vegan", "permute",
     "FactoMineR", "corrplot", "factoextra", "patchwork",
     "rstatix", "nlme", "multcomp", "emmeans", "lmtest",
     "qqplotr", "pracma"
   ))
   ```

3. Run either script from the repository directory, for example:

   ```r
   source("Code.R")
   # Or run the equilibration analysis separately:
   source("2Hincorporation.R")
   ```

**Important:** Both scripts use `file.choose()` to select input files interactively. When running `Code.R`, select the following CSV files in this order as each dialog appears:

1. `MultivariateMS_GROWTH.csv`
2. `Taxon-specific-bmcgrowth.csv`
3. `MultivariateABUNDANCE.csv`
4. `Resp_Growth_CUE.csv`
5. `Taxon-specific-bmcgrowth.csv` (again, for supplementary plots)

When running `2Hincorporation.R`, select `2Hincorporation.csv`.

Figure files are saved to the current R working directory; the scripts also display plots and print statistical results. The `file.choose()` interface requires an interactive environment, and the scripts are not configured as an automated or package-version-locked analysis pipeline. Figure filenames in the scripts are working filenames; consult the accepted manuscript for final figure numbering.

## Statistical distance metrics

The supplied `Code.R` uses **Euclidean distances** for the PLFA analyses: `adonis2(..., method = "euclidean")` for PERMANOVA and `vegdist(..., method = "euclid")` followed by `betadisper(..., type = "centroid")` for distances to group centroids. See the manuscript for the interpretation of the multivariate and dispersion results.

## Reference and citation

**Associated manuscript:** Canarini, A., Marmasse, G., Fuchslueger, L., Richter, A., and Schmidt, H. *Deuterium isotope probing of microbial growth in intact soil cores*. Accepted manuscript; final journal citation and DOI to be added when available.

For experimental design, isotope-labelling procedures, PLFA extraction, and interpretation of the results, please refer to the associated manuscript and its supplementary materials. A version-specific repository citation/DOI can be added after the code and data are archived.
