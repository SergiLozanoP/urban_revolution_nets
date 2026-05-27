# The Network Structure of the Urban Revolution

**Authors:** Giacomo Benati & Sergi Lozano

This repository contains the R code and datasets associated with the manuscript *"The Network Structure of the Urban Revolution"*, currently under consideration for publication.

---

## Overview

This project investigates the political organisation of Mesopotamia during the **Uruk period** (~3800–3000 BCE) by analysing the structure of transportation networks among archaeological sites with Uruk material culture. Networks are reconstructed for three consecutive chronological phases — LC3, LC4, and LC5 — and their structural properties are compared against null models to assess whether observed connectivity patterns deviate from random expectations.

---

## Chronological Phases

| Phase | Period | Date range |
|-------|--------|------------|
| LC3 | Late Chalcolithic 3 | ~3800–3600 BCE |
| LC4 | Late Chalcolithic 4 | ~3600–3300 BCE |
| LC5 | Late Chalcolithic 5 | ~3300–3000 BCE |

---

## Repository Structure

```
/
├── Uruk_nets_submission.R   # Main analysis script
└── Data/
    ├── LC_lat-long_v4.csv   # Geographic coordinates of all sites
    ├── LC3_all_v6.csv        # Edge list: all connections (LC3)
    ├── LC3_water_v6.csv      # Edge list: river connections only (LC3)
    ├── LC4_all_v6.csv        # Edge list: all connections (LC4)
    ├── LC4_water_v6.csv      # Edge list: river connections only (LC4)
    ├── LC5_all_v6.csv        # Edge list: all connections (LC5)
    ├── LC5_water_v6.csv      # Edge list: river connections only (LC5)
    ├── discarded.csv         # Sites excluded due to low chronological precision
    └── Non-Uruk.csv          # Sites excluded due to absence of Uruk assemblage
```

---

## Data Description

### Nodes
Each node represents an **archaeological site** where Uruk material culture has been documented in the existing literature. Geographic coordinates (latitude and longitude) are provided in `LC_lat-long_v4.csv`.

### Edges
Two types of transportation connections are modelled:

- **River connections** (`*_water_v6.csv`): links between sites along waterways, forming the backbone of the transportation network.
- **All connections** (`*_all_v6.csv`): river connections plus overland routes.

### Excluded sites
Two sets of sites were excluded from the main analysis and are provided for transparency:

- `discarded.csv`: sites with insufficient chronological precision to be assigned to a specific LC phase.
- `Non-Uruk.csv`: sites where no Uruk evidence has been found to date.

---

## Analysis Pipeline

The script `Uruk_nets_submission.R` performs the following steps:

### 1. Network construction
- Loads edge lists and builds `igraph` objects for each phase and connection type.
- Adds geographic coordinates as vertex attributes.
- Exports networks as `.GML` files for visualisation in external tools (e.g. Gephi).

### 2. Spatial validation
- Compares the geographic distribution of included vs. excluded sites using the `check_nodes()` function.
- Produces overlaid maps to assess whether excluded sites occupy spatial regions not covered by included ones.

### 3. Network robustness tests
- Implements the robustness testing procedure proposed by Peeples & Brughmans for archaeological networks.
- Iteratively removes nodes or edges and measures the impact on degree and betweenness centrality rankings (Spearman correlation).
- Results are presented as boxplots across removal percentages.

### 4. Basic network characterisation
- Computes key structural metrics: number of nodes, density, average clustering coefficient, average path length, degree centralisation, and betweenness centralisation.

### 5. Computational experiments: null models
Two null models are implemented to benchmark empirical network properties:

- **Random model** (`random_nm()`): rewires non-river edges randomly, subject to a maximum geographic distance constraint.
- **Nearest-neighbour model** (`nearneigh_nm()`): rewires non-river edges by connecting nodes to their geographically closest available neighbour.

Both models preserve the river-connection backbone and are run for 100 replications. Results are compared against empirical values for clustering coefficient, average path length, degree centralisation, and betweenness centralisation.

The **small-worldness index** (Humphries & Gurney, 2008) is also computed from the random model output.

---

## Dependencies

The following R packages are required:

```r
install.packages(c("readr", "dplyr", "igraph", "ggplot2", "geosphere", "sf"))
```

---

## Usage

1. Clone or download this repository.
2. Set your working directory at the top of `Uruk_nets_submission.R`:
   ```r
   setwd("path/to/your/local/folder")
   ```
3. Make sure all CSV files are placed inside a `Data/` subfolder.
4. Run the script section by section, following the comments.

---

## Citation

If you use this code or data, please cite the associated manuscript:

> Benati, G. & Lozano, S. (*under review*). The Network Structure of the Urban Revolution. *[Journal name to be added upon acceptance]*
>
> DOI: *[to be added upon acceptance]*

---

## License

- **Code** (`Uruk_nets_submission.R`): [MIT License](https://opensource.org/licenses/MIT) — free to use, modify and redistribute with attribution.
- **Data** (CSV files): [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — free to use and redistribute, provided appropriate credit is given to the authors.
