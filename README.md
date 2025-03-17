# Spatial-Autoregressive-Geometric-Mean-Model

This repository contains the R code to fit the **spatial autoregressive geometric mean model** to spatial count data, such as epidemiological data. The methodology is described in the paper: **_Spatial Autoregressive Modelling of Epidemiological Data: Geometric Mean Model Proposal._**, published in SORT (https://www.idescat.cat/sort/) in 2025.

## Script structure 
1. The dataset is simulated using the geometry of the Columbus shapefile from the `spData` library  

2. A mobility matrix is simulated, representing the mean proportion of time that people from region *i* spend in region *j* during a given time period  

3. Population counts are simulated

4. A spatially autocorrelated variable representing infectious disease counts is generated using a Gibbs sampler

4. The geometric mean spatial conditional model is fitted to the dataset. Different spatial weight matrices are considered:  
     - First-order contiguity 
     - Inverse distance 
     - Simulated mobility matrix
The model is estimated using the Integrated Nested Laplace Approximation (INLA) method via the `INLA` R package.  

## Dependencies  
This code has been tested with the following R version and package versions:  

- R version: 4.4.0  
- Required R packages: 
  - `spData` **2.3.3**  
  - `spdep` **1.3-3**  
  - `INLA` **24.05.10**  
  - `TruncExpFam` **1.2.0**  
