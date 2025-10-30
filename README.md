# Reproduction Material for Weylandt and Michailidis (2025+)

This repository contains materials to reproduce the experiments and figures
appearing in 

*Multivariate Analysis for Multiple Network Data via Semi-Symmetric Tensor PCA*
by Michael Weylandt and George Michailidis. [ArXiv 2202.04719](https://arxiv.org/abs/2202.04719)

To reproduce all analysis, simply run `make` at the top level of this directory.
To reproduce a particular figure or experiment, navigate to the appropriate
sub-directory and run `make` there. 

The organization of this repository mirrors the organization of the paper: 

```
├── LICENSE
├── README.md
├── Makefile
├── main_paper
│   ├── Makefile
│   ├── case_study_42
│   │   ├── figure_6.R
│   │   └── Makefile
│   └── simulation_studies_41
│       ├── Makefile
│       ├── simulation_01
│       │   ├── Makefile
│       │   ├── figure_2.R
│       │   ├── simulation_methods_comparison_results.rds
│       │   └── simulation_methods_figure_2.R
│       ├── simulation_02
│       │   ├── Makefile
│       │   ├── figure_3.R
│       │   ├── simulation_computational_convergence_figure_3.R
│       │   └── simulation_computational_convergence_results.rds
│       └── simulation_03
│           ├── Makefile
│           ├── figure_4-5.R
│           ├── simulation_1.R
│           ├── simulation_1.rds
│           ├── simulation_2.R
│           └── simulation_2.rds
└── supplements
    ├── additional_case_studies_f2
    │   ├── Makefile
    │   ├── scotus_case_study_f21
    │   │   ├── Makefile
    │   │   ├── README
    │   │   ├── figure_A5.R
    │   │   └── scotus_networks.csv
    │   └── finance_case_study_f22
    │       ├── Makefile
    │       └── figure_A6.R
    ├── additional_simulation_studies_f3
    │   ├── README
    │   ├── simulate_f3_data.R
    │   ├── multifactor_sstpca_f31
    │   │   ├── Makefile
    │   │   └── figure_A7.R
    │   ├── deflation_f32
    │   │   ├── Makefile
    │   │   └── figure_A8.R
    │   └── tuning_parameters_f33
    │       ├── Makefile
    │       ├── figure_A9.R
    │       ├── figure_A10.R
    │       ├── simulation_select_K.R
    │       └── simulation_select_rank_BIC.R
    └── Makefile
```

Beyond the top-level directory, all sub-directories end with the name of the 
corresponding section of the paper where the figure appears: *e.g.*, 
`case_study_42` contains code necessary to reproduce the figures appearing in 
Section 4.2 of the main manuscript, while `tuning_parameters_f33` contains
code corresponding to Section F.3.3 of the Supplementary Materials.

No code is provided for Figure 1 as it is a schematic drawing of our method
and not the result of any analysis or simulation. 

# Structure

## Figures

Each `figure_*.R` file contains the code to reproduce that figure. It should
be possible to run each file directly, though using that directory's 
`Makefile` is recommended as it will run any needed simulations first. 

The figures reproduced here are: 

- Figure 1: Schematic Diagram of the SST-PCA Decomposition (not included)
- [Figure 2](main_paper/simulation_studies_41/simulation_01/figure_2.R): A
  comparison of SST-PCA with some other non-network specific tensor decomposition
  algorithms under a set of popular low-rank graph models.
- [Figure 3](main_paper/simulation_studies_41/simulation_02/figure_3.R): 
  Demonstration of the rapid computational convergence of the alternating 
  SST-PCA algorithm.
- [Figure 4](main_paper/simulation_studies_41/simulation_03/figure_4-5.R):
  A 'statistical convergence' simulation demonstrating that SST-PCA i) has
  error decreasing inversely to signal strength and ii) scales
  relatively well with the number of nodes in a network. This plot also
  demonstrates that SST-PCA is rather robust to the choice of initialization.
- [Figure 5](main_paper/simulation_studies_41/simulation_03/figure_4-5.R):
  A 'statistical convergence' simulation demonstrating that SST-PCA i) has
  error decreasing inversely to number of networks observed. Like Figure 4, this 
  plot also demonstrates that SST-PCA is rather robust to the choice of 
  initialization.
- [Figure 6](main_paper/case_study_42/figure_6.R): A demonstration on SST-PCA
  on networks constructed from the school-day contacts of a set of French high
  school students. SST-PCA is able to accurately recover the low-rank structure
  (corresponding to different tracks of study) in this data. The 'time factor'
  estimated by SST-PCA also reveals intra-day patterns of interest. 
  
The Supplemental Materials have several additional figures: 

- [Figure A1](supplements/hs_case_study_additional_details_f1/figure_A1.R):
  Application of SST-PCA to the _first moment_ of the French High School Data. 
  While PCA is typically applied to estimate the variance structure of data, it
  can also be used to compute a low-rank approximation to the mean of 
  high-dimensional data. We demonstrate the ability of SST-PCA to be used in this
  fashion. 
- [Figure A2](supplements/hs_case_study_additional_details_f1/figure_A2.R): 
  A more detailed study of the time factor estimated by SST-PCA on the French 
  high school data. Here, we demonstrate the SST-PCA is able to flexibly adapt 
  to periods of different signal or edge density. 
- [Figure A3](supplements/hs_case_study_additional_details_f1/figure_A3-A4.R):
  A comparison of SST-PCA against other methods on the French High School data
  clustering task. Each method is applied to estimate a cluster structure of
  the students, which is then compared against a 'ground-truth' based on the
  students' educational specializations. SST-PCA is compared against 
  the COSIE framework of Arroyo *et al*, the JEG method of
  Wang *et al*, and the TWIST approach of Jing *et al*, as well as several
  simple forms of temporal aggregation (i.e., counting or averaging edges over 
  time). 
- [Figure A4](supplements/hs_case_study_additional_details_f1/figure_A3-A4.R)
  Visualization and comparisons of the cluster structure inferred by SST-PCA, 
  TWIST, COSIE, and JEG, each of which attempt to estimate some form of low-rank 
  'principal network.'
- [Figure A5](supplements/additional_case_studies_f2/scotus_case_study_f21/figure_A5.R):
  Application of SST-PCA to a data set extracted from the voting patterns of US
  Supreme Court Justices. While this data set is small (only $p=9$ Justices), 
  SST-PCA is able to find interesting "first moment" and "second moment" structure.
  Additionally, we combine SST-PCA with Wang and Samworth's method for time series
  changepoint detection to identify the most important change in Supreme Court
  dynamics in the 1995-2020 study period. 
- [Figure A6](supplements/additional_case_studies_f2/finance_case_study_f22/figure_A6.R):
  Application of SST-PCA to a data set extracted from the correlation of various
  national stock indices. Here, we apply SST-PCA to identify structural 
  correlations ("market beta"), regional effects, and the impact of the European
  Debt Crisis on the interconnectivity patterns of financial markets. 
- [Figure A7](supplements/additional_simulation_studies_f3/multifactor_sstpca_f31/figure_A7.R):
  A simulation study showing how SST-PCA can be iteratively re-applied to 
  estimate multiple SST-PCA factors in a greedy fashion. This simulation also
  demonstrates a moderate degree of robustness of SST-PCA to non-orthogonality
  of the underlying factors. 
- [Figure A8](supplements/additional_simulation_studies_f3/deflation_f32/figure_A8.R):
  A simulation study demonstrating the impact of different 
  deflation strategies. Depending on the orthogonality (or lack of orthogonality)
  of the underlying low-rank model, different deflation methods may be preferred
  when estimating multiple SST-PCA factors. 
- [Figure A9](supplements/additional_simulation_studies_f3/tuning_parameters_f33/figure_A9.R):
  A simulation study illustrating the effectiveness of a singular-value thresholding
  inspired approach to selecting the number of SST-PCA factors. 
- [Figure A10](supplements/additional_simulation_studies_f3/tuning_parameters_f33/figure_A10.R):
  A simulation study illustrating the effectiveness of BIC-based selection of the
  rank of the low-rank component of the SST-PCA principal network.
  While BIC-based rank selection does not quite match oracle performance, the
  difference is minor in moderate signal-to-noise regimes.
  
## Data

Three non-simulated data sets are used in this paper: 

- [French High School Students](http://www.sociopatterns.org/datasets/high-school-contact-and-friendship-networks/): 
  This data set, published by [SocioPatterns](http://www.sociopatterns.org)
  captures the interactions of a set of French high school students over a week
  in 2013. The students' fields of specialization are used as ground truth for
  a clustering task. This data set is used for the main case study of the paper,
  appearing in Figure 6 and Figures A1-A4.

  Reproduction code relying on this data set automatically downloads it. See
  [the code generating Figure 6](main_paper/case_study_42/figure_6.R) for an
  example. 

- [SCOTUS Voting Data Set](supplements/additional_case_studies_f2/scotus_case_study_f21/scotus_networks.csv)
  This data set was manually extracted by the authors, using data from 
  [SCOTUSBlog](https://www.scotusblog.com). Each edge captures the frequency
  of cases in which a pair of US Supreme Court Justices concurred in the judgement
  during a single October-to-September term.
  
  In practice, two Justices may reach the same decision with different reasoning
  and express their alternative rationale in a concurrence. For simplicity, we
  focus only on the actual decision and do not attempt to parse the differences
  in different judicial opinions.
  
- Stock Index Correlations: This data set is obtained by downloading the price
  history of the following Exchange-Traded Funds (ETFs) that track various
  important stock indices. Daily adjusted (total) returns are calculated and used
  to estimate the strength of connection between different stock markets (and by
  extension, national economies). The data is downloaded from Yahoo! Finance
  using the [`quantmod` `R` package](https://www.quantmod.com/). This data
  is downloaded automatically as part of 
  [replicating Figure A6](supplements/additional_case_studies_f2/finance_case_study_f22/figure_A6.R).
  
  The following ETFs are used: 
  
| Country         | ETF |
|:---------------:|:---:|
|Japan            |EWJ  |
|Taiwan           |EWT  |
|South Korea      |EWY  |
|China            |FXI  |
|Brazil           |EWZ  |
|Canada           |EWC  |
|United Kingdom   |EWU  |
|Germany          |EWG  |
|Switzerland      |EWL  |
|Australia        |EWA  |
|Mexico           |EWW  |
|Hong Kong        |EWH  |
|India 50         |INDY |
|Sweden           |EWD  |
|France           |EWQ  |
|Spain            |EWP  |
|Russia           |ERUS |
|Singapore        |EWS  |
|Italy            |EWI  |
|Indonesia        |EIDO |
|Chile            |ECH  |
|Thailand         |THD  |
|Netherlands      |EWN  |
|Poland           |EPOL |
|South Africa     |EZA  |
|Malaysia         |EWM  |
|Turkey           |TUR  |
|Israel           |EIS  |
|New Zealand      |ENZL |
|Peru             |EPU  |
|Philippines      |EPHE |
|Austria          |EWO  |
|Brazil Small-Cap |EWZS |
|Ireland          |EIRL |
|Belgium          |EWK  |
|United States    |IVV  |

## Functions of Interest

While SST-PCA is straightforward to implement 'from scratch', authors
interested in applying SST-PCA may be interested in functionality contained in
[`tensor_factorizations.cpp`](tensor_factorizations.cpp). This file, which
can be easily integrated into `R` using the 
[`RcppArmadillo` package](https://doi.org/10.1016/j.csda.2013.02.005) defines
two principal functions of interest: 

1. `ss_tpm`: This function implements the core SST-PCA decomposition, with
   the following arguments: 
  
   - `X`: the semi-symmetric tensor to decompose. This should be a p-by-p-by-T
     tensor symmetric along each slice: *i.e.*, `all(X[,,t] == t(X[,,t]))` for
     all `t`. Integer or binary tensors will be silently coerced to doubles for
     internal computations.
   - `u_init_strategy`: a flag controlling the initialization strategy used for
     $u^{(0)}$: 
      - `u_init_strategy = 0`: [Default] "Stable initialization" with 
        $u^{(0)} \propto \mathbf{1}$, suitable for network series in which the
        principal network is expected to be present in each "slice", with different
        but consistently positive, weights
      - `u_init_strategy = 1`: "Random initialization": each element of $u^{(0)}$
        is independently standard Gaussian, giving $u^{(0)}$ selected uniformly
        at random (Haar measure) over the unit sphere
      - `u_init_strategy = 2`: $u^{(0)}$ is initialized to a user-provided vector
        passed as the optional argument `u_init`
      - Other values are not allowed. 
   - `rank`: the rank of the principal network $\hat{V}$ estimated by SST-PCA.
      Default 1. Must be a positive integer. 
   - `eps`: the tolerance used to check convergence, stopping when 
      $|u^{(k)} - u^{(k-1)}| + |v^{(k)} - v^{(k-1)}|$ falls below `eps`. Default `1e-6`.
      Must be a positive real number. 
   - `max_iter`: the maximum number of SST-PCA iterations before stopping.
      Default 1000. Must be a positive integer.
   - `u_init`: An initial value to be used for $u^{(0)}$. Required when 
      `u_init_strategy=2`. If passed, this must be a vector of length `T`, where
      `T` is the third dimension of `X`. 

2. `ss_tpm_large`: An experimental alternative to SST-PCA which uses an iterative
   power method for the $V$-update, rather than a full eigendecomposition. For 
   very large networks, this may be faster, though the convergence properties of
   this approach are less well understood.
   
Both functions return a list containing: 

- `u_hat`: the estimated SST-PCA time factor (T-vector)
- `v_hat`: the orthogonal matrix factor of the SST-PCA principal network 
  (p-by-rank matrix)
- `V_hat`: the SST-PCA principal network, formed by taking the outer product of
  `v_hat` with itself (p-by-p matrix of specified rank)
- `d`: the SST-PCA scaling factor, roughly equivalent to a singular value 
  (positive scalar)
- `X_hat`: the SST-PCA approximation to the input tensor `X`, formed by the 
  outer product of `u_hat` and `V_hat`, scaled by `d` (semi-symmetric tensor
  of dimensions matching the original input)
   
## Version Information

The following version of `R` and associated packages were used to generate
these figures: 

```
─ Session info ──────────────────────────────────────────────────────
 setting  value
 version  R version 4.5.1 (2025-06-13)
 os       macOS Sequoia 15.6.1
 system   aarch64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/New_York
 date     2025-10-29
 rstudio  2025.09.0+387 Cucumberleaf Sunflower (desktop)
 pandoc   3.7.0.2 @ /opt/homebrew/bin/ (via rmarkdown)
 quarto   1.7.33 @ /usr/local/bin/quarto

─ Packages ──────────────────────────────────────────────────────────
 package       * version  date (UTC) lib source
 askpass         1.2.1    2024-10-04 [1] CRAN (R 4.5.0)
 backports       1.5.0    2024-05-23 [1] CRAN (R 4.5.0)
 base64enc       0.1-3    2015-07-28 [1] CRAN (R 4.5.0)
 bit             4.6.0    2025-03-06 [1] CRAN (R 4.5.0)
 bit64           4.6.0-1  2025-01-16 [1] CRAN (R 4.5.0)
 blob            1.2.4    2023-03-17 [1] CRAN (R 4.5.0)
 broom           1.0.10   2025-09-13 [1] CRAN (R 4.5.0)
 bslib           0.9.0    2025-01-30 [1] CRAN (R 4.5.0)
 cachem          1.1.0    2024-05-16 [1] CRAN (R 4.5.0)
 callr           3.7.6    2024-03-25 [1] CRAN (R 4.5.0)
 cellranger      1.1.0    2016-07-27 [1] CRAN (R 4.5.0)
 cli             3.6.5    2025-04-23 [1] CRAN (R 4.5.0)
 clipr           0.8.0    2022-02-22 [1] CRAN (R 4.5.0)
 conflicted      1.2.0    2023-02-01 [1] CRAN (R 4.5.0)
 cpp11           0.5.2    2025-03-03 [1] CRAN (R 4.5.0)
 crayon          1.5.3    2024-06-20 [1] CRAN (R 4.5.0)
 curl            7.0.0    2025-08-19 [1] CRAN (R 4.5.0)
 data.table      1.17.8   2025-07-10 [1] CRAN (R 4.5.0)
 DBI             1.2.3    2024-06-02 [1] CRAN (R 4.5.0)
 dbplyr          2.5.1    2025-09-10 [1] CRAN (R 4.5.0)
 digest          0.6.37   2024-08-19 [1] CRAN (R 4.5.0)
 dplyr         * 1.1.4    2023-11-17 [1] CRAN (R 4.5.0)
 dtplyr          1.3.2    2025-09-10 [1] CRAN (R 4.5.0)
 evaluate        1.0.5    2025-08-27 [1] CRAN (R 4.5.0)
 farver          2.1.2    2024-05-13 [1] CRAN (R 4.5.0)
 fastmap         1.2.0    2024-05-15 [1] CRAN (R 4.5.0)
 fontawesome     0.5.3    2024-11-16 [1] CRAN (R 4.5.0)
 forcats       * 1.0.0    2023-01-29 [1] CRAN (R 4.5.0)
 fs              1.6.6    2025-04-12 [1] CRAN (R 4.5.0)
 gargle          1.6.0    2025-09-03 [1] CRAN (R 4.5.0)
 generics        0.1.4    2025-05-09 [1] CRAN (R 4.5.0)
 ggforce         0.5.0    2025-06-18 [1] CRAN (R 4.5.0)
 ggnewscale      0.5.2    2025-06-20 [1] CRAN (R 4.5.0)
 ggplot2       * 4.0.0    2025-09-11 [1] CRAN (R 4.5.0)
 ggraph          2.2.2    2025-08-24 [1] CRAN (R 4.5.0)
 ggrepel         0.9.6    2024-09-07 [1] CRAN (R 4.5.0)
 glue            1.8.0    2024-09-30 [1] CRAN (R 4.5.0)
 googledrive     2.1.2    2025-09-10 [1] CRAN (R 4.5.0)
 googlesheets4   1.1.2    2025-09-03 [1] CRAN (R 4.5.0)
 graphlayouts    1.2.2    2025-01-23 [1] CRAN (R 4.5.0)
 gridExtra       2.3      2017-09-09 [1] CRAN (R 4.5.0)
 gtable          0.3.6    2024-10-25 [1] CRAN (R 4.5.0)
 haven           2.5.5    2025-05-30 [1] CRAN (R 4.5.0)
 highr           0.11     2024-05-26 [1] CRAN (R 4.5.0)
 hms             1.1.3    2023-03-21 [1] CRAN (R 4.5.0)
 htmltools       0.5.8.1  2024-04-04 [1] CRAN (R 4.5.0)
 httr            1.4.7    2023-08-15 [1] CRAN (R 4.5.0)
 ids             1.0.1    2017-05-31 [1] CRAN (R 4.5.0)
 igraph          2.1.4    2025-01-23 [1] CRAN (R 4.5.0)
 isoband         0.2.7    2022-12-20 [1] CRAN (R 4.5.0)
 jquerylib       0.1.4    2021-04-26 [1] CRAN (R 4.5.0)
 jsonlite        2.0.0    2025-03-27 [1] CRAN (R 4.5.0)
 knitr           1.50     2025-03-16 [1] CRAN (R 4.5.0)
 labeling        0.4.3    2023-08-29 [1] CRAN (R 4.5.0)
 lattice         0.22-7   2025-04-02 [1] CRAN (R 4.5.1)
 lifecycle       1.0.4    2023-11-07 [1] CRAN (R 4.5.0)
 lubridate     * 1.9.4    2024-12-08 [1] CRAN (R 4.5.0)
 magrittr        2.0.4    2025-09-12 [1] CRAN (R 4.5.0)
 MASS            7.3-65   2025-02-28 [1] CRAN (R 4.5.1)
 Matrix          1.7-4    2025-08-28 [1] CRAN (R 4.5.0)
 memoise         2.0.1    2021-11-26 [1] CRAN (R 4.5.0)
 mime            0.13     2025-03-17 [1] CRAN (R 4.5.0)
 modelr          0.1.11   2023-03-22 [1] CRAN (R 4.5.0)
 openssl         2.3.3    2025-05-26 [1] CRAN (R 4.5.0)
 patchwork       1.3.2    2025-08-25 [1] CRAN (R 4.5.0)
 pillar          1.11.0   2025-07-04 [1] CRAN (R 4.5.0)
 pkgconfig       2.0.3    2019-09-22 [1] CRAN (R 4.5.0)
 plyr            1.8.9    2023-10-02 [1] CRAN (R 4.5.0)
 polyclip        1.10-7   2024-07-23 [1] CRAN (R 4.5.0)
 prettyunits     1.2.0    2023-09-24 [1] CRAN (R 4.5.0)
 processx        3.8.6    2025-02-21 [1] CRAN (R 4.5.0)
 progress        1.2.3    2023-12-06 [1] CRAN (R 4.5.0)
 ps              1.9.1    2025-04-12 [1] CRAN (R 4.5.0)
 purrr         * 1.1.0    2025-07-10 [1] CRAN (R 4.5.0)
 quantmod      * 0.4.28   2025-06-19 [1] CRAN (R 4.5.0)
 R6              2.6.1    2025-02-15 [1] CRAN (R 4.5.0)
 ragg            1.5.0    2025-09-02 [1] CRAN (R 4.5.0)
 rappdirs        0.3.3    2021-01-31 [1] CRAN (R 4.5.0)
 RColorBrewer    1.1-3    2022-04-03 [1] CRAN (R 4.5.0)
 Rcpp            1.1.0    2025-07-02 [1] CRAN (R 4.5.0)
 RcppArmadillo   15.0.2-1 2025-09-08 [1] CRAN (R 4.5.0)
 readr         * 2.1.5    2024-01-10 [1] CRAN (R 4.5.0)
 readxl          1.4.5    2025-03-07 [1] CRAN (R 4.5.0)
 rematch         2.0.0    2023-08-30 [1] CRAN (R 4.5.0)
 rematch2        2.1.2    2020-05-01 [1] CRAN (R 4.5.0)
 reprex          2.1.1    2024-07-06 [1] CRAN (R 4.5.0)
 reshape2        1.4.4    2020-04-09 [1] CRAN (R 4.5.0)
 rlang           1.1.6    2025-04-11 [1] CRAN (R 4.5.0)
 rmarkdown       2.29     2024-11-04 [1] CRAN (R 4.5.0)
 rstudioapi      0.17.1   2024-10-22 [1] CRAN (R 4.5.0)
 rTensor         1.4.9    2025-08-25 [1] CRAN (R 4.5.0)
 rvest           1.0.5    2025-08-29 [1] CRAN (R 4.5.0)
 S7              0.2.0    2024-11-07 [1] CRAN (R 4.5.0)
 sass            0.4.10   2025-04-11 [1] CRAN (R 4.5.0)
 scales          1.4.0    2025-04-24 [1] CRAN (R 4.5.0)
 selectr         0.4-2    2019-11-20 [1] CRAN (R 4.5.0)
 stringi         1.8.7    2025-03-27 [1] CRAN (R 4.5.0)
 stringr       * 1.5.2    2025-09-08 [1] CRAN (R 4.5.0)
 sys             3.4.3    2024-10-04 [1] CRAN (R 4.5.0)
 systemfonts     1.3.1    2025-10-01 [1] CRAN (R 4.5.0)
 textshaping     1.0.3    2025-09-02 [1] CRAN (R 4.5.0)
 tibble        * 3.3.0    2025-06-08 [1] CRAN (R 4.5.0)
 tidygraph       1.3.1    2024-01-30 [1] CRAN (R 4.5.0)
 tidyr         * 1.3.1    2024-01-24 [1] CRAN (R 4.5.0)
 tidyselect      1.2.1    2024-03-11 [1] CRAN (R 4.5.0)
 tidyverse     * 2.0.0    2023-02-22 [1] CRAN (R 4.5.0)
 timechange      0.3.0    2024-01-18 [1] CRAN (R 4.5.0)
 tinytex         0.57     2025-04-15 [1] CRAN (R 4.5.0)
 TTR           * 0.24.4   2023-11-28 [1] CRAN (R 4.5.0)
 tweenr          2.0.3    2024-02-26 [1] CRAN (R 4.5.0)
 tzdb            0.5.0    2025-03-15 [1] CRAN (R 4.5.0)
 utf8            1.2.6    2025-06-08 [1] CRAN (R 4.5.0)
 uuid            1.2-1    2024-07-29 [1] CRAN (R 4.5.0)
 vctrs           0.6.5    2023-12-01 [1] CRAN (R 4.5.0)
 viridis         0.6.5    2024-01-29 [1] CRAN (R 4.5.0)
 viridisLite     0.4.2    2023-05-02 [1] CRAN (R 4.5.0)
 vroom           1.6.5    2023-12-05 [1] CRAN (R 4.5.0)
 withr           3.0.2    2024-10-28 [1] CRAN (R 4.5.0)
 xfun            0.53     2025-08-19 [1] CRAN (R 4.5.0)
 xml2            1.4.0    2025-08-20 [1] CRAN (R 4.5.0)
 xts           * 0.14.1   2024-10-15 [1] CRAN (R 4.5.0)
 yaml            2.3.10   2024-07-26 [1] CRAN (R 4.5.0)
 zoo           * 1.8-14   2025-04-10 [1] CRAN (R 4.5.0)

 [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
 * ── Packages attached to the search path.

─────────────────────────────────────────────────────────────────────
```
