# Reproduction Material for Weylandt and Michailidis (2025+)

This repository contains materials to reproduce the experiments and figures
appearing in 

*Multivariate Analysis for Multiple Network Data via Semi-Symmetric Tensor PCA*
by Michael Weylandt and George Michailidis. [ArXiv 2202.04719](https://arxiv.org/abs/2202.04719)

To reproduce all analysis, simply run `make` at the top level of this directory.
To reproduce a particular figure or experiment, navigate to the appropriate
subdirectory and run `make` there. 

The organization of this repository mirrors the organization of the paper: 

```
├── LICENSE
├── README.md
├── Makefile
├── main_paper
│   ├── Makefile
│   ├── case_study_sec_42
│   │   ├── figure_6.R
│   │   └── Makefile
│   └── simulation_studies_sec_41
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
    │   ├── scotus_case_study_sec_f21
    │   │   ├── Makefile
    │   │   ├── README
    │   │   ├── figure_A5.R
    │   │   └── scotus_networks.csv
    │   └── finance_case_study_sec_f22
    │       ├── Makefile
    │       └── figure_A6.R
    ├── additional_simulation_studies_f3
    │   ├── README
    │   ├── simulate_f3_data.R
    │   ├── multifactor_sstpca_sec_f31
    │   │   ├── Makefile
    │   │   └── figure_A7.R
    │   ├── deflation_sec_f32
    │   │   ├── Makefile
    │   │   └── figure_A8.R
    │   └── tuning_parameters_sec_f33
    │       ├── Makefile
    │       ├── figure_A9.R
    │       ├── figure_A10.R
    │       ├── simulation_select_K.R
    │       └── simulation_select_rank_BIC.R
    └── Makefile
```

No code is provided for Figure 1 as it is a schematic drawing of our method
and not the result of any analysis or simulation. 

# Structure

## Figures

Each `figure_*.R` file contains the code to reproduce that figure. It should
be possible to run each file directly, though using that directory's 
`Makefile` is recommended as it will run any needed simulations first. 

The figures reproduced here are: 

- Figure 1: Schematic Diagram of the SST-PCA Decomposition (not included)
- [Figure 2](main_paper/simulation_studies_sec41/simulation_01/figure_2.R): A
  comparison of SST-PCA with some other non-network specific tensor decomposition
  algorithms under a set of popular low-rank graph models.
- [Figure 3](main_paper/simulation_studies_sec41/simulation_02/figure_3.R): 
  Demonstration of the rapid computational convergence of the alternating 
  SST-PCA algorithm.
- [Figure 4](main_paper/simulation_studies_sec41/simulation_03/figure_4-5.R):
  A 'statistical convergence' simulation demonstrating that SST-PCA i) has
  error decreasing inversely to signal strength and ii) scales
  relatively well with the number of nodes in a network. This plot also
  demonstrates that SST-PCA is rather robust to the choice of initialization.
- [Figure 5](main_paper/simulation_studies_sec41/simulation_03/figure_4-5.R):
  A 'statistical convergence' simulation demonstrating that SST-PCA i) has
  error decreasing inversely to number of networks observed. Like Figure 4, this 
  plot also demonstrates that SST-PCA is rather robust to the choice of 
  initialization.
- [Figure 6](main_paper/case_study_sec42/figure_6.R): A demonstration on SST-PCA
  on networks constructed from the school-day contacts of a set of French high
  school students. SST-PCA is able to accurately recover the low-rank structure
  (corresponding to different tracks of study) in this data. The 'time factor'
  estimated by SST-PCA also reveals intra-day patterns of interest. 
  
The Supplemental Materials have several additional figures: 

- [Figure A1](supplements/hs_case_study_additional_details_sec_f1/figure_A1.R):
  Application of SST-PCA to the _first moment_ of the French High School Data. 
  While PCA is typically applied to estimate the variance structure of data, it
  can also be used to compute a low-rank approximation to the mean of 
  high-dimensional data. We demonstrate the ability of SST-PCA to be used in this
  fashion. 
- [Figure A2](supplements/hs_case_study_additional_details_sec_f1/figure_A2.R): 
  A more detailed study of the time factor estimated by SST-PCA on the French 
  high school data. Here, we demonstrate the SST-PCA is able to flexibly adapt 
  to periods of different signal or edge density. 
- [Figure A3](supplements/hs_case_study_additional_details_sec_f1/figure_A3-A4.R):
  A comparison of SST-PCA against other methods on the French High School data
  clustering task. Each method is applied to estimate a cluster structure of
  the students, which is then compared against a 'ground-truth' based on the
  students' educational specializations. SST-PCA is compared against 
  the COSIE framework of Arroyo *et al*, the JEG method of
  Wang *et al*, and the TWIST approach of Jing *et al*, as well as several
  simple forms of temporal aggregation (i.e., counting or averaging edges over 
  time). 
- [Figure A4](supplements/hs_case_study_additional_details_sec_f1/figure_A3-A4.R)
  Visualization and comparisons of the cluster structure inferred by SST-PCA, 
  TWIST, COSIE, and JEG, each of which attempt to estimate some form of low-rank 
  'principal network.'
- [Figure A5](supplements/additional_case_studies_sec_f2/scotus_case_study_sec_f21/figure_A5.R):
  Application of SST-PCA to a data set extracted from the voting patterns of US
  Supreme Court Justices. While this data set is small (only $p=9$ Justices), 
  SST-PCA is able to find interesting "first moment" and "second moment" structure.
  Additionally, we combine SST-PCA with Wang and Samworth's method for time series
  changepoint detection to identify the most important change in Supreme Court
  dynamics in the 1995-2020 study period. 
- [Figure A6](supplements/additional_case_studies_sec_f2/finance_case_study_sec_f22/figure_A6.R):
  Application of SST-PCA to a data set extracted from the correlation of various
  national stock indices. Here, we apply SST-PCA to identify structural 
  correlations ("market beta"), regional effects, and the impact of the European
  Debt Crisis on the interconnectivity patterns of financial markets. 
- [Figure A7](supplements/additional_simulation_studies_sec_f3/multifactor_sstpca_sec_f31/figure_A7.R):
  A simulation study showing how SST-PCA can be iteratively re-applied to 
  estimate multiple SST-PCA factors in a greedy fashion. This simulation also
  demonstrates a moderate degree of robustness of SST-PCA to non-orthogonality
  of the underlying factors. 
- [Figure A8](supplements/additional_simulation_studies_sec_f3/deflation_sec_f32/figure_A8.R):
  A simulation study demonstrating the impact of different 
  deflation strategies. Depending on the orthogonality (or lack of orthogonality)
  of the underlying low-rank model, different deflation methods may be preferred
  when estimating multiple SST-PCA factors. 
- [Figure A9](supplements/additional_simulation_studies_sec_f3/tuning_parameters_sec_f33/figure_A9.R):
  A simulation study illustrating the effectiveness of a singular-value thresholding
  inspired approach to selecting the number of SST-PCA factors. 
- [Figure A10](supplements/additional_simulation_studies_sec_f3/tuning_parameters_sec_f33/figure_A10.R):
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
  [the code generating Figure 6](main_paper/case_study_sec42/figure_6.R) for an
  example. 

- [SCOTUS Voting Data Set](supplements/additional_case_studies_f2/scotus_case_study_sec_f21/scotus_network.csv)
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
  [replicating Figure A6](supplements/additional_case_studies_f2/finance_case_study_sec_f22/figure_A6.R).
  
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
   
