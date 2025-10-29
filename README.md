This repository contains materials to reproduce the experiments and figures
appearing in 

*Multivariate Analysis for Multiple Network Data via Semi-Symmetric Tensor PCA*
by Michael Weylandt and George Michailidis. 

[ArXiv 2202.04719](https://arxiv.org/abs/2202.04719)

---

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
    │       ├── Makefile
    │       ├── README
    │       ├── figure_A5.R
    │       └── scotus_networks.csv
    │   └── finance_case_study_sec_f22
    │   │   ├── figure_A6.R
    │   │   └── Makefile
    ├── additional_simulation_studies_f3
    │   ├── README
    │   ├── multifactor_sstpca_sec_f31
    │   │   ├── Makefile
    │   │   └── figure_A7.R
    │   ├── deflation_sec_f32
    │   │   ├── Makefile
    │   │   └── figure_A8.R
    │   ├── simulate_f3_data.R
    │   └── tuning_parameters_sec_f33
    │   │   ├── Makefile
    │       ├── figure_A9.R
    │       ├── figure_A10.R
    │       ├── simulation_select_K.R
    │       └── simulation_select_rank_BIC.R
    └── Makefile
```
