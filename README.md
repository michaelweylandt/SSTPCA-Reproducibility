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
.
├── LICENSE
├── README.md
├── main_paper
│   ├── case_study_sec_42
│   └── simulation_studies_sec_41
│       ├── simulation_01
│       ├── simulation_02
│       └── simulation_03
└── supplements
    ├── additional_case_studies_f2
    │   ├── finance_case_study_sec_f22
    │   └── scotus_case_study_sec_f21
    ├── additional_simulation_studies_f3
    │   ├── deflation_sec_f32
    │   ├── multifactor_sstpca_sec_f31
    │   └── tuning_parameters_sec_f33
    └── hs_case_study_additional_details_sec_f1
```
