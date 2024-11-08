# DSA42288S Final Year Project

This repository contains the code and results for my FYP that aims to synthesize joint distributions from marginal densities with controlled statistical properties, such as desired correlations and mutual information.


## Abstract

In many statistical and machine learning applications, controlling the relationships between random variables is essential for tasks such as data simulation, model validation, and experimental design. This project focuses on constructing joint distributions from marginal densities with specific statistical properties, including desired levels of correlation and mutual information, while preserving the integrity of the marginal distributions. 

First, by utilizing elliptical copulas, we construct joint distributions that meet specified correlation targets. We show that these results are sensitive to the methods used for approximating inverse cumulative distribution functions (CDFs). Thus, we also investigate how results vary with the choice of inverse CDF approximations such as linear interpolation, Akima interpolation, polynomial fitting, and quantile-based methods to accurately map random samples to the desired marginal distributions.

Next, this research delves into mutual information as a measure of dependency for categorical variables. An adaptive simulated annealing algorithm is developed to adjust mutual information between variables, offering precise control over their dependency structure. This approach could be useful in applications where fine-tuned associations are necessary, such as synthetic data generation for education and controlled experimentation.

Finally, the project explores how to intentionally engineer Simpson’s Paradox in continuous datasets. By manipulating correlations between subgroups of data, we induce scenarios where aggregated data reveals a different relationship than is apparent within individual groups. This has potential implications for data analysis, highlighting the dangers of misinterpretation when data is improperly aggregated.

Overall, this project provides a framework for synthesizing joint distributions from marginal densities, with applications across data science, simulation, and statistical analysis. The ability to control correlation and mutual information opens new avenues for model testing and the generation of realistic datasets with predetermined statistical properties.

## Directory Structure

The project is organized as follows:

```plaintext
.
├── README.md
├── correlation_continuous
│   ├── final_eval.R
│   ├── final_results
│   └── old
├── mutual_information_categorical
│   ├── code_SL_modified_MAX_MIN.R
│   ├── final_results
│   ├── get_desired_stepwise.R
│   └── old
├── simpsons_paradox_cont
│   ├── GA_Fns.R
│   ├── attempt_01_copulas_composition.R
│   ├── attempt_02_SA_Composition.R
│   ├── attempt_03_GA_Composition.R
│   ├── final_results
│   ├── old
│   └── param_grid_model.R
└── simpsons_paradox_disc
    ├── final_results
    ├── get_desired_SL.R
    ├── log_odds_comapre.R
    └── old
```

## Overview of Main Directories and Files

This repository is divided into three primary modules:

1. Copulas and Correlation (`correlation_continuous`) This module explores the creation of joint distributions that meet specified correlation targets using elliptical copulas. It also examines how different methods for approximating inverse cumulative distribution functions affect the results.  

2. Mutual Information for Contingency Tables (`mutual_information_categorical`): This module focuses on adjusting mutual information between categorical variables, useful for tasks where precise dependency control is essential, such as synthetic data generation and educational examples.

3. Simpson’s Paradox: A Case Study, Divided into two sub-modules:
- Continuous Variables (`simpsons_paradox_cont`): Engineering Simpson’s Paradox in continuous datasets faceted by a categorical variable.  
- Categorical Variables (`simpsons_paradox_disc`): Engineering Simpson’s Paradox in categorical datasets faceted by a categorical variable.  

---


## Detailed Description of Each Directory and File

### 1. correlation_continuous

This module contains the code for constructing joint distributions with specific correlation structures using copulas.

- `final_eval.R`: This file includes all evaluation functions and implementations for various inverse CDF approximation methods, such as linear and Akima interpolation, polynomial fitting, and quantile-based methods. It also contains metrics for evaluating these methods.
- `final_results/`: Contains images summarizing the results of experiments, including correlation vs. smoothness trade-offs, error plots, and metric-specific plots (e.g., KS, CVM, SM, TV).
- `old/`: Stores legacy code and previous results for reference, including implementations of various copulas (e.g., Archimedean copulas), density estimation methods, and early experiments.

### 2. mutual_information_categorical

This module focuses on maximizing or minimizing mutual information between categorical variables using different optimization techniques.

- `code_SL_modified_MAX_MIN.R`: Contains results from maximization and minimization experiments on Simulated Annealing using various approaches, including small optimization steps and Sinkhorn distance calculations.
- `final_results/`: Contains final images summarizing the results of these experiments, such as mutual information increase and decrease plots, as well as maximization plots for different contingency table sizes.  
- `get_desired_stepwise.R`: Implements a stepwise algorithm to achieve the desired mutual information between variables, allowing for precise control over dependency structures.  
- `old/`: Stores older code versions and experiment results, including separate maximization and minimization files, modified Sinkhorn implementations, and initial experiments.  

### 3. simpsons_paradox_cont

This module investigates Simpson’s Paradox within continuous datasets. It employs Genetic Algorithms (GA) and Simulated Annealing (SA) to engineer datasets where aggregated and subgroup-level relationships differ.

- `GA_Fns.R`: Contains all Genetic Algorithm functions, including mutation, crossover, and selection, which are used in `attempt_03_GA_Composition.R`.  
- `attempt_01_copulas_composition.R`: The initial approach to engineering Simpson’s Paradox by modifying correlation levels for each category, inspired by the composition method.  
- `attempt_02_SA_Composition.R`: Builds on the first attempt by introducing Simulated Annealing to reduce separation between categories.  
- `attempt_03_GA_Composition.R`: Adds Genetic Algorithm-based techniques to further refine the separation-reduction process.  
- `final_results/`: Contains images with final results, such as plots demonstrating Simpson’s Paradox by adjusting correlations within subgroups.  
- `old/`: Stores older files and results for reference, including various experiments and parameter adjustments.  
- `param_grid_model.R`: Code for experimenting with a parameteric grids model to optimize data point allocations for preserving initial correlation structure.  

## 4. simpsons_paradox_disc

This module focuses on creating Simpson’s Paradox with categorical variables. It employs a log-odds approach to compare different aggregation models.

- `final_results/`: Contains final images summarizing the results of experiments, including comparison plots for weighted log-odds and normalized-count weighted effects.  
- `get_desired_SL.R`: Aims to generate Simpson’s Paradox using Agresti’s log-odds measure, showing differences in effects when data is aggregated.  
- `log_odds_comapre.R`: Compares Agresti’s log-odds against weighted log-odds using various aggregation methods, offering insights into how different data aggregation affects the paradox.  
- `old/`: Stores older files and initial experiments related to Simpson’s Paradox with categorical data.  



## Additional Notes

- Final Results Folders: Each module contains a final_results folder where images and plots summarizing the experimental results are stored.

- Old Code: Each module has an old subdirectory containing legacy code and preliminary experiment results. These may serve as references for understanding the evolution of methods used in this project.  

## Acknowledgments

This project was developed as part of the DSA42288S Final Year Project.  I would like to express my  gratitude to my supervisor, Dr. Vikneswaran Gopal, for his invaluable guidance, support, and mentorship throughout this project. I am also grateful to the faculty and staff at NUS for their continuous support, and to my family for their encouragement along the way.


## Contact

For any questions or inquiries, please contact me at `naman.agr03@gmail.com`.


