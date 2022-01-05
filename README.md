# Computational Framework for L$_{2}$E Structured Regression Problems

We introduce a user-friendly computational framework for implementing robust versions of a wide variety of structured regression methods with the L$_{2}$ criterion.  In addition to introducing an algorithm for performing L$_{2}$E regression, our framework enables robust regression with the L$_{2}$ criterion for additional structural constraints, works without requiring complex tuning procedures on the precision parameter, can be used to identify heterogeneous subpopulations, and can incorporate readily available non-robust structured regression solvers.  We provide convergence guarantees for the framework and demonstrate its flexibility with some examples.  Supplementary materials for this article are available online.

## Installation

To install the latest stable version from CRAN:

  ```{r}
install.packages('L2E-package-demo')
```

To install the latest development version from GitHub:

  ```{r}
# install.packages("devtools")
devtools::install_github('jocelynchi/L2E-package-demo')
```

## Getting Started

We've included an introductory [demo](https://jocelynchi.github.io/L2E-package-demo/articles/l2e-intro.html) on how to use the `L2E` framework with examples from the accompanying journal manuscript.

## Citing L2E

A pdf of the accompanying journal manuscript for `L2E` can be found at [arXiv:2105.03228](https://arxiv.org/abs/2010.04133).  To cite the `L2E` framework, please use the following BibTeX entry.

```
@article{l2e,
  author = {Jocelyn T. Chi and Eric C. Chi},
  title = {A User-Friendly Computational Framework for Robust Structured Regression with the L$_2$ Criterion},
  journal = {Journal of Computational and Graphical Statistics, In press.},
  year = {2022+},
  url = {https://arxiv.org/abs/2010.04133}
}
```

