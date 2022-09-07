# Computational Framework for L$_{2}$E Structured Regression Problems

The `L2E` package (version 2.0) implements the computational framework for L$_2$E regression in Liu, Chi, and Lange (2022+), which was built on the previous work in Chi and Chi (2022). Both works employ the block coordinate descent  strategy to solve a nonconvex optimization problem but utilize different methods for the inner block descent updates. We refer to the method in Liu, Chi, and Lange (2022+) as "MM" and the one in Chi and Chi (2022) as "PG" in our package. This package provides code to replicate some examples illustrating the usage of the frameworks in both manuscripts.


## Installation

To install the latest stable version from CRAN:

  ```{r}
install.packages('L2E')
```

To install the latest development version from GitHub:

  ```{r}
# install.packages("devtools")
devtools::install_github('jocelynchi/L2E-package-demo')
```

## Getting Started

We've included an introductory [demo](https://jocelynchi.github.io/L2E-package-demo/articles/l2e-intro.html) on how to use the `L2E` framework with examples from the accompanying journal manuscripts.

## Citing the package

Please reference the following manuscripts when citing this package.  Thank you!

```

@article{L2E-Chi,
  title={A User-Friendly Computational Framework for Robust Structured Regression with the L$_2$ Criterion},
  author={Chi, Jocelyn T. and Chi, Eric C.},
  journal={Journal of Computational and Graphical Statistics},
  pages={1--12},
  year={2022},
  publisher={Taylor \& Francis}
}

```

```
@article{L2E-Liu,
  title={A Sharper Computational Tool for L$_2$E  Regression},
  author={Liu, Xiaoqian and Chi, Eric C. and Lange, Kenneth},
  journal={arXiv preprint arXiv:2203.02993},
  year={2022}
}
```