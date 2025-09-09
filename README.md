# EmuCP
Enhancing approximate modular Bayesian inference by emulating the conditional posterior
https://www.sciencedirect.com/science/article/pii/S0167947325001112

This R package implements the Emulating the Conditional Posterior (ECP) algorithm. This algorithm improves on the multiple imputation algorithm for sampling from a cut posterior distribution by using an emulator to learn a mapping from cut parameters to a parametric approximation of the conditional posterior. The emulator serves to increase the number of imputations for the approximation to the cut distribution with minimal computation.

The package includes two toy examples showing how to use the code.

