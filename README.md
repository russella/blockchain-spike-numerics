# blockchain-spike-numerics
Stability event probability estimates for PoW and PoS blockchains in a setting with adversarial spikes

This repository contains two directories of interest "pow" and "pos"; they contain independent codebases. To use either of them, navigate to the directory and "make all". The results will be two executables with names poX and poXthr (where X = "a" or "w").

- The poX executable computes the resulting settlement errors (as a function of time or slots) for the appropriate setting with a particular, given, adversarial budget.
- The poXthr executable computes the number of timesteps (or slots) necessary to achieve a particular error value.

Remarks.
- In the pos case, the program takes the resulting Binomial distribution parameter and the spike budget.
- In the pow case, honest and adversarial rates are described as Poisson parameters (giving the expected number of successes per unit time). A step-to-step approximation guarantee is requested which determines how long the distribution is evolved prior to the double spend attack--it stops when two back-to-back distributions are within this error bound in total deviation.
