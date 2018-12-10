# LZeroSpikeInference: A package for estimating spike times from calcium imaging data using an L0 penalty [![Build Status](https://travis-ci.com/jewellsean/LZeroSpikeInference.svg?token=oixVftbrq2TkrSApRKn2&branch=master)](https://travis-ci.com/jewellsean/LZeroSpikeInference)

This package implements an algorithm for deconvolving calcium imaging data for a single neuron in order to estimate the times at which the neuron spikes. See [https://jewellsean.github.io/fast-spike-deconvolution/](https://jewellsean.github.io/fast-spike-deconvolution/) for tutorials and additional information. 

This algorithm solves the optimization problems
### AR(1) model
minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_{c_t neq gamma c_{t-1} }

for the global optimum, where y_t is the observed fluorescence at the tth timepoint. We also solve the above problem with the constraint that c_t >= 0 (hardThreshold = T). 

### AR(1) with intercept
minimize_{c1,...,cT,b1,...,bT} 0.5 sum_{t=1}^T (y_t - c_t - b_t)^2 + lambda sum_{t=2}^T 1_{c_t neq gamma c_{t-1}, b_t neq b_{t-1} }

where the indicator variable 1_{(A,B)} equals 1 if the event A cup B holds, and equals zero otherwise.

Install 
-----

In R, if ``devtools`` is installed type 

```r
devtools::install_github("jewellsean/LZeroSpikeInference")
```

Usage
----

Once installed type 
```{r}
library(LZeroSpikeInference)
?LZeroSpikeInference
```

Python
----

This package can be called from Python using the [py2](https://rpy2.bitbucket.io/) package. To install LZeroSpikeInference and rpy2 for use in Python first

1. Install R (for example `apt-get install r-base`)

and then from within R install this package (as above). Then pip install rpy2

2. pip install --user rpy2

The following example illustrates use of the LZeroSpikeInference package from python 

```python
from numpy import array
import rpy2.robjects.packages
lzsi = rpy2.robjects.packages.importr("LZeroSpikeInference")
d = lzsi.simulateAR1(n = 500, gam = 0.998, poisMean = 0.009, sd = 0.15, seed = 8)
fit = lzsi.estimateSpikes(d[1], **{'gam':0.998, 'lambda':8, 'type':"ar1"})
spikes = array(fit[0])
fittedValues = array(fit[1])
```

Thanks to Luke Campagnola for suggesting this approach!  

Reference
-----

See Jewell and Witten, [Exact Spike Train Inference Via L0 Optimization (2017)](https://arxiv.org/abs/1703.08644)
