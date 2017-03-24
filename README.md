# LZeroSpikeInference: A package for estimating spike times from calcium imaging data using an L0 penalty

This package implements an algorithm for deconvolving calcium imaging data for a single neuron in order to estimate the times at which the neuron spikes.

This algorithm solves the optimization problems
## AR(1) model
minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_{c_t neq gamma c_{t-1} }
for the global optimum, where y_t is the observed fluorescence at the tth
timepoint.

## AR(1) with intercept
minimize_{c1,...,cT,b1,...,bT} 0.5 sum_{t=1}^T (y_t - c_t - b_t)^2 + lambda sum_{t=2}^T 1_{c_t neq gamma c_{t-1}, b_t neq b_{t-1} }
where the indicator variable 1_{(A,B)} equals 1 if the event A cup B holds, and equals zero otherwise.

## Difference of Exponentials
minimize_{c1,...,cT,d1,...,dT} 0.5 sum_{t=1}^{T} ( y_t - (c_t - d_t) )^2 + lambda sum_{t =2}^{T} 1_{ c_t neq gamma _c c_{t-1}, d_t neq gamma_d d_{t-1} }.

See Jewell and Witten, Exact Spike Train Inference Via L0 Optimization (2017)

Install 
-----

If ``devtools`` is installed type 

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


