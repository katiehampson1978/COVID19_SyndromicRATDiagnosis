# Model fitting config values
# number of chains
N_Chains <- 2
# restrict the initial range to limit numerical overflow 
# (and thus inits being rejected)
Init_R <- 0.1

## Improve mixing
# Think of like step length in Gibbs, this will slow things down
AdaptDelta <- 0.99
# Only effects performance if reached
TreeDepth <- 100