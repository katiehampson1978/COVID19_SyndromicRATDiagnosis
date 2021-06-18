# Model fitting config values
# restrict the initial range to limit numerical overflow (inits being rejected)
Init_R <- 0.1
# iter includes warmup - default is half of iter will be warmup
Iter <- 2000
# improve mixing
# Think of like step length in Gibbs, this will slow things down
AdaptDelta <- 0.99
# Recommended by D Pascall, only effects performance if reached
TreeDepth <- 100