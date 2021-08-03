# Gutmicro KITP project Summer 2021 day 7 tasks
# Now we are going to detect genetic differences between samples (Quasi-phasable samples)
# We calculate alt freq for each site for each population (fi1,fi2).
# Because of sampling noise that has a Binomial distribution around the true frequency
#   and because sampling noise is independent for each sample, we will get differences in the consensus 
#   that are false positives.
# We want to find sites where truely consensus 1 != consensus 2, but we can't expect full sweeps so
#   a threshold is chosen (80%)
# To distinguish between genotype 1 and 2, we check the freq for each site i. If fi > 80%, we assign it genotype 2
#   and when fi < 20 we assign it genotype 1 -> results for each bacteria in within_host_changes.txt 
