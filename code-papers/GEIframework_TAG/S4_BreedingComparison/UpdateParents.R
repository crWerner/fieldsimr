# Script removes 10 oldest parents each year and replaces them 
# with most recent genotypes from the EYT stage to form a new crossing block
cat(" Update parents \n")

# Drop 10 parents
Parents <- Parents[11:nParents]

# Update with new 10 parents from the EYT stage
Parents <- c(Parents, EYT)

# Report parameters
gv_tpe <- multi_asr_output(pop = Parents, ntraits = 1, nenvs = nEnvTPE, nreps = nRepEYT, cov.mat = covs_tpe)
output <- save_output(output, Parents, "Parents", year, gv_tpe = gv_tpe)
