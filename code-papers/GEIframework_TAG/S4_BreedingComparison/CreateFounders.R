# Create founders
cat(" Simulating founder parents \n")

# Generate initial haplotypes
founderPop <- runMacs(nInd     = nParents, 
                      nChr     = 10, 
                      segSites = nQtl,
                      inbred   = TRUE, 
                      species  = "WHEAT")
SP = SimParam$new(founderPop)

# Create an input object with FieldSimR for AlphaSimR
# to obtain desired genotype slopes 
input_asr <- multi_asr_input(nenvs  = nEnvTPE, 
                             mean   = MeanG,
                             var    = diag(De), 
                             corA   = Ce, 
                             nterms = k)

# Add additive trait
SP$addTraitA(nQtlPerChr = nQtl,
             mean       = input_asr$mean,
             var        = input_asr$var,
             corA       = input_asr$corA)

# Create founder parents
Parents <- newPop(founderPop)

#-------------
rm(founderPop)

