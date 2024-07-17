# Global parameters script contains parameters for simulation of
# a breeding program

#=======================================================================
##-- Simulation parameters 
#=======================================================================
nReps     = 1   # Number of simulation replicates
nCycles   = 20  # Number of simulation cycles

# Trait parameters
MeanG     = 0   # Overall trait mean
nQtl      = 300 # Number of QTLs per chromosome
k         = 7   # Number of environmental covariate terms

#=======================================================================
##-- Stage-specific parameters
#=======================================================================

nEnvTPE = 1000 # Target population of environments
nStages = 5

# Parents
nParents  = 30  # Number of parents (and founders)
nCrosses  = 100 # Number of crosses per year
nColParents  = 6
nRowParents  = 10

# HDRW
nDH       = 100 # Number of DH lines produced per cross
nEnvHDRW  = 1   # h2=0.07
nRepHDRW  = 1
nColHDRW  = 100
nRowHDRW  = 100
varE_HDRW = 8

# PYT
nPYT     = 500
famMax   = 10 # Maximum number of DH lines per cross
nEnvPYT  = 2  # h2=0.2
nRepPYT  = 1
nColPYT  = 25
nRowPYT  = 20
varE_PYT = 4 

# AYT
nAYT     = 50
nEnvAYT  = 5 # h2=0.66
nRepAYT  = 2
nColAYT  = 10
nRowAYT  = 10
varE_AYT = 4

# EYT
nEYT     = 10
nEnvEYT  = 20 # h2=0.90
nRepEYT  = 2
nColEYT  = 4
nRowEYT  = 5
varE_EYT = 2
