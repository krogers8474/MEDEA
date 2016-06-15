# MEDEA
R code written to run a stochastic, finite population model with MEDEA incorporated.

First step is to initialize all the paramters on the simulation (sim_2-loci.R) page. Parameters are:
max = generations
N = population size
sims = simulations to run
loss1 = mutation rate against M1 factor
loss4 = mutation rate against M4 factor
bottleneck = precentage of population removed 
max.pop.inc = factor population increases per generation after a bottleneck occurs.

Also need to set up the initial population vecotr for all 16 genotypes.
  If population is fixed for M1 and M4 with 1000 individuals total:
  pop <- (1000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
Then name the genotypes for each element:
  names(pop) <- c("ABAB", "ABAb", "ABaB", "ABab", "AbAB", "AbAb", "AbaB", "Abab", 
                  "aBAB", "aBAb", "aBaB", "aBab", "abAB", "abAb", "abaB", "abab")
                  
The simulation page calls the functions from fnx1_M.R which has all the functions listed to run through the generations.

