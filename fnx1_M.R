# this function takes the base population and returns a new population table 
# that has the number of contributions each genotype makes to the gamete pool
# with this implementation we actually currently dont do anything with the 
# selection fucntion.  However, I can imagine exploring how much of a fitness
# must be present to stop the spread of Medea so lets keep this framework handy
selection <- function(pop, Ne) {
  # create an empty vector for absolute fitness
  abs.fit <- vector(mode = "numeric", length = 16)
  #If wanted to add selection against MEDEA
  # here we could decrement abs.fit by a value to represent the cost of
  # being medea positive.
  #abs.fit[1:16] <- 1 - s * c(4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0)
  abs.fit[1:16] <- 1
  
  # convert absolute fitness to relative or marginal fitness
  w <- abs.fit / max(abs.fit[pop != 0]) * (pop / Ne)
  
  # use marginal fitness to draw from population
  as.vector(rmultinom(1:16, size = Ne, prob = w))
}

# this fnx takes a population table and generates the gamete pool
makeGametes <- function(pop, r, Ne, U.a, U.b){
  # gamete pool: AB,  Ab, aB, ab, Ab', a'B, ab', a'b, a'b'
  # in these cases a prime symbol indicates "poisoned" for that 
  # allele.
  #
  #create an empty vector for gametes
  x <- vector(length=9, mode="numeric")
  # AB 
  x[1] <- (sum(pop[1], pop[2]*(1/2), pop[3]*(1/2), pop[4]*(1/2)*(1-r), pop[5]*(1/2), 
               pop[7]*r*(1/2), pop[9]*(1/2), pop[10]*r*(1/2), pop[13]*(1/2)*(1-r)))/Ne
  # Ab 
  x[2] <- (sum(pop[2]*(1/4), pop[4]*r*(1/4), pop[5]*(1/4), pop[6], pop[7]*(1/4)*(1-r),
               pop[8]*(1/2), pop[10]*(1/4)*(1-r), pop[13]*r*(1/4), pop[14]*(1/2)))/Ne
  # aB 
  x[3] <- (sum(pop[3]*(1/4), pop[4]*r*(1/4), pop[7]*(1/4)*(1-r), pop[9]*(1/4), pop[10]*(1/4)*(1-r),
               pop[11], pop[12]*(1/2), pop[13]*r*(1/4), pop[15]*(1/2)))/Ne
  # ab 
  x[4] <- (sum(pop[4]*(1/4)*(1-r), pop[7]*r*(1/4), pop[8]*(1/4), pop[10]*r*(1/4), pop[12]*(1/4),
               pop[13]*(1/4)*(1-r), pop[14]*(1/4), pop[15]*(1/4), pop[16]))/Ne

  # Add the 'poisened' gametes 
  #Ab'
  x[5] <- (sum(pop[2]*(1/4), pop[4]*r*(1/4), pop[5]*(1/4), pop[7]*(1/4)*(1-r), pop[10]*(1/4)*(1-r),
               pop[13]*r*(1/4))) / Ne 

  # a'B 
  x[6] <- (sum(pop[3]*(1/4), pop[4]*r*(1/4), pop[7]*(1/4)*(1-r), pop[9]*(1/4), pop[10]*(1/4)*(1-r),
                 pop[13]*r*(1/4))) / Ne
  
  # ab'
  x[7] <- (sum(pop[12]*(1/4), pop[15]*(1/4)))/ Ne

  # a'b 
  x[8] <- (sum(pop[8]*(1/4), pop[14]*(1/4)))/ Ne
  
  # a'b' 
  x[9] <- (sum(pop[4]*(1/4)*(1-r), pop[7]*r*(1/4), pop[10]*r*(1/4), pop[13]*(1/4)*(1-r)))/ Ne

  # Incorporation of mutation per gamete
  y <- vector(length=9, mode="numeric")
  #AB
  y[1] <- x[1] - U.a*x[1] - U.b*x[1] - U.a*U.b*x[1]
  #Ab
  y[2] <- x[2] - U.a*x[2] + U.b*x[1]
  #aB
  y[3] <- x[3] - U.b*x[3] + U.a*x[1]
  #ab
  y[4] <- x[4] + U.a*x[2] + U.b*x[3] + U.a*U.b*x[1]
  #Ab'
  y[5] <- x[5] - U.a*x[5]
  #a'B
  y[6] <- x[6] - U.b*x[6]
  #ab'
  y[7] <- x[7] + U.a*x[5]
  #a'b
  y[8] <- x[8] + U.b*x[8]
  #a'b' Nothing can mutate into this allele, as mutation into a poison allele is impossible
  y[9] <- x[9]
  ## now we return this vector
  return(y)
  
}

# this fnx takes a gamete pool and reconstitutes a
makePop <- function(gametes, Ne){
  # first we calculate probs based on frequency of gametes in pool
  # then we just draw from a multinomial distribution
  g <- gametes
  z.freq <- c(#1 ABxAB
              g[1] * g[1],
              #2 ABxAb + ABxAb'
              g[1] * g[2] + g[1] * g[5], 
              #3 ABxaB + ABxa'B
              g[1] * g[3] + g[1] * g[6],
              #4 ABxab + ABxab’ + ABxa’b + ABxa’b’
              g[1] * g[4] + g[1] * g[7] + g[1] * g[8] + g[1] * g[9],
              #5 AbxAB 
              g[2] * g[1], 
              #6 AbxAb
              g[2] * g[2], 
              #7 AbxaB + Abxa’B 
              g[2] * g[3] + g[2] * g[6], 
              #8 Abxab + Abxa’b
              g[2] * g[4] + g[2] * g[8],
              #9 aBxAB 
              g[3] * g[1], 
              #10 aBxAb + aBxAb’
              g[3] * g[2] + g[3] * g[5], 
              #11 aBxaB
              g[3] * g[3], 
              #12 aBxab + aBxab’
              g[3] * g[4] + g[3] * g[7],
              #13 abxAB 
              g[4] * g[1], 
              #14 abxAb 
              g[4] * g[2], 
              #15 abxaB 
              g[4] * g[3], 
              #16 abxab
              g[4] * g[4])
  z <- (rmultinom(1:16, Ne, prob = z.freq))
  
  return(z)
}


#Tracking allele frequency of M4 (B)
mon4.fnx <- function(results, max){
  count<-1
  x<-vector()
  for(j in seq.int(from=1, to=max)){
    x[count] <- ((2*sum(results[c(1,3,9,11), j]) + sum(results[c(2,4,5,7,10,12,13,15), j]))) / (2*sum(results[,j]))
    count<-count+1
  }
  return(x)
}

#Tracking allele frequency of M1 (A)
mon1.fnx <- function(results, max){
  count<-1
  y<-vector()
  for(s in seq.int(from=1, to=max)){
    y[count] <- ((2*sum(results[c(1,2,5,6), s]) + sum(results[c(3,4,7,8,9,10,13,14), s]))) / (2*sum(results[,s]))
    count<-count+1
  }
  return(y)
}


mon4.fnx_pop <- function(results, max){
  count<-1
  z<-vector()
  for(w in seq.int(from=1, to=max)){
    z[count] <- ((sum(results[c(1,2,3,5,7,9,10,11,12,13,15), w]))) 
    count<-count+1
  }
  return(z)
}

mon1.fnx_pop <- function(results, max){
  count<-1
  n<-vector()
  for(m in seq.int(from=1, to=max)){
    n[count] <- ((sum(results[c(1,2,3,4,5,6,7,8,9,10,13,14), m]))) 
    count<-count+1
  }
  return(n)
}


# # this function monitors the allele of interest
# trackAllele <- function(results, allele){
#   x<-vector()
#   for(j in 1:ncol(results)){
#     if(allele == "A" | allele == "a"){
#       x[j] <- ((2*sum(results[c(1,2,5,6), j]) + sum(results[c(3, 4, 7, 8, 9, 10, 13, 14), j]))) / (2*sum(results[,j]))
#     }
#     if(allele == "B" | allele == "b"){
#       x[j] <- ((2*sum(results[c(1,3,9,11), j]) + sum(results[c(2, 4, 5, 7, 10, 12, 13, 15), j]))) / (2*sum(results[,j]))
#     }
#   }
#   if(allele == "b" | allele == "a") x <- 1-x
#   return(x)
# }
# 
# # this fnx check to see if we A,a,B, or b are present at end of sim
# checkAllele <- function(results, allele){
#   # get final result
#   final <- results[, ncol(results)]
#   # get all alleles
#   x <- unlist(strsplit(as.character(row.names(results))[final != 0],split=""))
#   # test
#   return(allele %in% x)
# }
