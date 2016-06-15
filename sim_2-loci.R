max <- 1000
N <- 1000000
sims<- 1000
loss1 <- 0.095
loss4 <- 0.095


#Bottleneck is in percentage
bottleneck <- 0.99990

# this is a new argument that says by what factor your 
# species can increase its population in a single generation
# setting it to 4 would mean that if the current population 
# were 50 that the next generation could increase only to 200
# in the following generation this could increase to 800 this would
# continue till the value N set above is reached
max.pop.inc <- 100

source("fnx1_M.R")

#result for all simulations
major.result1 <- data.frame(matrix(, sims, 1000))
major.result4 <- data.frame(matrix(, sims, 1000))



for(i in 1:sims){
  ## SETUP FOR SIMULATION
  # we will assume that initially are population isn't in a bottleneck
  Ne <- N
  pop <- vector(length=16, mode="numeric")
  pop <- c(1000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  names(pop) <- c("ABAB", "ABAb", "ABaB", "ABab", "AbAB", "AbAb", "AbaB", "Abab", 
                  "aBAB", "aBAb", "aBaB", "aBab", "abAB", "abAb", "abaB", "abab")

  # gamete table
  gametes <- vector(length=9, mode="numeric")
  names(gametes) <- c("AB", "Ab", "aB", "ab", "Ab'", "a'B", "ab'", "a'b", "a'b'")
  
  # a matrix to store results of a single simulation
  results <- data.frame(matrix(, 16, max))
  row.names(results) <- names(pop)
  
  # increment variable
  counter <- 1
  
  # record starting conditions
  results[, counter] <- pop
  
  ## Simulation
  while(counter < max){
    counter <- counter + 1
    pop[1:16] <- selection(pop, 
                           N = Ne)
    gametes[1:9] <- makeGametes(pop, 
                                r = 0.5, 
                                N = Ne, 
                                U.a = loss1, 
                                U.b = loss4
                                )
    
    # this is the bottleneck occuring every 7 and 10 generations
    if(counter %% 10 == 7 | counter %% 10 == 0 ){
      Ne <- N - (bottleneck * N)
    }else{
      if(Ne < N) Ne <- Ne * max.pop.inc
      if(Ne > N) Ne <- N
    }
    # perform the bottleneck
    pop <- makePop(gametes, N = Ne)
    results[, counter] <- pop
  }
  # storing the results
  major.result1[i, 1:1000] <- mon1.fnx(results=results, max=max)
  major.result4[i, 1:1000] <- mon4.fnx(results=results, max=max)

}


#Plotting frequencies of M1 and M4
par(mfrow=c(2,2))
matplot(t(major.result1), type="l", col = rainbow(10), lty = 1, ylim = range(0:1), xlab="Generations", ylab="Freq. of M1")
matplot(t(major.result4), type="l", col = rainbow(10), lty = 1, ylim = range(0:1), xlab="Generations", ylab="Freq. of M4")

#col="#ff000010"
#col = "#0000ff10"
#col = rainbow(10)

# checking to make sure trackAllele works
plot(trackAllele(results, "B"), type="l", ylim = range(0:1), col="blue", xlab ="generations", ylab = "Frequency of M4")
plot(trackAllele(results, "b"))
plot(trackAllele(results, "A"), type="l", ylim = range(0:1), col="red", xlab ="generations", ylab = "Frequency of M1")
plot(trackAllele(results, "a"))

# check to make sure checkAllele works
checkAllele(results, "A")
checkAllele(results, "a")
checkAllele(results, "B")
checkAllele(results, "b")

# using single simulation result of interest
plot(x=1:1000, y=major.result1[262,], type="l", ylim = range(0:1), col="red", xlab ="generations", ylab = "Frequency of M1 & M4")
lines(x=1:1000, y=major.result4[262,], type="l", ylim = range(0:1), col="blue", xlab ="generations", ylab = "Frequency ofM1 & M4")
points(x=c(850,850),y=c(0.6,0.5),pch=17, col=c("red","blue"))
text(x=c(850,850), y=c(0.6,0.5), labels=c("M1","M4"),pos=2)

