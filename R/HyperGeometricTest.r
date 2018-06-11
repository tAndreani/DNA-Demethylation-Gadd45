#Population size (All the Genes): 28409
#Sample size (Genes): 131
#Number of items in the  Population size that are classified as successes:
#Number of items in the sample that are classified as successes 

#Rep1 Up and Down
totalpop <- 37993 
sample1 <- 131    #(Genes Diff Expressed)
sample2 <- 2115    #(Genes with HyperMe)
sum(dhyper(14:sample2, sample1, totalpop-sample1, sample2))

#Rep1 only Down
totalpop <- 37993 
sample1 <- 91    #(Genes Down regulated Diff Expressed)
sample2 <- 2115    #(Genes with HyperMe)
sum(dhyper(10:sample2, sample1, totalpop-sample1, sample2))

#Rep2 Up and Down
totalpop <- 37993 
sample1 <- 131    #(Genes Diff Expressed)
sample2 <- 2120    #(Genes with HyperMe)
sum(dhyper(10:sample2, sample1, totalpop-sample1, sample2))

#Rep2 only Down
totalpop <- 37993 
sample1 <- 91    #(Genes Diff Expressed)
sample2 <- 2120    #(Genes with HyperMe)
sum(dhyper(9:sample2, sample1, totalpop-sample1, sample2))


#Comm Rep1 and Rep2 Down
totalpop <- 37993 
sample1 <- 91    #(Genes Diff Expressed)
sample2 <- 865    #(Genes with HyperMe)
sum(dhyper(6:sample2, sample1, totalpop-sample1, sample2))


#Comm Rep1 and Rep2 Up and Down
totalpop <- 37993 
sample1 <- 131    #(Genes Diff Expressed)
sample2 <- 865    #(Genes with HyperMe)
sum(dhyper(7:sample2, sample1, totalpop-sample1, sample2))


#Rep1 and Rep2 Up and Down
totalpop <- 37993 
sample1 <- 131    #(Genes Diff Expressed)
sample2 <- 3370     #(Genes with HyperMe)
sum(dhyper(17:sample2, sample1, totalpop-sample1, sample2))


#Genes with HyperDMRs obtained with GREAT Up and Down
totalpop <- 37993 
sample1 <- 131    #(Genes Diff Expressed)
sample2 <- 7389     #(Genes with HyperMe)
sum(dhyper(34:sample2, sample1, totalpop-sample1, sample2))


#Previous approach 1kb
totalpop <- 37993 
sample1 <- 130    #(Genes Diff Expressed)
sample2 <- 510     #(Genes with HyperMe)
sum(dhyper(6:sample2, sample1, totalpop-sample1, sample2))

#Extension 10 Kb Genes Diff Up and Down
totalpop <- 37993 
sample1 <- 130    #(Genes Diff Expressed)
sample2 <- 3651     #(Genes with HyperMe)
sum(dhyper(22:sample2, sample1, totalpop-sample1, sample2))

#Extension 10 Kb Genes Diff Down
totalpop <- 37993 
sample1 <- 91    #(Genes Diff Expressed)
sample2 <- 3651     #(Genes with HyperMe)
sum(dhyper(16:sample2, sample1, totalpop-sample1, sample2))


#Extension 20 Kb Genes Diff Up nad Down
totalpop <- 37993 
sample1 <- 130    #(Genes Diff Expressed)
sample2 <-  6065    #(Genes with HyperMe)
sum(dhyper(30:sample2, sample1, totalpop-sample1, sample2))

#Extension 20 Kb Genes Diff Down
totalpop <- 37993 
sample1 <- 91    #(Genes Diff Expressed)
sample2 <- 6065     #(Genes with HyperMe)
sum(dhyper(22:sample2, sample1, totalpop-sample1, sample2))
