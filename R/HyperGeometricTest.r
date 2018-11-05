#Hyper DMRs in G45-TKO genes associated with 2C like genes
totalpop <-  24096  #Background 
sample1 <-  6817 #(Genes with HyperMe using GREAT standard paramenter)
sample2 <-  525   #(Genes 2C)
fisher.test(matrix(c(6817-178,24096-6817-525,178,525-178), 2, 2), alternative='l')

#Down in 2c G45 DKO and Up at 2C state
totalpop <- 25256  #Total genes background
sample1 <- 50    #(Genes Diff Expressed Down in 2C Gadd DKO)
sample2 <- 5265    #(Genes with 2C up)
fisher.test(matrix(c(5262-22,25256-5265-50,22,50-22), 2, 2), alternative='l')


#Up in 2c G45 DKO and Up at 2C state
totalpop <- 25256 #Total genes background
sample1 <- 54    #(Genes Diff Expressed Up in 2C Gadd DKO)
sample2 <- 5265    #(Genes with 2C up)
fisher.test(matrix(c(5262-3,25256-5265-54,3,54-3), 2, 2), alternative='l')


#Down in 2c G45 DKO and Down at 2C state
totalpop <- 25256 #Total genes background
sample1 <- 50    #(Genes Diff Expressed Down in 2C Gadd DKO)
sample2 <- 2325    #(Genes with 2C down)
fisher.test(matrix(c(2325-7,25256-2325-50,7,50-7), 2, 2), alternative='l')


#Up in 2c G45 DKO and Down at 2C state
totalpop <- 25256 #Total genes background
sample1 <- 54    #(Genes Diff Expressed Up in 2C Gadd DKO)
sample2 <- 2325    #(Genes with 2C down)
fisher.test(matrix(c(2325-17,25256-2325-54,17,54-17), 2, 2), alternative='l')

#Down in 2c G45 DKO and 2C-like genes (Fig. 6D)
totalpop <- 25256 #Total genes background
sample1 <- 50    #(Genes Diff Expressed Down in 2C Gadd DKO)
sample2 <- 525    #(2C like genes G45-DKO)
fisher.test(matrix(c(525-7,25256-525-50,7,50-7), 2, 2), alternative='l')


