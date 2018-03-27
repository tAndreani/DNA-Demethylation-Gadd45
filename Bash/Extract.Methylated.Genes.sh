#From http://www.geneimprint.com/site/genes-by-species.Mus+musculus copy the genes and select only the imprinted ones

grep Imprinted test2.save | grep Not -v > Imprinted.Genes.Mouse.txt
