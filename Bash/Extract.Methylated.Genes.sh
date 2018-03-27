#From http://www.geneimprint.com/site/genes-by-species.Mus+musculus copy the genes and select only the imprinted ones


grep Imprinted imprinted.and.not.imprinted.txt | grep Not -v > Imprinted.Genes.Mouse.txt
