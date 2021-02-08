library(data.table)
library(igraph)
library(randnet)
library(multiviewtest)
library(doParallel)

minDegree <- 1

binary.RAW <- fread("./data/HomoSapiens_binary_hq.txt")
cocomp.RAW <- fread("./data/HomoSapiens_cocomp_hq.txt")

binary.RAW[, id_A := paste(Uniprot_A, Gene_A, sep="-")]
binary.RAW[, id_B := paste(Uniprot_B, Gene_B, sep="-")]
cocomp.RAW[, id_A := paste(Uniprot_A, Gene_A, sep="-")]
cocomp.RAW[, id_B := paste(Uniprot_B, Gene_B, sep="-")]

# Remove self-loops
binary <- binary.RAW[id_A != id_B, ]
cocomp <- cocomp.RAW[id_A != id_B, ]

binary.gr <- graph_from_data_frame(binary[, .(id_A, id_B)], directed=F)
cocomp.gr <- graph_from_data_frame(cocomp[, .(id_A, id_B)], directed=F)

# Filter to nodes which appear in both graphs 
binary.ids <- V(binary.gr)$name
cocomp.ids <- V(cocomp.gr)$name
union.ids <- union(binary.ids, cocomp.ids)
intersect.ids <- intersect(binary.ids, cocomp.ids)

protein.attr <- data.frame(protein=union.ids, stringsAsFactors=F) 
protein.attr <- merge(protein.attr, 
                      data.frame(protein=binary.ids, 
                                 binary.degree=degree(binary.gr), 
                                 stringsAsFactors=F), 
                      all.x=T, sort=F)
protein.attr$binary.degree[is.na(protein.attr$binary.degree)] <- 0
protein.attr <- merge(protein.attr, 
                      data.frame(protein=cocomp.ids, 
                                 cocomp.degree=degree(cocomp.gr), 
                                 stringsAsFactors=F), 
                      all.x=T, sort=F)
protein.attr$cocomp.degree[is.na(protein.attr$cocomp.degree)] <- 0

protein.attr$filter <- (protein.attr$binary.degree >= minDegree) & (protein.attr$cocomp.degree >= minDegree)
new.binary.gr <- induced_subgraph(binary.gr, 
                                  vids=protein.attr$protein[protein.attr$filter])
new.cocomp.gr <- induced_subgraph(cocomp.gr, vids=protein.attr[protein.attr$filter, ]$protein)

# Pick number of communities in each view, with max = 100 
adj.bin <- as.matrix(as_adjacency_matrix(new.binary.gr))*1 
rownames(adj.bin) <- NULL 
colnames(adj.bin) <- NULL 
adj.cocomp <- as.matrix(as_adjacency_matrix(new.cocomp.gr))*1
rownames(adj.cocomp) <- NULL 
colnames(adj.cocomp) <- NULL 
hintdat <- list(adj.bin, adj.cocomp)

K1 <- BHMC.estimate(adj.bin, K.max=100)
K2 <- BHMC.estimate(adj.cocomp, K.max=100)

# Run P2LRT
RNGkind(sample.kind = "Rounding")
set.seed(123, kind="L'Ecuyer-CMRG")
registerDoParallel(cores=12)
p2lrt.res <- test_indep_com(hintdat, K1=K1, K2=K2, nperm=10^4, step=1e-5, parallel=TRUE)

save(p2lrt.res, file="final-results-spectral.Rdata")
