library(igraph)
library(NetCoMi)

set.seed(10)

# -----------------------
# 1) Create networks
# -----------------------
n <- 20
p <- 0.15
m <- 2

g_er <- sample_gnp(n = n, p = p, directed = FALSE, loops = FALSE)
g_ba <- sample_pa(n = n, m = m, directed = FALSE)

g_er <- simplify(as_undirected(g_er), remove.multiple = TRUE, remove.loops = TRUE)
g_ba <- simplify(as_undirected(g_ba), remove.multiple = TRUE, remove.loops = TRUE)



# -----------------------
# 2) Adjacency matrices (binary, diag=0)
# -----------------------
A1 <- as.matrix(igraph::as_adjacency_matrix(g_er, sparse = FALSE))
A2 <- as.matrix(igraph::as_adjacency_matrix(g_ba, sparse = FALSE))

A1[A1 != 0] <- 1; diag(A1) <- 0
A2[A2 != 0] <- 1; diag(A2) <- 0

# -----------------------
# 3) NetCoMi GCD object
# -----------------------
gcd <- calcGCD(A1, A2)
cat("NetCoMi GCD =", gcd$GCD, "\n\n")

# -----------------------
# 4) GDV (orbit-count vectors) from NetCoMi
# IMPORTANT: use the FULL matrices, do NOT slice [1:20,]
# -----------------------
gdv1 <- gcd$ocount1
gdv2 <- gcd$ocount2

dim(gdv1)
dim(gdv2)

head(gdv1)
head(gdv2)

# -----------------------
# 5) GCM = correlation of GDV columns (NetCoMi definition)
# -----------------------
gcm1_manual <- cor(log(gdv1 + 1), method = "spearman", use = "pairwise.complete.obs")
gcm2_manual <- cor(log(gdv2 + 1), method = "spearman", use = "pairwise.complete.obs")

gcm1_manual[is.na(gcm1_manual)] <- 0; diag(gcm1_manual) <- 1
gcm2_manual[is.na(gcm2_manual)] <- 0; diag(gcm2_manual) <- 1

round(gcm1_manual, 2)
round(gcm2_manual, 2)

# -----------------------
# 6) Verify you did it right (compare to NetCoMi’s stored GCMs)
# -----------------------
max(abs(gcm1_manual - gcd$gcm1))
max(abs(gcm2_manual - gcd$gcm2))

# If both are ~0 (or extremely tiny like 1e-12), you matched NetCoMi exactly.

gcd <- calcGCD(A1, A2)
gcd
# Orbit counts per node (Graphlet Degree Vectors)
oc1 <- gcd$ocount1[1:20,]
oc2 <- gcd$ocount2[1:20,]

dim(oc1)
dim(oc2)

# Show first 6 nodes (rows) and all orbit-columns used by NetCoMi
head(oc1)
head(oc2)

gcm1 <- gcd$gcm1
gcm2 <- gcd$gcm2

round(gcm1, 2)
round(gcm2, 2)


gcmtest <- testGCM(gcd)

gcm1_heat <- plotHeat(gcmtest$gcm1, pmat = gcmtest$pAdjust1, type = "mixed")
gcm2_heat <- plotHeat(gcmtest$gcm2, pmat = gcmtest$pAdjust2, type = "mixed")
gcm1_manuel_heat <- plotHeat(round(gcm1_manual, 2), pmat = gcmtest$pAdjust1, type = "mixed")
gcm2_manuel_heat <- plotHeat(round(gcm2_manual, 2), pmat = gcmtest$pAdjust2, type = "mixed")




par(mfrow = c(1, 1))

plot(
  g_er,
  vertex.size  = 15,
  vertex.label = V(g_er)$label,
  vertex.label.cex = 0.7,
  vertex.label.color = "black",
  edge.arrow.mode = 0,
  main = "Erdős–Rényi (nodes numbered)"
)

plot(
  g_ba,
  vertex.size  = 15,
  vertex.label = V(g_ba)$label,
  vertex.label.cex = 0.7,
  vertex.label.color = "black",
  edge.arrow.mode = 0,
  main = "Barabási–Albert (nodes numbered)"
)

gdv1_df <- as.data.frame(gdv1)
gdv1_df$Node <- 1:21


table1 <- knitr::kable(
  gdv1_df[1:7,1:11],
  caption = "GDV with explicit node indices"
)
gdv2_df <- as.data.frame(gdv2)
gdv2_df$Node <- 1:21


table2 <- knitr::kable(
  gdv2_df[1:7,c(12,1:11)],
  caption = "GDV with explicit node indices"
)

