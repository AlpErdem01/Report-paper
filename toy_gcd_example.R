install.packages("Matrix")
library(igraph)
# ------------------------------------------------------------
# Toy networks for GCD seminar
# Step 1: create simple example graphs
# ------------------------------------------------------------

# Packages
library(igraph)
library(orca)
library(Matrix)
library(NetCoMi)
set.seed(1)

# -----------------------------
# Network 1: Erdos–Renyi (random)
# -----------------------------
n <- 20
p <- 0.15

g_er <- sample_gnp(
  n = n,
  p = p,
  directed = FALSE,
  loops = FALSE
)

g_er <- simplify(g_er)

# -----------------------------
# Network 2: Barabasi–Albert
# -----------------------------
m <- 2

g_ba <- sample_pa(
  n = n,
  m = m,
  directed = FALSE
)

g_ba <- simplify(g_ba)

# -----------------------------
# Sanity checks
# -----------------------------
cat("ER network:\n")
cat(" Nodes:", vcount(g_er), "\n")
cat(" Edges:", ecount(g_er), "\n")
cat(" Connected:", is.connected(g_er), "\n\n")

cat("BA network:\n")
cat(" Nodes:", vcount(g_ba), "\n")
cat(" Edges:", ecount(g_ba), "\n")
cat(" Connected:", is_connected(g_ba), "\n\n")

# -----------------------------
# Simple visual check
# -----------------------------
par(mfrow = c(1, 2))

plot(
  g_er,
  main = "Erdos–Renyi (Random)",
  vertex.size = 15,
  vertex.label = NA
)

plot(
  g_ba,
  main = "Barabasi–Albert (Preferential Attachment)",
  vertex.size = 15,
  vertex.label = NA
)

par(mfrow = c(1, 1))




as_undirected_safe <- function(g) {
  if (exists("as_undirected", where = asNamespace("igraph"), inherits = FALSE)) {
    igraph::as_undirected(g)
  } else {
    igraph::as.undirected(g)
  }
}

adjmat_safe <- function(g, sparse = TRUE) {
  # Prefer igraph's adjacency matrix converter, but provide fallbacks.
  if (exists("as_adjacency_matrix", where = asNamespace("igraph"), inherits = FALSE)) {
    igraph::as_adjacency_matrix(g, sparse = sparse)
  } else if (exists("get.adjacency", where = asNamespace("igraph"), inherits = FALSE)) {
    igraph::get.adjacency(g, sparse = sparse)
  } else {
    stop("No adjacency-matrix function found in your igraph. Update igraph.")
  }
}

prep_for_orca <- function(g) {
  g <- as_undirected_safe(g)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  # ORCA expects an unweighted 0/1 adjacency
  A <- adjmat_safe(g, sparse = TRUE)
  A <- Matrix::drop0(A)
  A[A != 0] <- 1
  list(g = g, A = A)
}

# ---------- GDV (Graphlet Degree Vectors) ----------
# For GCD-11 you typically use 2–4 node graphlet orbits (11 orbits).
# ORCA returns more orbits; we compute all available, then take first 11.
# ---------- GDV (GCD-11) via ORCA (works even when orca() is not exported) ----------

compute_gdv_gcd11 <- function(g) {
  
  # --- undirected + simple ---
  if (exists("as_undirected", where = asNamespace("igraph"), inherits = FALSE)) {
    g <- igraph::as_undirected(g)
  } else {
    g <- igraph::as.undirected(g)
  }
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  # --- IMPORTANT: drop isolates (ORCA cannot handle them) ---
  g <- igraph::delete_vertices(g, degree(g) == 0)
  
  # --- force clean 1..N indexing ---
  g <- igraph::induced_subgraph(g, V(g))
  igraph::V(g)$name <- as.character(seq_len(vcount(g)))
  
  # --- edge list (strictly positive, contiguous) ---
  E <- igraph::as_edgelist(g, names = FALSE)
  
  if (min(E) <= 0) {
    stop("Non-positive node ids detected after reindexing (should never happen).")
  }
  
  # --- ORCA GDVs ---
  gdv_all <- orca::count4(E)
  
  if (ncol(gdv_all) < 11) {
    stop("ORCA returned fewer than 11 orbits; cannot compute GCD-11.")
  }
  
  gdv_11 <- gdv_all[, 1:11, drop = FALSE]
  colnames(gdv_11) <- paste0("orbit_", 1:11)
  
  gdv_11
}




# ---------- GCM (Graphlet Correlation Matrix) ----------
# Standard choice: Spearman correlation of log(GDV+1)

compute_gcm <- function(gdv) {
  X <- log(gdv + 1)
  gcm <- suppressWarnings(cor(X, method = "spearman", use = "pairwise.complete.obs"))
  
  # Handle orbits with zero variance -> NA correlations
  gcm[is.na(gcm)] <- 0
  diag(gcm) <- 1
  gcm
}

# ---------- GCD (Graphlet Correlation Distance) ----------
# Euclidean distance between upper triangles (excluding diagonal)

gcd_from_gcms <- function(gcm1, gcm2) {
  if (!all(dim(gcm1) == dim(gcm2))) stop("GCM dimensions do not match.")
  idx <- upper.tri(gcm1, diag = FALSE)
  v1 <- gcm1[idx]
  v2 <- gcm2[idx]
  sqrt(sum((v1 - v2)^2))
}

# ---------- OPTIONAL: heatmap plotting ----------
plot_gcm_base <- function(gcm, main = "GCM") {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mar = c(6, 6, 3, 2))
  image(
    1:ncol(gcm), 1:nrow(gcm), t(gcm)[, nrow(gcm):1],
    axes = FALSE, xlab = "", ylab = "", main = main
  )
  axis(1, at = 1:ncol(gcm), labels = colnames(gcm), las = 2, cex.axis = 0.8)
  axis(2, at = 1:nrow(gcm), labels = rev(rownames(gcm)), las = 2, cex.axis = 0.8)
  box()
}

plot_gcm_gg <- function(gcm, title = "GCM") {
  df <- as.data.frame(as.table(gcm))
  colnames(df) <- c("Orbit1", "Orbit2", "rho")
  ggplot(df, aes(Orbit1, Orbit2, fill = rho)) +
    geom_tile() +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ============================================================
# USE YOUR EXISTING NETWORKS: g_er and g_ba (from your Step 1)
# ============================================================

# 1) GDV
gdv_er <- compute_gdv_gcd11(g_er)
gdv_ba <- compute_gdv_gcd11(g_ba)

# 2) GCM
gcm_er <- compute_gcm(gdv_er)
gcm_ba <- compute_gcm(gdv_ba)

# 3) GCD
gcd_er_ba <- gcd_from_gcms(gcm_er, gcm_ba)
cat("GCD (ER vs BA) =", gcd_er_ba, "\n")

# Optional quick checks / visuals
# plot_gcm_base(gcm_er, main = "GCM - Erdos-Renyi")
# plot_gcm_base(gcm_ba, main = "GCM - Barabasi-Albert")
# print(plot_gcm_gg(gcm_er, "GCM - Erdos-Renyi"))
# print(plot_gcm_gg(gcm_ba, "GCM - Barabasi-Albert"))
