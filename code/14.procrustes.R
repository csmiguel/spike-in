# mds plot from bray-curtis
library(tidyverse)
library(phyloseq)
library(vegan)

# 16s
# load phyloseq with absolute counts
ps16s_abs <-
  readRDS("data/intermediate/ps16s_abs.rds")
ps16s_rel <-
  readRDS("data/intermediate/ps16s_ktu_filt_rel.rds")

# MDS with weighted UNIFRAC distances
bray_16s_abs <-
  phyloseq::ordinate(ps16s_abs,
                     method = "MDS",
                     distance = "bray")

bray_16s_rel <-
  phyloseq::ordinate(ps16s_rel,
                     method = "MDS",
                     distance = "bray")

pro1 <- vegan::procrustes(bray_16s_abs$vectors[,1:2],
                          bray_16s_rel$vectors[,1:2],
                          scores = "sites")
pdf(file = "output/procrustes_16s.pdf",
    width = 7,
    height = 5)
plot(pro1, kind = 1, type = "p")
dev.off()

# test signficance
pro1_test <-
  protest(X = bray_16s_abs$vectors[,1:2],
        Y = bray_16s_rel$vectors[,1:2],
        scores = "sites", permutations = 999)

# `plot histogram with distances
bray_abs <- phyloseq::distance(ps16s_abs, method ="bray")
bray_rel <- phyloseq::distance(ps16s_rel, method ="bray")
dists <- bray_abs - bray_rel
pdf(file = "output/hist_bray_dist_16s.pdf",
    width = 5, height = 4)
hist(dists, main = "bray abs - bray rel")
dev.off()

# its
# load phyloseq with absolute counts
psits_abs <-
  readRDS("data/intermediate/psITS_abs.rds")
psits_rel <-
  readRDS("data/intermediate/psITS_ktu_filt_rel.rds")

# MDS with weighted UNIFRAC distances
bray_its_abs <-
  phyloseq::ordinate(psits_abs,
                     method = "MDS",
                     distance = "bray")

bray_its_rel <-
  phyloseq::ordinate(psits_rel,
                     method = "MDS",
                     distance = "bray")

pro1 <- vegan::procrustes(bray_its_abs$vectors[,1:2],
                          bray_its_rel$vectors[,1:2],
                          scores = "sites")
pdf(file = "output/procrustes_its.pdf",
    width = 7,
    height = 5)
plot(pro1, kind = 1, type = "p")
dev.off()

# test signficance
pro1_test <-
  protest(X = bray_its_abs$vectors[,1:2],
          Y = bray_its_rel$vectors[,1:2],
          scores = "sites", permutations = 999)

# `plot histogram with distances
bray_abs <- phyloseq::distance(psits_abs, method ="bray")
bray_rel <- phyloseq::distance(psits_rel, method ="bray")
dists <- bray_abs - bray_rel
pdf(file = "output/hist_bray_dist_its.pdf", width = 5, height = 4)
hist(dists, main = "bray abs - bray rel")
dev.off()
