# filter ktu phyloseq
library(phyloseq)
library(tidyverse)
library(vegan)

its <- readRDS("data/intermediate/psITS_ktu_filt.rds")
its2 <- readRDS("data/intermediate/ITSx/psITS2_ktu_filt.rds")

#load functions for filtering phyloseq objects
dir("code/functions", "ps_filter", full.names = T) %>% sapply(source)
source("code/functions/otu_table2df.r")
source("code/functions/ps_merge.r")

# get sequences from phyloseq
its_seqs <-
  its %>%
  refseq() %>%
  as.character

its2_seqs <-
  its2 %>%
  refseq() %>%
  as.character

its2_seqs_noPh <-
  its2 %>%
  ps_filter_phylum_is_NA() %>%
  refseq() %>%
  as.character

# histogram read length--------
df_hist <- data.frame(
  nt = c(its2_seqs %>% nchar, its_seqs %>% nchar),
  dataset = c(rep("its2", length(its2_seqs)), rep("its", length(its_seqs)))
  )
                      
# procrustes analysis to compare relative versus absolute counts
library(tidyverse)
library(phyloseq)
library(vegan)

# 16s
# absolute counts
ps16s_abstemp <-
  readRDS("data/intermediate/ps16s_abs.rds")
# geometric mean of absolute counts
geomM_bact <-
  sample_sums(ps16s_abstemp) %>%
  {exp(mean(log(.)))}
ps16s_abs <-
  ps16s_abstemp %>%
  transform_sample_counts(function(x) log(1 + x))

# relative
ps16s_rel <-
  readRDS("data/intermediate/ps16s_ktu_filt_rel.rds") %>%
  transform_sample_counts(function(x) x / sum(x)) %>% # proportions
  transform_sample_counts(function(x) x * geomM_bact) %>%
  transform_sample_counts(function(x) log(1 + x))

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
protest_16s <-
  protest(X = bray_16s_abs$vectors[,1:2],
          Y = bray_16s_rel$vectors[,1:2],
          scores = "sites", permutations = 999)

# plot histogram with distances
bray_abs <- phyloseq::distance(ps16s_abs, method ="bray")
bray_rel <- phyloseq::distance(ps16s_rel, method ="bray")
dist16s <- bray_abs - bray_rel
pdf(file = "output/hist_bray_dist_16s.pdf",
    width = 5, height = 4)
hist(dist16s, main = "bray abs - bray rel")
dev.off()

# its
# absolute counts
psits_abstemp <-
  readRDS("data/intermediate/psITS_abs.rds")
# geometric mean of absolute counts
geomM_its <<-
  sample_sums(psits_abstemp) %>%
  {exp(mean(log(.)))}
psits_abs <-
  psits_abstemp %>%
  transform_sample_counts(function(x) log(1 + x))

# relative quantification
psits_rel <-
  readRDS("data/intermediate/psITS_ktu_filt_rel.rds") %>%
  transform_sample_counts(function(x) x / sum(x)) %>% # proportions
  transform_sample_counts(function(x) x * geomM_its) %>%
  transform_sample_counts(function(x) log(1 + x))

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
protest_its <-
  protest(X = bray_its_abs$vectors[,1:2],
          Y = bray_its_rel$vectors[,1:2],
          scores = "sites", permutations = 999)

# plot histogram with distances
bray_abs <- phyloseq::distance(psits_abs, method ="bray")
bray_rel <- phyloseq::distance(psits_rel, method ="bray")
distits <- bray_abs - bray_rel
pdf(file = "output/hist_bray_dist_its.pdf", width = 5, height = 4)
hist(distits, main = "bray abs - bray rel")
dev.off()

# print protest results to text file
sink("output/protest.txt")
cat("\nProtest 16s\n")
protest_16s
cat("\n\nProtest ITS\n\n")
protest_its
sink()

# print dist results to file
sink("output/dist_beta_div.txt")
cat(
  paste("\n16S data abs-rel dist. Average",
        mean(dist16s),
        "standard deviation",
        sd(dist16s))
)
cat(
  paste("\nITS data abs-rel dist. Average",
        mean(distits),
        "standard deviation",
        sd(distits))
)
sink()
<-
  ggplot(df_hist, aes(x = nt)) + 
  geom_histogram(data = subset(df_hist, dataset == "its"), fill = "grey30", alpha = 0.5) +
  geom_histogram(data = subset(df_hist, dataset == "its2"), fill = "grey70", alpha = 0.5) +
  xlab("nucleotides") +
  ylab("frequency") +
  theme_classic()

# assigned taxonomies-------

temp_df1 <-
  its %>% tax_table %>% as.data.frame %>%
  apply(2, function(x) is.na(x) %>% sum) %>%
  data.frame(nas = .) %>% rownames_to_column("tax") %>%
  mutate(total = ntaxa(its),
         no_nas = total - nas,
         dataset = "ITS_full") %>%
  select(-nas)

temp_df2 <-
  its2 %>% tax_table %>% as.data.frame %>%
  apply(2, function(x) is.na(x) %>% sum) %>%
  data.frame(nas = .) %>% rownames_to_column("tax") %>%
  mutate(total = ntaxa(its2),
         no_nas = total - nas,
         dataset = "ITS2") %>%
  select(-nas)

p_barplot <-
  rbind(temp_df1, temp_df2) %>%
    pivot_longer(
      cols = c(total, no_nas), values_to = "nas") %>%
    mutate(tax = factor(tax, levels = rank_names(its))) %>%
    filter(name == "no_nas") %>%
    ggplot(aes(x = tax, y = nas, fill = dataset)) +
    geom_bar(stat = "identity", position=position_dodge()) +
    scale_fill_manual(values = c("darkseagreen", "coral1")) +
    xlab(NULL) +
    ylab("Number of taxa") +
    theme_minimal()

# procrustes ---------------------
# MDS with weighted UNIFRAC distances
bray_its <-
  phyloseq::ordinate(its,
                     method = "MDS",
                     distance = "bray")

bray_its2 <-
  phyloseq::ordinate(its2,
                     method = "MDS",
                     distance = "bray")

pro1 <- vegan::procrustes(bray_its$vectors[,1:2],
                          bray_its2$vectors[,1:2],
                          scores = "sites")


# save figures-----------
ggsave(filename = "output/itsx_hist_itsVSits2.pdf",
       plot = phist,
       width = 5,
       height = 4)
ggsave(filename = "output/itsx_tax_assignations.pdf",
       plot = p_barplot,
       width = 6,
       height = 5)

pdf(file = "output/itsx_procrustes_16s.pdf",
    width = 7,
    height = 5)
plot(pro1, kind = 1, type = "p")
dev.off()
