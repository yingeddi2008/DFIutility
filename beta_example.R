# example PCoA with bray-curtis and Unifrac

library(phyloseq)
library(gridExtra)

# PCoA --------------------------------------------------------------------

phy <- readRDS("<your phyloseq file name>")

# filter out sample with less than 200 total reads, you can change this cut off to 100 if necessary
phy <- subset_samples(phy, sample_sums(phy) > 200 )
phy <- subset_samples(phy, taxa_sums(phy) > 0 )

# normalization + rarefy
phynor <- phyloseq::transform_sample_counts(phy, function(x) x / sum(x)*2000 )
phyrr <- phyloseq::rarefy_even_depth(phy)

# bray-curtis
bcpcoa <- ordinate(phynor, distance = "bray", method = "PCoA")
bplt <- plot_ordination(phynor,bcpcoa, color = "group") +
  labs(title = "PCoA with Bray-Curtis distance")

# unifrac
unipcoa <- ordinate(phyrr, distance = "unifrac", weighted = F, method = "PCoA")
uplt <- plot_ordination(phyrr, unipcoa, color = "group") +
  labs(title = "PCoA with UniFrac distance")

# put them side by side
gridExtra::grid.arrange(bplt, uplt, ncol = 2)
