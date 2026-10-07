
#Respiratory and blood mcicrobiome analysis


#Load packages

library(ggplot2)
library(phyloseq)
library(maaslin3)
library(dplyr)
library(tidyr)
library(patchwork)
library(microViz)
library(ggtext)
library(stringr)
library(readr)
library(tibble)
library(ggsignif)
library(lme4)
library(emmeans)
library(ggrepel)


# A. PRE ANALYSIS

 #Step 1: Create phyloseq object

#read files
otu_po    <- read.csv("input_data/otu_table.csv", row.names = 1, check.names = FALSE)
tax_po    <- read.csv("input_data/tax_table.csv", row.names = 1, check.names = FALSE)
sample_po <- read.csv("input_data/sample_data.csv", row.names = 1, check.names = FALSE)


tax_po$Phyla <- gsub("Firmicutes_.*", "Firmicutes", tax_po$Phyla)

#Convert to matrix/dataframe structures
otu_mat   <- as.matrix(otu_po)
tax_mat   <- as.matrix(tax_po)
sample_df <- as.data.frame(sample_po)

#build phyloseq object
phylobject <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(sample_df)
)



#Step 2: turn low abundant taxa to abundance of 0 

# 1. Transform to relative abundance (0 to 1 scale)
phy_rel <- transform_sample_counts(phylobject, function(x) x / sum(x))

# 2. Create a logical mask: TRUE for counts >= 0.1% (0.001) within each sample
mask <- otu_table(phy_rel) >= 0.001

# 3. Extract original absolute counts as a matrix
otu_counts <- as(otu_table(phylobject), "matrix")

# 4. Set any taxon count with < 0.1% relative abundance in that sample to 0
otu_counts[!mask] <- 0

# 5. Create the new filtered phyloseq object
phylobject_filtered <- phylobject
otu_table(phylobject_filtered) <- otu_table(otu_counts, taxa_are_rows = TRUE)


#Step 3: remove Cutibacterium. acnes, Staphylococcus epidermidis, staphylococcus saprophyticus from dataset as they might be contaminants.

phylobject_filtered <- subset_taxa(phylobject_filtered, 
                                   !(Species %in% c("Cutibacterium acnes", 
                                                    "Staphylococcus epidermidis", 
                                                    "Staphylococcus saprophyticus")))


##################################
# Section B:

#CREATE FIGURES AND PERFORM STATISTICAL ANALYSIS

# Assign colours to phyla

Phyla_colors <- c( "Firmicutes"  = "#A6CEE3",
                   "Proteobacteria" = "#B2DF8A",
                   "Bacteroidota"= "#E7298A",
                   "Actinobacteriota" =  "#7570B3",
                   "Fusobacteriota" = "#666666",
                   "Patescibacteria" = "#83d5af")

# Assign colours to species

species_colors <- c(
  "Streptococcus pneumoniae" = "#A6CEE3",
  "Haemophilus influenzae/aegyptius"= "#1F78B4",
  "Moraxella catarrhalis" = "#B2DF8A",
  "Moraxella nonliquefaciens"  = "#E31A1C",
  "Escherichia coli" = "#FDBF6F",
  "Pseudomonas PSLC41"  = "#FF7F00",
  "Klebsiella pneumoniae" = "#CAB2D6",
  "Streptococcus PSLC27" = "#6A3D9A",
  "Streptococcus PSLC49" = "#FFFF99",
  "Neisseria meningitidis" = "#B15928",
  "Suttonella indologenes" = "#1ff8ff",
  "Staphylococcus aureus" = "#1B9E77",
  "Neisseria lactamica"  = "#D95F02",
  "Corynebacterium accolens" = "#7570B3",
  "Dolosigranulum pigrum" = "#E7298A",
  "Veillonella nakazawae" = "#66A61E",
  "Prevotella melaninogenica" = "#E6AB02",
  "Corynebacterium pseudodiphtheriticum"  = "#A6761D",
  "Ornithobacterium hominis" = "#666666",
  "Actinobacillus ureae"  = "#4b6a53",
  "Haemophilus parahaemolyticus"  = "#b249d5",
  "Prevotella histicola" = "#7edc45",
  "Dolosigranulum PSLC40"  = "#5c47b8",
  "Streptococcus agalactiae" = "#cfd251",
  "Pseudomonas fluorescens"  = "#ff69b4",
  "Moraxella lincolnii" = "#69c86c",
  "Neisseria polysaccharea"  = "#cd3e50",
  "Limosilactobacillus PSLC31"  = "#83d5af",
  "Other"  = "lightgrey")

#########

##Figure 1, sample availability


# 1. Load data
Sample_availability <- read.csv("input_data/Sample_availability.csv")


# 2. Define order for specimen columns
desired_order <- c("NPS", "NP.OP", "IS", "Blood", "LA")

# 3. set the corrected column name
long_data_stratified <- Sample_availability %>%
  mutate(profile = paste(NPS, NP.OP, IS, Blood, LA, sep = "_")) %>%
  arrange(profile) %>% 
  pivot_longer(cols = all_of(desired_order), names_to = "Sample_Type", values_to = "Available") %>%
  mutate(
    Sample_Type = recode(Sample_Type, "NP.OP" = "NP-OP"),
    Sample_Type = factor(Sample_Type, levels = c("NPS", "NP-OP", "IS", "Blood", "LA")),
    
    Participant = factor(Participant, levels = rev(unique(Participant)))
  )

# 4. Create plot
figure_1 <- ggplot(long_data_stratified, aes(x = Sample_Type, y = Participant, fill = Available)) +
  geom_tile(color = "white", linewidth = 0.3) + 
  
  scale_fill_manual(values = c("Yes" = "#2c7bb6", "No" = "#e0e0e0")) + 
  
  theme_minimal(base_size = 14) +
  labs(x = "Specimen", y = "Participant", fill = "Available") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12), 
    axis.text.y = element_text(size = 14),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 17)
  )

ggsave("Figure 1_sample availability.png", 
      plot = figure_1,
      width = 9,     
      height = 7,   
      dpi = 300)     

###############################################


######## Figure 2, MAGs and AMR

#Figure 2A

class_colors <- c(
  "Sulfonamide" = "#0075DC" ,
  "Aminoglycoside" = "#993F00",
  "Macrolide" =  "#4C005C" ,
  "Beta-lactam" = "#005C31",
  "Tetracycline" = "#be0032",
  "MLS" = "#FFA8BB",
  "Fluoroquinolone" = "#FFE100", 
  "Peptide antibiotic" = "#00998F" 
)


amr_prevalence_table <- read.csv("amr_prevalence_table.csv")

#Panel a: resistance gene prevalence

f2a <- ggplot(amr_prevalence_table,
              aes(x = Prevalence_Percent, y = best_hit_aro, fill = drug_class)) +
  geom_col() +
  theme_minimal() +
  scale_fill_manual(values = class_colors, guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(
    x    = "Prevalence (%) in high-quality MAGs",
    y    = "Antimicrobial resistance gene",
    fill = "Class"
  ) +
  theme(
    axis.text.y    = element_markdown(size = 10, face = "italic"),
    axis.text.x    = element_text(size = 10),
    axis.title     = element_text(size = 12, face = "bold"),
    legend.title   = element_text(size = 11),
    legend.position = "bottom",
    legend.text    = element_text(size = 10)
  ) +
  guides(fill = guide_legend(ncol = 3))


# Panel b: per-MAG AMR gene presence/absence heatmap

amr_data <- read.csv ("input_data/amr_data.csv")

amr_data <- amr_data %>% 
  mutate(
    specimen = factor(
      specimen,
      levels = c("NP-OP", "NPS", "IS", "LA", "Blood")
    )
  )

f2b <- ggplot(
  amr_data,
  aes(x = mag, y = gene, fill = presence)
) +
  geom_tile(
    color = "white",
    linewidth = 0.4
  ) +
  ggh4x::facet_wrap2(
    ~ specimen,
    scales = "free",
    ncol = 3,
    axes = "all",
    strip.position = "bottom"
  ) +
  scale_fill_manual(
    values = c(
      "Present" = "darkred",
      "Absent" = "beige"
    )
  ) +
  scale_x_discrete(position = "bottom") +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x.bottom = element_text(
      angle = 45,
      hjust = 1,
      face = "bold.italic"
    ),
    axis.ticks.x.bottom = element_line(),
    axis.text.y = element_text(
      face = "bold.italic",
      size = 11
    ),
    strip.placement = "outside",
    strip.text = element_text(
      face = "bold",
      size = 14,
      vjust = 0
    ),
    panel.grid = element_blank(),
    panel.spacing = unit(2, "lines"),
    legend.position = "bottom"
  )


# --- Combine

figure_2 <- f2a + f2b +
  plot_layout(widths = c(3, 7)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ')')


ggsave("figure_2_MAGs_AMR.png",
       plot = figure_2, width = 16, height = 8, dpi = 300)

#######################################

###Figure 3, composition barplot by specimen type

Phylobject_figure3 <- phylobject_filtered

sample_data(Phylobject_figure3)$specimen_type <- factor(
  sample_data(Phylobject_figure3)$specimen_type,
  levels = c("NP-OP", "NPS", "IS", "Blood", "LA"),
  labels = c("NP-OP (n=12)", "NPS (n=29)", "IS (n=12)", "Blood (n=11)", "LA (n=7)")
)
desired_order <- c("LA (n=7)", "Blood (n=11)", "IS (n=12)", "NP-OP (n=12)", "NPS (n=29)")

#species plot

species_figure3 <- Phylobject_figure3 %>%
  ps_select(specimen_type , case_control_sub_classification) %>% 
  phyloseq::merge_samples(group = "specimen_type") %>%
  comp_barplot(tax_level = "Species", palette = species_colors, 
            n_taxa = 14, bar_width = 0.8, sample_order = desired_order) +
  coord_flip() + 
  labs(y = "Relative abundance") +
  guides(fill = guide_legend(title = NULL, ncol = 2), 
         legend.key.size = unit(0.3, "cm")) +
  theme(axis.text.y = element_blank(), legend.position = "bottom", 
        legend.text = element_text(face = "italic"))

#phylum plot

phyla_figure3 <- Phylobject_figure3 %>%
  ps_select(specimen_type, case_control_sub_classification) %>% 
  phyloseq::merge_samples(group = "specimen_type") %>%
  comp_barplot(tax_level = "Phyla", palette = Phyla_colors, 
               bar_width = 0.8, sample_order = desired_order) +
  coord_flip() + 
  labs(y = "Relative abundance", x = "Specimen") +
  guides(fill = guide_legend(title = NULL, ncol = 2), 
         legend.key.size = unit(0.3, "cm")) +
  theme(legend.position = "bottom")

#combine species and phylum plots

figure_3 <- (phyla_figure3 + theme(plot.tag.position = c(0.14, 1.0), plot.tag = element_text(face = "bold", size = 14))) + 
  (species_figure3 + theme(plot.tag.position = c(0.05, 1.0), plot.tag = element_text(face = "bold", size = 14))) + 
  plot_annotation(tag_levels = 'a')


ggsave("Figure 3_compositional barplot.png", 
       plot = figure_3,
       width = 10,     
       height = 7,    
       dpi = 300)     


#######################################################

#Figure 4, Shannon diversity and statistical analysis

# Building analysis dataframe 
alpha_df <- phylobject_filtered %>%
  ps_calc_diversity(rank = "Species", index = "shannon") %>%
  samdat_tbl() %>%
  mutate(specimen_type = factor(specimen_type,
                                levels = c("NP-OP", "IS", "NPS", "Blood", "LA")))


# checking Sample and participant counts per specimen type (for figure annotation) 
n_summary <- alpha_df %>%
  group_by(specimen_type) %>%
  summarise(n_samples = n(),
            n_participants = n_distinct(patid), .groups = "drop")

#Mixed-effects model:
mixed_model <- lmer(shannon_Species ~ specimen_type + (1 | patid), data = alpha_df)


# Specimen Pairwise comparisons 
emm <- emmeans(mixed_model, ~ specimen_type)
pairwise_result <- pairs(emm, adjust = "tukey")
pairwise_df <- as.data.frame(pairwise_result)

pairwise_df <- pairwise_df %>%
  mutate(
    group1 = str_trim(str_split_fixed(contrast, " - ", 2)[, 1]),
    group2 = str_trim(str_split_fixed(contrast, " - ", 2)[, 2]),
    group1 = str_remove_all(group1, "[()]"),
    group2 = str_remove_all(group2, "[()]"),
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(stars))

sig_comparisons <- purrr::map2(pairwise_df$group1, pairwise_df$group2, c)
sig_annotations <- pairwise_df$stars

# Create Figure with sample sizes on x-axis 
axis_labels <- n_summary %>%
  mutate(label = paste0(specimen_type, "\n(n=", n_samples, ")"))

figure_4 <- alpha_df %>%
  ggplot(aes(x = specimen_type, y = shannon_Species, fill = specimen_type)) +
  geom_boxplot(width = 0.3, color = "black") +
  theme_bw() +
  scale_fill_manual(values = c(
    "IS"    = "#E41A1C",
    "NP-OP" = "#377EB8",
    "NPS"   = "#4DAF4A",
    "Blood" = "#984EA3",
    "LA"    = "#FF7F00"
  )) +
  scale_x_discrete(labels = setNames(axis_labels$label, axis_labels$specimen_type)) +
  labs(y = "Shannon Diversity (Species)", x = "Specimen") +
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  geom_signif(
    comparisons  = sig_comparisons,
    annotations  = sig_annotations,
    step_increase = 0.08,
    tip_length    = 0.01,
    textsize      = 5
  )


ggsave("figure_4_shannon_diversity.png",
       plot = figure_4, width = 12, height = 7, dpi = 300)


################


# Figure 5, Jaccard index with permutation-test significance 

# Calculate species-level Jaccard similarity
jaccard <- phylobject_filtered %>%
  tax_agg(rank = "Species") %>%
  tax_transform("binary") %>%
  dist_calc(dist = "jaccard") %>%
  dist_get()

jaccard_index <- 1 - as.matrix(jaccard)


#Prepare sample metadata

sample_lookup <- samdat_tbl(phylobject_filtered) %>%
  select(.sample_name, patid, specimen_type)


# descriptive data used for Figure 5

# 'combined_jaccard_index_data'- This file contains the selected pairwise comparisons to calculate the descriptive Jaccard means shown
# in Figure 5.

combined_jaccard_index_data <- read.csv(
  "combined_jaccard_index_data.csv"
)

# Calculate mean, SD and SE for each specimen-pair comparison
jaccard_summary <- combined_jaccard_index_data %>%
  group_by(pair) %>%
  summarise(
    mean_jaccard = mean(jaccard_index),
    sd_jaccard = sd(jaccard_index),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    se_jaccard = sd_jaccard / sqrt(n)
  )


# Keep only blood-respiratory comparisons
jaccard_summary_blood <- jaccard_summary %>%
  filter(grepl("Blood", pair)) %>%
  mutate(
    pair_label = gsub(" vs ", " and ", pair),
    mean_jaccard_percent = mean_jaccard * 100,
    se_jaccard_percent = se_jaccard * 100
  )

# Permutation-test function

permutation_test_jaccard <- function(
    type_a,
    type_b,
    sample_lookup,
    jaccard_matrix,
    n_perm = 9999,
    seed = 123) {
  
  set.seed(seed)
  
# Identify samples belonging to each specimen type
  ids_a <- sample_lookup %>%
    filter(specimen_type == type_a)
  
  ids_b <- sample_lookup %>%
    filter(specimen_type == type_b)
  
  
# Identify participants who have both specimen types
  paired_patids <- intersect(
    ids_a$patid,
    ids_b$patid
  )
  
# Identify the corresponding samples for each participant
  samples_a <- ids_a$.sample_name[
    match(paired_patids, ids_a$patid)
  ]
  
  samples_b <- ids_b$.sample_name[
    match(paired_patids, ids_b$patid)
  ]
  
  
# Observed mean Jaccard similarity
  
# Calculate Jaccard similarity for the true
# participant-matched pairs
  observed_values <- mapply(
    function(a, b) {
      jaccard_matrix[a, b]
    },
    samples_a,
    samples_b
  )
  
  observed_mean <- mean(observed_values)
  

# Generate null distribution

# Keep specimen type A fixed and randomly shuffle
# specimen type B across participants.
  null_means <- replicate(
    n_perm,
    {
      
      shuffled_b <- sample(samples_b)
      
      shuffled_values <- mapply(
        function(a, b) {
          jaccard_matrix[a, b]
        },
        samples_a,
        shuffled_b
      )
      
      mean(shuffled_values)
    }
  )
  
 
  # Tests whether the true participant-matched pairing produces greater similarity than random pairing.
  p_value <- mean(
    null_means >= observed_mean
  )
  
  list(
    pair = paste(type_a, "vs", type_b),
    n_pairs = length(paired_patids),
    observed_mean = observed_mean,
    null_mean = mean(null_means),
    null_sd = sd(null_means),
    p_value = p_value
  )
}


# Run permutation tests

# NPS versus blood
result_nps_blood <- permutation_test_jaccard(
  type_a = "NPS",
  type_b = "Blood",
  sample_lookup = sample_lookup,
  jaccard_matrix = jaccard_index,
  n_perm = 9999,
  seed = 123
)

# Induced sputum versus blood
result_is_blood <- permutation_test_jaccard(
  type_a = "IS",
  type_b = "Blood",
  sample_lookup = sample_lookup,
  jaccard_matrix = jaccard_index,
  n_perm = 9999,
  seed = 123
)

# NP/OP versus blood
result_npop_blood <- permutation_test_jaccard(
  type_a = "NP-OP",
  type_b = "Blood",
  sample_lookup = sample_lookup,
  jaccard_matrix = jaccard_index,
  n_perm = 9999,
  seed = 123
)

# Combine permutation-test results

results_summary <- bind_rows(
  lapply(
    list(
      result_nps_blood,
      result_is_blood,
      result_npop_blood
    ),
    function(x) {
      
      data.frame(
        pair = x$pair,
        n_pairs = x$n_pairs,
        observed_mean = x$observed_mean,
        null_mean = x$null_mean,
        null_sd = x$null_sd,
        p_value = x$p_value
      )
    }
  )
)


# Adjust the three permutation-test P values using the Benjamini-Hochberg procedure
results_summary$q_value <- p.adjust(
  results_summary$p_value,
  method = "BH"
)

results_summary$significant <- results_summary$q_value < 0.05

# Add permutation-test results to Figure 5 data

jaccard_summary_blood <- jaccard_summary_blood %>%
  left_join(
    results_summary %>%
      select(
        pair,
        p_value,
        q_value,
        significant
      ),
    by = "pair"
  ) %>%
  mutate(
    
    # Significance label based on BH-adjusted q value
    sig_label = case_when(
      q_value < 0.001 ~ "***",
      q_value < 0.01  ~ "**",
      q_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    
    # Add sample size to x-axis labels
    pair_label = paste0(
      pair_label,
      "\n(n=",
      n,
      ")"
    )
  )


#  Figure 5


figure_5 <- ggplot(
  jaccard_summary_blood,
  aes(
    x = reorder(
      pair_label,
      -mean_jaccard_percent
    ),
    y = mean_jaccard_percent
  )
) +
  
  # Mean Jaccard similarity
  geom_bar(
    stat = "identity",
    fill = "lightblue",
    width = 0.3
  ) +
  
  # Standard error
  geom_errorbar(
    aes(
      ymin = mean_jaccard_percent -
        se_jaccard_percent,
      ymax = mean_jaccard_percent +
        se_jaccard_percent
    ),
    width = 0.1,
    color = "black"
  ) +
  
  # Statistical significance
  geom_text(
    aes(
      label = sig_label,
      y = mean_jaccard_percent +
        se_jaccard_percent +
        4
    ),
    size = 7
  ) +
  
  # Y-axis
  scale_y_continuous(
    name = "Mean Jaccard similarity (%)",
    breaks = c(
      0,
      20,
      40,
      60,
      80,
      100
    ),
    limits = c(0, 100)
  ) +
  
  # X-axis
  labs(
    x = "Blood and respiratory specimen pair"
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    axis.title = element_text(
      size = 16
    ),
    axis.text = element_text(
      size = 14
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )


ggsave(
  "figure_5_jaccard.png",
  plot = figure_5,
  width = 10,
  height = 6,
  dpi = 300
)


############################################################


###Figure_6, NPS composition by severity

f6 <- phylobject_filtered %>%
  ps_mutate(case_control_sub_classification = recode(case_control_sub_classification,
                                                     "severe case" = "Severe cases (n=15)",
                                                     "very severe case" = "Very severe cases (n=14)")) %>%
  ps_select(specimen_type, case_control_sub_classification) %>%
  ps_filter(specimen_type == "NPS") %>%
  phyloseq::merge_samples(group = "case_control_sub_classification") %>%
  comp_barplot(tax_level = "Species", n_taxa = 14, bar_width = 0.8, other_name = "Other") +
  coord_flip() +
  labs(x = "Disease severity", y = "Relative abundance")


figure_6 <- f6 +
  scale_fill_manual(
    values = species_colors,
    labels = function(x) parse(text = paste0("italic('", x, "')"))
  ) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(face = "italic", size = 10),
    legend.key.size = unit(0.8, "lines")
  )


ggsave("figure_6_NPS_severity.png",
       plot = figure_6, width = 10, height = 6, dpi = 300)


##############

#Comparing S. pneumoniae relative abundance and Hi relative abundance by vaccination status.

# --- S. pneumoniae abundance vs PCV vaccination status, NPS only ---
sp_abundance <- phylobject_filtered %>%
  ps_filter(specimen_type == "NPS") %>%
  tax_transform("compositional", rank = "Species") %>%
  ps_get() %>%
  psmelt() %>%
  filter(Species == "Streptococcus pneumoniae") %>%
  mutate(vaccination_collapsed = case_when(
    pcv_vaccination_status == "full" ~ "full",
    pcv_vaccination_status %in% c("none", "partial") ~ "incomplete",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(vaccination_collapsed))


sp_wilcox <- wilcox.test(Abundance ~ vaccination_collapsed, data = sp_abundance)

#  H. influenzae/aegyptius abundance vs Hib vaccination status, NPS only

hi_abundance <- phylobject_filtered %>%
  ps_filter(specimen_type == "NPS") %>%
  tax_transform("compositional", rank = "Species") %>%
  ps_get() %>%
  psmelt() %>%
  filter(Species == "Haemophilus influenzae/aegyptius") %>%
  mutate(vaccination_collapsed = case_when(
    hib_vaccination_status == "full" ~ "full",
    hib_vaccination_status %in% c("none", "partial") ~ "incomplete",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(vaccination_collapsed))


hi_wilcox <- wilcox.test(Abundance ~ vaccination_collapsed, data = hi_abundance)


# --- BH correction across the two tests ---

sp_wilcox_p <- sp_wilcox$p.value
hi_wilcox_p <- hi_wilcox$p.value

p_adjust_result <- p.adjust(c(sp_wilcox_p, hi_wilcox_p), method = "BH")
names(p_adjust_result) <- c("S. pneumoniae vs PCV", "H. influenzae vs Hib")


sp_abundance %>%
  group_by(vaccination_collapsed) %>%
  summarise(median_abundance = median(Abundance), mean_abundance = mean(Abundance), n = n())

hi_abundance %>%
  group_by(vaccination_collapsed) %>%
  summarise(median_abundance = median(Abundance), mean_abundance = mean(Abundance), n = n())


#####################################################################

#Supplementary 

##Figure s1 barplots by participant

##order specimens
specimen_levels <- c("NPS", "NP-OP", "IS", "Blood", "LA")

#function to build individual plots

make_pair <- function(spec_type) {
  
  ps_sub <- phylobject_filtered %>%
    ps_filter(specimen_type == spec_type)
  
  samp_df <- data.frame(sample_data(ps_sub)) %>%
    rownames_to_column("sample_id") %>%
    mutate(participant_num = as.numeric(participant_label)) %>%
    arrange(participant_num)
  
  ordered_samples <- samp_df$sample_id
  
  phyla_plot <- ps_sub %>%
    comp_barplot("Phyla", label = "participant_label", sample_order = ordered_samples) +
    coord_flip() +
    labs(y = "Relative abundance", x = "Participant") +
    scale_fill_manual(values = Phyla_colors) +
    guides(fill = "none") +
    theme(axis.ticks.y = element_blank(),
          legend.position = "none",
          axis.text.x  = element_text(size = 11),   
          axis.title.x = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(face = "bold"))
  
  species_plot <- ps_sub %>%
    tax_fix() %>%
    comp_barplot("Species", label = "participant_label", n_taxa = 14,
                 sample_order = ordered_samples) +
    coord_flip() +
    labs(y = "Relative abundance", x = NULL) +
    scale_fill_manual(values = species_colors,
                      labels = function(x) parse(text = paste0("italic('", x, "')"))) +
    guides(fill = "none") +
    theme(axis.ticks.y = element_blank(),
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 11))  
  
  phyla_label   <- wrap_elements(full = grid::textGrob("Phylum",   gp = grid::gpar(fontsize = 14)))
  species_label <- wrap_elements(full = grid::textGrob("Species", gp = grid::gpar( fontsize = 14)))
  
  phyla_col   <- phyla_label   / phyla_plot   + plot_layout(heights = c(0.05, 1))
  species_col <- species_label / species_plot + plot_layout(heights = c(0.05, 1))
  
  combo <- phyla_col | species_col
  
  title_grob <- wrap_elements(full = grid::textGrob(
    spec_type, gp = grid::gpar(fontface = "bold", fontsize = 16)
  ))
  
  title_grob / combo + plot_layout(heights = c(0.07, 1))
}

pairs_list <- lapply(specimen_levels, make_pair)
names(pairs_list) <- specimen_levels

# buid shared legends 
phyla_legend_plot <- phylobject_filtered %>%
  comp_barplot("Phyla", label = "participant_label") +
  scale_fill_manual(values = Phyla_colors, name = "Phyla") +
  theme(legend.text  = element_text(size = 12),
        legend.title = element_text(face = "bold"))

species_legend_plot <- phylobject_filtered %>%
  tax_fix() %>%
  comp_barplot("Species", label = "participant_label", n_taxa = 14) +
  scale_fill_manual(values = species_colors, name = "Species",
                    labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  theme(legend.text  = element_text(size = 12, face = "italic"),
        legend.title = element_text(face = "bold"))

get_leg <- function(p) {
  tryCatch(
    cowplot::get_legend(p),
    error = function(e) cowplot::get_plot_component(p, "guide-box", return_all = TRUE)[[1]]
  )
}

phyla_legend   <- get_leg(phyla_legend_plot)
species_legend <- get_leg(species_legend_plot)

legend_panel <- wrap_plots(
  wrap_elements(full = phyla_legend),
  wrap_elements(full = species_legend),
  ncol = 2
)

figure_s1 <- wrap_plots(
  pairs_list[["NPS"]], pairs_list[["NP-OP"]], pairs_list[["IS"]],
  pairs_list[["Blood"]], pairs_list[["LA"]], legend_panel,
  ncol = 3, nrow = 2
)


ggsave("figure_s1_individual_barplots.png", 
       plot = figure_s1, width = 18, height = 14, 
       dpi = 300)

#########################


# Figure S2, pie charts of dominant species by specimen type

dominant_species <- read.csv ("input_data/dominant_species.csv")

dominant_summary <- dominant_species %>%
  mutate(
    Specimen_Type = case_match(
      Specimen_Type,
      "NP-OP" ~ "NP-OP (n=12)",
      "IS"    ~ "IS (n=12)",
      "NPS"   ~ "NPS (n=29)",
      "Blood" ~ "Blood (n=11)",
      "LA"    ~ "LA (n=7)"
    )
  ) %>%
  group_by(Specimen_Type, Species) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  ) %>%
  group_by(Specimen_Type) %>%
  mutate(
    Percentage = round((Count / sum(Count)) * 100, 1)
  ) %>%
  arrange(
    Specimen_Type,
    desc(Percentage)
  ) %>%
  mutate(
    ymax = cumsum(Percentage),
    ymin = lag(ymax, default = 0),
    ymid = (ymax + ymin) / 2
  ) %>%
  ungroup() %>%
  mutate(
    Specimen_Type = factor(
      Specimen_Type,
      levels = c(
        "NPS (n=29)",
        "NP-OP (n=12)",
        "IS (n=12)",
        "Blood (n=11)",
        "LA (n=7)"
      )
    ),
    repel = Specimen_Type == "NPS (n=29)" &
      Percentage == 3.4,
    x_text = ifelse(repel, NA, 0.6),
    x_label = ifelse(repel, 1.1, NA)
  )

figure_s2 <- ggplot(dominant_summary,
                    aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = Species)) +
  geom_rect(color = "black") +
  coord_polar(theta = "y") +
  facet_wrap(~ Specimen_Type, ncol = 3) +
  scale_fill_manual(values = species_colors) +
  geom_label_repel(
    data = subset(dominant_summary, repel),
    aes(x = x_label, y = ymid, label = paste0(Percentage, "%")),
    fill = "white", size = 3.5, box.padding = 0.3,
    direction = "y", segment.color = "grey40", show.legend = FALSE
  ) +
  geom_text(
    data = subset(dominant_summary, !repel),
    aes(x = x_text, y = ymid, label = paste0(Percentage, "%")),
    size = 4, fontface = "bold", color = "black"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.text = element_text(face = "italic"),
    strip.text  = element_text(face = "bold", size = 12)
  ) +
  labs(fill = "Dominant Species")

figure_s2

ggsave("figure_s2_dominant_species_pies.png",
       plot = figure_s2, width = 12, height = 10, dpi = 300)

######################
# Figure S3, alpha diversity in NPS by case severity and statistical analysis

phylobject_filtered_NPS <- phylobject_filtered %>%
  ps_filter(specimen_type == "NPS")

richness_NPS <- estimate_richness(phylobject_filtered_NPS)

wilcox.test(richness_NPS$Shannon ~ 
              sample_data(phylobject_filtered_NPS)$case_control_sub_classification)


figure_s3 <- phylobject_filtered_NPS %>% 
  ps_mutate(case_control_sub_classification = recode(case_control_sub_classification,
                                                     "severe case" = "Severe cases (n=15)",
                                                     "very severe case" = "Very severe cases (n=14)")) %>%
  ps_calc_diversity(rank = "Species", index = "shannon") %>%
  samdat_tbl() %>%
  ggplot(aes(y = case_control_sub_classification, x = shannon_Species,
             fill = case_control_sub_classification)) +
  geom_boxplot(width = 0.3, color = "black") +
  theme_bw() +
  scale_fill_manual(values = c("Severe cases (n=15)" = "#377EB8", "Very severe cases (n=14)" = "#E41A1C")) +
  coord_flip() +
  labs(y = "Case status", x = "Shannon Diversity (Species)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14, face = "bold"), # Increases "Case status" and "Shannon Diversity"
    axis.text = element_text(size = 13)                 
)


ggsave("figure_s3_NPS_diversity_by_severity.png",
       plot = figure_s3, width = 10, height = 6, dpi = 300)



#NPS and shannon diversity multivariable model
nps_df <- alpha_df %>% filter(specimen_type == "NPS")

# Multivariable model with collapsed vaccination strategy

nps_df <- nps_df %>%
  mutate(vaccination_status_collapsed = case_when(
    pcv_vaccination_status == "full" ~ "full",
    pcv_vaccination_status %in% c("none", "partial") ~ "incomplete",
    TRUE ~ NA_character_
  ))

nps_severity_model_v2 <- lm(shannon_Species ~ very_severe_pneumonia + age_months + vaccination_status_collapsed,
                            data = nps_df)



###############################

# Figure S4: beta diversity (PCoA) and PERMANOVA, NPS by case severity

figure_s4 <- phylobject_filtered %>%
  ps_filter(specimen_type == "NPS") %>%
  ps_mutate(case_control_sub_classification = recode(case_control_sub_classification,
                                                     "severe case" = "Severe cases (n=15)",
                                                     "very severe case" = "Very severe cases (n=14)")) %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 3, color = "case_control_sub_classification") +
  theme_classic(12) +
  coord_fixed(0.7) +
  stat_ellipse(aes(color = case_control_sub_classification), level = 0.95) +
  scale_color_manual(
    values = c("Severe cases (n=15)" = "#377EB8", "Very severe cases (n=14)" = "#E41A1C"),
    name = NULL
  ) +
  guides(color = guide_legend(title = NULL)) +
  labs(caption = NULL) +
  theme (legend.text = element_text(size = 12))


ggsave("figure_s4_NPS_PCoA.png",
       plot = figure_s4, width = 8, height = 7, dpi = 300)

# PERMANOVA: severity + age + Vaccination 

new_sample_data <- phylobject_filtered %>%
  samdat_tbl() %>%
  mutate(vaccination_status_collapsed = case_when(
    pcv_vaccination_status == "full" ~ "full",
    pcv_vaccination_status %in% c("none", "partial") ~ "incomplete",
    TRUE ~ NA_character_
  ))

sample_data(phylobject_filtered) <- new_sample_data %>%
  column_to_rownames(".sample_name") %>%
  sample_data()

permanova_result <- phylobject_filtered %>%
  ps_filter(specimen_type == "NPS") %>%
  ps_filter(!is.na(vaccination_status_collapsed)) %>%   # match complete-case set used in the lm
  dist_calc(dist = "aitchison") %>%
  dist_permanova(variables = "case_control_sub_classification + age_months + vaccination_status_collapsed",
                 n_perms = 999, seed = 123) %>%
  perm_get()

packageVersion("microViz")


####################################################################

#Maaslin3 analysis, differences in severity using NPS samples only


abundance_data <- otu_table(phylobject_filtered_NPS)
abundance_data_df <- as.data.frame(abundance_data)

# Get the species names from the tax_table
specie_names <- as.character(tax_table(phylobject_filtered_NPS)[,"Species"])

# Set the column names to the specie names
abundance_data_df <- t(abundance_data_df)

colnames(abundance_data_df) <- specie_names

# Extract metadata from phyloseq object
metadata_df <- as(sample_data(phylobject_filtered_NPS), "data.frame")


output <- maaslin3(input_data = abundance_data_df,
                   input_metadata = metadata_df,
                   output = "output/maaslin3_abundance",
                   formula = '~ case_control_sub_classification + age_months',
                   normalization = 'CLR',
                   transform = 'NONE',
                   augment = TRUE,
                   standardize = TRUE,
                   min_abundance = 0,        
                   min_prevalence = 0.15,    
                   max_significance = 0.1,
                   median_comparison_abundance = TRUE,
                   median_comparison_prevalence = FALSE,
                   warn_prevalence = FALSE,
                   evaluate_only = 'abundance',   
                   max_pngs = 250,
                   cores = 1)


writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)

#*********************************END***************************************
