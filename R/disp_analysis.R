# Load library
library(ggplot2)
library(cowplot) 
library(dplyr)
library(sp)
library(tidyverse)
library(ggeffects)
library(rotl)
library(ape)
library(brms)
library(rstan)
library(MCMCglmm)

mytheme <- function() {
  theme_bw() + 
    theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.8), # Border around plotting area.
          panel.grid.major = element_blank(), # Major grid lines blank
          panel.grid.minor = element_blank(), # Minor grid lines blank
          axis.line = element_blank(), # axis line size
          axis.ticks = element_line(colour = "black", size = 0.8),
          axis.text = element_text(size = 10, colour = "black"), # axis text size
          axis.title = element_text(size = 10), #axis title size
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel bg)
} # set up plot theme

# Set directory

# Load data
ind_raw_data <- read.csv("ind_disp_raw_data.csv")
pop_raw_data <- read.csv("pop_disp_raw_data.csv")
str(ind_raw_data)

# Cleaning individual raw data
ind_data_clean <- ind_raw_data %>% 
  tibble::rowid_to_column("es_ID") %>% 
  dplyr::mutate(Zr          = (1 / 2) * log((1 + corr_coeff) / (1 - corr_coeff)), # Fisher-transformed correlation coefficient
                Zr_v        = 1 / (sample_size - 3), # Sampling variance (v)
                Zr_sei      = 1 / sqrt(sample_size - 3), # Standard error (SE)
                Zr_inv      = 1 / Zr_sei, # Precision (inverse of SE) 
                Zr_z        = Zr / Zr_sei, # Egger - z score 
                Zr_w        = 1 / Zr_v, # Weight (inverse of Zr_v)
                year_centre = year - mean(year), # Mean-centring year of publication
                disp_trait  = factor(disp_trait, levels = c("Activity", "Exploration", "Dispersal")),
                lnMass_g    = log(mean_mass_g)
  )

# Cleaning population raw data
pop_data_clean <- pop_raw_data %>% 
  tibble::rowid_to_column("es_ID") %>% # add effect size id
  dplyr::mutate(n_core_sc    = n_core / shared_control,
                abslnRR      = abs(lnRR),
                v            = sd_front ^ 2 / (n_front * mean_front ^ 2) + sd_core ^ 2 / (n_core * mean_core ^ 2), # Sampling variance (v)
                sei          = sqrt(v), # Standard error (SE)
                sei_inv      = 1 / sei, # Precision (inverse of SE)
                eff_n        = (4 * n_core * n_front) / (n_core + n_front), # Effective sample size
                inv_eff_n    = (n_core + n_front) / (n_core * n_front), # Inverse of eff_n
                sqrt_inv_eff = sqrt(inv_eff_n), # Square root of the inverse of eff_n
                year_centre  = year - mean(year), # Mean-centring year of publication
                lnMass_g     = log(mean_mass_g),
                lnRate       = log(dispersal_rate_km_y),
                lnTime       = log(time_diff_y)
  )


## DATA SUMMARY ##--------------------------------------------------------------------------------
# Effect size summary
ind_data_clean %>% dplyr::summarise(count = n())
pop_data_clean %>% dplyr::summarise(count = n())

# Study numbers
ind_data_clean %>% dplyr::summarise(count = length(unique(study_ID)))
pop_data_clean %>% dplyr::summarise(count = length(unique(study_ID)))

# Species summary
ind_data_clean %>% dplyr::summarise(count = length(unique(species)))
pop_data_clean %>% dplyr::summarise(count = length(unique(species)))

# Species summary
ind_sum <- ind_data_clean %>% 
  dplyr::group_by(trait, response) %>% 
  dplyr::summarise(ef_n = n(),
                   study_n = length(unique(study_ID)),
                   species_n = length(unique(species)))

pop_sum <- pop_data_clean %>% 
  dplyr::group_by(trait, response) %>% 
  dplyr::summarise(ef_n = n(),
                   study_n = length(unique(study_ID)),
                   species_n = length(unique(species)))
  

## 1. INDIVIDUAL MOVEMENT ANALYSIS ##-------------------------------------------------------------------------------
# BUILD PHYLOGENY ALL #
ind_species <- sort(unique(as.character(ind_data_clean$species_OTL))) # generate list of species (as character format)
ind_taxa    <- rotl::tnrs_match_names(names = ind_species) # match taxonomic names to the OTL

# check if species list match OT identifier
ind_taxa[ind_taxa$approximate_match == TRUE,] # none so far

# retrieving phylogenetic relationships among taxa in the form of a trimmed sub-tree
ind_tree <- rotl::tol_induced_subtree(ott_ids = ott_id(ind_taxa), label_format = "name")

# Compute branch lengths
set.seed(1) 
ind_tree <- ape::compute.brlen(ind_tree, method = "Grafen", power = 1)
ind_tree <- ape::multi2di(ind_tree, random = TRUE) # use a randomization approach to deal with polytomies

# Check tree is ultrametric
is.ultrametric(ind_tree) # TRUE

# Create correlation matrix for analysis
ind_phylo_cor <- vcv(ind_tree, cor = T)

# Fig production
ind_tree_tip_label <- ind_tree$tip.label # extract tree tip names
ind_species_list   <- levels(ind_data_clean$species_OTL) # extract species name

# Check if lengths match for both data and tree
length(unique(ind_tree_tip_label)) # 74
length(unique(ind_species_list)) # 74

ind_data_clean$species <- gsub(" ", "_", ind_data_clean$species_OTL)
ind_data_clean$species <- factor(ind_data_clean$species, levels = ind_tree_tip_label) # relevel order by tree
ind_species_trait_data     <- ind_data_clean %>% 
  dplyr::distinct(species) %>% 
  dplyr::mutate(thermal_strategy = ind_data_clean$thermal_strategy[match(species, ind_data_clean$species)],
                taxa             = ind_data_clean$taxa[match(species, ind_data_clean$species)]) %>% 
  tibble::column_to_rownames(var = 'species')

mycol <- viridis::viridis(6) # set 6 discrete colours
diversitree::trait.plot(ind_tree, ind_species_trait_data,
                        cols = list(thermal_strategy = c("#023FA5", "#8E063B"), taxa = mycol),
                        type = 'p', cex.lab = 0.7, w = 0.05)

# OVERALL EFFECT #
# Overall effect analysis
rstan::rstan_options(auto_write = TRUE) # translate to STAN platform for running Bayesian model
options(mc.cores = parallel::detectCores()) # detects how many cores avaliable to use

ind_overall_model <- brms::brm(Zr | se(Zr_v) ~ disp_trait-1 + thermal_strategy + sex + age + origin + Zr_sei + year_centre + 
                                 (1 | es_ID) + (1 | study_ID) + (1 | species_OTL) + (1|gr(species, cov = phylo)),
                               data    = ind_data_clean,
                               family  = gaussian,
                               data2   = list(phylo = ind_phylo_cor),
                               iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                               control = list(adapt_delta = 0.999, max_treedepth = 20))

# Check convergence
brms::pp_check(ind_overall_model)

summary(ind_overall_model)
# Population-Level Effects = average effects
# Group-Level Effects = heterogeneity

brms::fixef(ind_overall_model)

# Heterogeneity
ind_overall_post <- posterior_samples(ind_overall_model) # extracting the posterior distributions from our models
parnames(ind_overall_post) # check parameter names

# WI = weight
Zr_WI <- na.omit(ind_data_clean$Zr_w)

# s2I = measurement error variance = sigma2m
s2I_Zr <- sum(Zr_WI[is.finite(Zr_WI)] * (length(Zr_WI) - 1)) / (sum(Zr_WI[is.finite(Zr_WI)]) ^ 2 - sum(Zr_WI[is.finite(Zr_WI)] ^ 2))

# Total variance, including measurement error variance
total_var_Zr <- ind_overall_post$sd_es_ID__Intercept +
  ind_overall_post$sd_species__Intercept +
  ind_overall_post$sd_species_OTL__Intercept +
  ind_overall_post$sd_study_ID__Intercept +
  s2I_Zr

# Total heterogeneity I2
I2_total_Zr     <- (total_var_Zr - s2I_Zr) / total_var_Zr
I2_total_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_total_Zr),
                           bayestestR::hdi(I2_total_Zr, ci = 0.95)$CI_low,
                           bayestestR::hdi(I2_total_Zr, ci = 0.95)$CI_high), 3) * 100

# Observational level I2
I2_esID_Zr     <- ind_overall_post$sd_es_ID__Intercept / total_var_Zr
I2_esID_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_esID_Zr),
                          bayestestR::hdi(I2_esID_Zr, ci = 0.95)$CI_low,
                          bayestestR::hdi(I2_esID_Zr, ci = 0.95)$CI_high), 3) * 100

# Study ID I2
I2_studyID_Zr     <- ind_overall_post$sd_study_ID__Intercept / total_var_Zr
I2_studyID_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_studyID_Zr),
                             bayestestR::hdi(I2_studyID_Zr, ci = 0.95)$CI_low,
                             bayestestR::hdi(I2_studyID_Zr, ci = 0.95)$CI_high), 3) * 100

# Phylogeny I2: notice that s2I is substracted from this calculation as phylogenetic
# relatedness is a "fixed random effect"
I2_phylo_Zr     <- ind_overall_post$sd_species__Intercept / (total_var_Zr - s2I_Zr)
I2_phylo_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_phylo_Zr),
                           bayestestR::hdi(I2_phylo_Zr, ci = 0.95)$CI_low,
                           bayestestR::hdi(I2_phylo_Zr, ci = 0.95)$CI_high), 3) * 100

# Species ID I2
I2_species_Zr     <- ind_overall_post$sd_species_OTL__Intercept / total_var_Zr
I2_species_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_species_Zr),
                             bayestestR::hdi(I2_species_Zr, ci = 0.95)$CI_low,
                             bayestestR::hdi(I2_species_Zr, ci = 0.95)$CI_high), 3) * 100

# Funnel plot
metafor::funnel(x = ind_data_clean$Zr, sei = ind_data_clean$Zr_sei, pch = 1)


# TRAIT - ACTVITY #
activity_trait_model <- brms::brm(Zr | se(Zr_v) ~ -1 + trait + thermal_strategy + Zr_sei + year_centre + (1 | es_ID) + (1 | study_ID) + (1 | species_OTL) + (1|gr(species, cov = phylo)),
                                  data    = ind_data_clean %>%
                                    dplyr::filter(disp_trait == "Activity") %>% 
                                    dplyr::group_by(trait) %>% 
                                    dplyr::filter(n() >= 5),
                                  family  = gaussian,
                                  data2   = list(phylo = ind_phylo_cor),
                                  iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                                  control = list(adapt_delta = 0.999, max_treedepth = 20))
# Check convergence
brms::pp_check(activity_trait_model)
summary(activity_trait_model)

act_trait_me <- marginal_effects(activity_trait_model, "trait")
act_trait_me <- as.data.frame(act_trait_me[[1]]) %>% 
  dplyr::mutate(disp_trait = "Activity",
                estimate = estimate__,
                ci.lb = lower__,
                ci.ub = upper__)

# TRAIT - EXPLORATION #
explore_trait_model <- brms::brm(Zr | se(Zr_v) ~ -1 + trait + thermal_strategy + Zr_sei + year_centre + (1 | es_ID) + (1 | study_ID) + (1 | species_OTL) + (1|gr(species, cov = phylo)),
                                 data    = ind_data_clean %>%
                                   dplyr::filter(disp_trait == "Exploration") %>% 
                                   dplyr::group_by(trait) %>% 
                                   dplyr::filter(n() >= 5),
                                 family  = gaussian,
                                 data2   = list(phylo = ind_phylo_cor),
                                 iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                                 control = list(adapt_delta = 0.999, max_treedepth = 20))
# Check convergence
brms::pp_check(explore_trait_model)
summary(explore_trait_model)

exp_trait_me <- marginal_effects(explore_trait_model, "trait")
exp_trait_me <- as.data.frame(exp_trait_me[[1]]) %>% 
  dplyr::mutate(disp_trait = "Exploration",
                estimate = estimate__,
                ci.lb = lower__,
                ci.ub = upper__)

# TRAIT - DISPERSAL - Only one trait with more than 5 effects size #
disp_trait_model <- brms::brm(Zr | se(Zr_v) ~ 1 + thermal_strategy + Zr_sei + year_centre + (1 | es_ID) + (1 | study_ID) + (1 | species_OTL) + (1|gr(species, cov = phylo)),
                              data    = ind_data_clean %>%
                                dplyr::filter(disp_trait == "Dispersal") %>% 
                                dplyr::group_by(trait) %>% 
                                dplyr::filter(n() >= 5),
                              family  = gaussian,
                              data2   = list(phylo = ind_phylo_cor),
                              iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                              control = list(adapt_delta = 0.999, max_treedepth = 20))

# Check convergence
brms::pp_check(disp_trait_model)
summary(disp_trait_model)

disp_trait_me <- as.data.frame(fixef(disp_trait_model)) %>% 
  dplyr::mutate(disp_trait = "Dispersal",
                trait = "Condition",
                estimate = Estimate,
                ci.lb = Q2.5,
                ci.ub = Q97.5)
disp_trait_me <- head(disp_trait_me, -3)

trait_me <- dplyr::bind_rows(act_trait_me, exp_trait_me, disp_trait_me) %>%
  dplyr::mutate(disp_trait  = factor(disp_trait, levels = c("Activity", "Exploration", "Dispersal")))

trait_me %>% 
  ggplot(aes(x = trait, y = estimate, group = disp_trait)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(data = ind_data_clean, aes(x = trait, y = Zr, size = Zr_inv), position = position_jitter(0.1), colour = "grey", alpha = 0.5) +
  geom_point(aes(colour = disp_trait), size = 4, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub, colour = disp_trait), size = 0.8, width = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("#7D1D67", "#E54787", "#FFA076")) +
  xlab(NULL) + ylab(expression(Effect~size~(italic(Z)[r]))) +
  coord_flip() +
  facet_wrap( ~ disp_trait) +
  mytheme() + theme(legend.position = "bottom")

# TAXA - ACTVITY (non-phylo) #
activity_taxa_model <- brms::brm(Zr | se(Zr_v) ~ -1 + taxa + Zr_sei + year_centre + (1 | es_ID) + (1 | study_ID) + (1 | species_OTL),
                                 data    = ind_data_clean %>%
                                   dplyr::filter(disp_trait == "Activity") %>% 
                                   dplyr::group_by(taxa) %>% 
                                   dplyr::filter(n() >= 5),
                                 family  = gaussian,
                                 iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                                 control = list(adapt_delta = 0.999, max_treedepth = 20))

brms::pp_check(activity_taxa_model)
brms::fixef(activity_taxa_model)

act_taxa_me <- marginal_effects(activity_taxa_model, "taxa")
act_taxa_me <- as.data.frame(act_taxa_me[[1]]) %>% 
  dplyr::mutate(disp_trait = "Activity",
                estimate = estimate__,
                ci.lb = lower__,
                ci.ub = upper__)

# TAXA - EXPLORATION (non-phylo) #
explore_taxa_model <- brms::brm(Zr | se(Zr_v) ~ -1 + taxa + Zr_sei + year_centre + (1 | es_ID) + (1 | study_ID) + (1 | species_OTL),
                                data    = ind_data_clean %>%
                                  dplyr::filter(disp_trait == "Exploration") %>% 
                                  dplyr::group_by(taxa) %>% 
                                  dplyr::filter(n() >= 5),
                                family  = gaussian,
                                iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                                control = list(adapt_delta = 0.999, max_treedepth = 20))

brms::pp_check(explore_taxa_model)
summary(explore_taxa_model)

exp_taxa_me <- marginal_effects(explore_taxa_model, "taxa")
exp_taxa_me <- as.data.frame(exp_taxa_me[[1]]) %>% 
  dplyr::mutate(disp_trait = "Exploration",
                estimate = estimate__,
                ci.lb = lower__,
                ci.ub = upper__)

# TAXA - DISPERSAL (non-phylo) #
disp_taxa_model <- brms::brm(Zr | se(Zr_v) ~ -1 + taxa + Zr_sei + year_centre + (1 | es_ID) + (1 | study_ID) + (1 | species_OTL),
                             data    = ind_data_clean %>%
                               dplyr::filter(disp_trait == "Dispersal") %>% 
                               dplyr::group_by(taxa) %>% 
                               dplyr::filter(n() >= 5),
                             family  = gaussian,
                             iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                             control = list(adapt_delta = 0.999, max_treedepth = 20))

brms::pp_check(disp_taxa_model)
summary(disp_taxa_model)

disp_taxa_me <- marginal_effects(disp_taxa_model, "taxa")
disp_taxa_me <- as.data.frame(disp_taxa_me[[1]]) %>% 
  dplyr::mutate(disp_trait = "Dispersal",
                estimate = estimate__,
                ci.lb = lower__,
                ci.ub = upper__)

taxa_me <- dplyr::bind_rows(act_taxa_me, exp_taxa_me, disp_taxa_me) %>%
  dplyr::mutate(disp_trait  = factor(disp_trait, levels = c("Activity", "Exploration", "Dispersal")))

# Plot figure
taxa_col <- c("#1A9E77", "#D76227", "#7672B2", "#E4AC25", "#E54988", "#66A643")
taxa_me %>% 
  ggplot(aes(x = taxa, y = estimate, group = disp_trait)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(data = ind_data_clean, aes(x = taxa, y = Zr, size = Zr_inv), position = position_jitter(0.1), colour = "grey", alpha = 0.3) +
  geom_point(aes(colour = taxa), size = 4, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub, colour = taxa), size = 0.8, width = 0.1, show.legend = FALSE) +
  scale_color_manual(values = taxa_col) +
  xlab(NULL) + ylab(expression(Effect~size~(italic(Z)[r]))) +
  coord_flip() +
  facet_wrap( ~ disp_trait) +
  mytheme() + theme(legend.position = "bottom")

# METABOLISM #
# Subset metabolism trait
mr_data <- ind_data_clean %>%
  dplyr::filter(trait == "Metabolism") %>%
  droplevels()

# Grouping MR
mr_data <- mr_data %>%
  dplyr::mutate(group = dplyr::recode(response, 
                                      "Routine MR"  = "Active MR", 
                                      "Field MR"    = "Active MR",
                                      "Standard MR" = "Inactive MR",
                                      "Resting MR"  = "Inactive MR",
                                      "Basal MR"    = "Inactive MR",
                                      "Muscle citrate synthase"      = "Metabolic enzyme",
                                      "Muscle cytochrome c oxidase"  = "Metabolic enzyme",
                                      "Muscle lactate dehydrogenase" = "Metabolic enzyme",
                                      "Liver citrate synthase"       = "Metabolic enzyme",
                                      "Liver cytochrome c oxidase"   = "Metabolic enzyme"))

mr_model <- brms::brm(Zr | se(Zr_v) ~ -1 + group + age + Zr_sei + year_centre + (1 | es_ID) + (1 | study_ID) + (1 | species_OTL) + (1|gr(species, cov = phylo)),
                      data    = mr_data %>%
                        dplyr::group_by(group) %>% 
                        dplyr::filter(n() >= 5),
                      family  = gaussian,
                      data2   = list(phylo = ind_phylo_cor),
                      iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                      control = list(adapt_delta = 0.999, max_treedepth = 20))
summary(mr_model)

mr_model_me <- marginal_effects(mr_model, "group")
mr_model_me <- as.data.frame(mr_model_me[[1]]) %>% 
  dplyr::mutate(estimate = estimate__,
                ci.lb = lower__,
                ci.ub = upper__)

mr_plot <- mr_model_me %>% 
  ggplot(aes(x = group, y = estimate)) +
  geom_jitter(data = mr_data, aes(x = group, y = Zr, size = Zr_inv), position = position_jitter(0.1), colour = "grey", alpha = 0.3) +
  geom_point(size = 4, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  xlab(NULL) + ylab(expression(Effect~size~(italic(Z)[r]))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  mytheme() + theme(legend.position = "bottom")

# LOCOMOTOR CAPACITY #
# Subset metabolism trait
loco_data <- ind_data_clean %>%
  dplyr::filter(trait == "Locomotor capacity") %>%
  droplevels()

loco_model <- brms::brm(Zr | se(Zr_v) ~ -1 + response + thermal_strategy + Zr_sei + year_centre + (1 | es_ID) + (1 | study_ID) + (1 | species_OTL) + (1|gr(species, cov = phylo)),
                        data    = loco_data %>%
                          dplyr::group_by(response) %>% 
                          dplyr::filter(n() >= 5),
                        family  = gaussian,
                        data2   = list(phylo = ind_phylo_cor),
                        iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                        control = list(adapt_delta = 0.999, max_treedepth = 20))

summary(loco_model)

loco_model_me <- marginal_effects(loco_model, "response")
loco_model_me <- as.data.frame(loco_model_me[[1]]) %>% 
  dplyr::mutate(estimate = estimate__,
                ci.lb = lower__,
                ci.ub = upper__)

loco_plot <- loco_model_me %>% 
  ggplot(aes(x = response, y = estimate)) +
  geom_jitter(data = loco_data, aes(x = response, y = Zr, size = Zr_inv), position = position_jitter(0.1), colour = "grey", alpha = 0.3) +
  geom_point(size = 4, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  xlab(NULL) + ylab(expression(Effect~size~(italic(Z)[r]))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(position = "top") +
  coord_flip() +
  mytheme() + theme(legend.position = "bottom")

plot_grid(mr_plot, loco_plot, ncol = 2, align = 'v', axis = 'lr')

## 2. POPULATION-LEVEL ANALYSIS ##-------------------------------------------------------------------------------
# TEMPERATURE & DISPERSAL MODEL #
temp_prior <- prior(normal(0, 5), class = "Intercept") + prior(normal(0, 5), class = "b")

# Regression model, let the residuals of standard deviation (sigma) vary with the predictor
temp_model <- brms::brm(bf(lnRate ~ temp_diff * disp_mode, sigma ~ temp_diff),
                        data    = pop_data_clean,
                        family  = gaussian,
                        prior   = temp_prior,
                        iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                        control = list(adapt_delta = 0.999, max_treedepth = 20))

# Model summary
summary(temp_model)

# Extract diserpsal mode model predictions
rate_marg_eff <- as.data.frame(ggeffects::ggpredict(temp_model, terms = c("disp_mode")))

# Plot
rate_plot <- ggplot(data = rate_marg_eff, aes(x = x, y = exp(predicted))) +
  geom_jitter(data = pop_data_clean, aes(x = disp_mode, y = exp(lnRate)), colour = "grey", size = 2, alpha = 0.3, position = position_jitter(0.02)) +
  geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high), colour = x), size = 0.8, width = 0.1) +
  geom_point(aes(colour = x), size = 3) +
  ylab(expression("Rate of dispersal (km y"^"-1"*")")) +
  xlab(NULL) +
  scale_color_manual(values = c("#FFA076", "#E54787", "#7D1D67")) +
  mytheme()

# Extract temp difference model predictions
temp_marg_eff <- as.data.frame(ggeffects::ggpredict(temp_model, terms = c("temp_diff[sample=50]")))
ggeffects::ggpredict(temp_model, terms = c("temp_diff[sample=50]"))

# Plot
temp_plot <- ggplot(data = temp_marg_eff, aes(x = x, y = exp(predicted))) +
  geom_ribbon(aes(x = x, ymin = exp(conf.low), ymax = exp(conf.high)), fill = "#7D1D67", alpha = 0.1) +
  geom_line(aes(), size = 1, linetype = "dashed", show.legend = FALSE) +
  geom_point(data = pop_data_clean, aes(x = temp_diff, y = exp(lnRate), colour = disp_mode), size = 2) +
  ylab(expression("Rate of dispersal (km y"^"-1"*")")) +
  xlab("Temperature difference between core and expansion front (°C)") +
  scale_color_manual(values = c("#FFA076", "#E54787", "#7D1D67")) +
  mytheme()

disp_plot <- cowplot::plot_grid(rate_plot + theme(legend.position = "none"), 
                                temp_plot + theme(legend.position = "none"), 
                                ncol = 2, align = 'h', axis = 'tb')

# RAINFALL MODEL #
rain_prior <- prior(normal(0, 5), class = "Intercept") + prior(normal(0, 5), class = "b")

# Regression model, let the residuals of standard deviation (sigma) vary with the predictor
rain_model <- brms::brm(bf(dispersal_rate_km_y ~ rainfall_diff * disp_mode),
                        data    = pop_data_clean,
                        family  = gaussian,
                        prior   = rain_prior,
                        iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                        control = list(adapt_delta = 0.999, max_treedepth = 20))

# Model summary
summary(rain_model)

# lnRR & TIME DIFFERENCE MODEL #
time_prior <- prior(normal(0, 5), class = "Intercept") + prior(normal(0, 1), class = "b")
time_model <- brms::brm(bf(abslnRR ~ lnTime * disp_mode + (1 | study_ID) + (1 | species), sigma ~ lnTime),
                        data    = pop_data_clean,
                        family  = hurdle_lognormal,
                        prior   = time_prior,
                        iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                        control = list(adapt_delta = 0.999, max_treedepth = 20))

# Check model convergence
brms::pp_check(time_model)

# Model summary
summary(time_model) # no difference between dispersal mode. Therefore visualise overall model

# Extract model predictions
time_marg_eff <- as.data.frame(ggeffects::ggpredict(time_model, terms = c("lnTime[sample=40]"))) 

# Plot 
time_plot <- ggplot(data = time_marg_eff, aes(x = x, y = predicted)) +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "#7D1D67",alpha = 0.1) +
  geom_line(aes(), size = 1, linetype = "dashed", show.legend = FALSE) +
  geom_point(data = pop_data_clean, aes(x = lnTime, y = abslnRR, colour = disp_mode), size = 2) +
  ylab("Absolute effect size (lnRR)") +
  xlab("Time since divergence (years)") +
  scale_x_continuous(breaks = c(0, 1.609438, 2.302585, 3.218876, 4.317488), 
                     labels = round(c(exp(0), exp(1.609438), exp(2.302585), exp(3.218876), exp(4.317488)), 1)) +
  scale_color_manual(values = c("#FFA076", "#E54787", "#7D1D67")) +
  mytheme()

cowplot::plot_grid(disp_plot + theme(legend.position = "none"), 
                   time_plot + theme(legend.position = "none"), 
                   nrow = 2, labels = "AUTO")

# BUILD PHYLOGENY ALL #
pop_species <- sort(unique(as.character(pop_data_clean$species_OTL))) # generate list of species (as character format)
pop_taxa    <- rotl::tnrs_match_names(names = pop_species) # match taxonomic names to the OTL

# check if species list match OT identifier
pop_taxa[pop_taxa$approximate_match == TRUE,] # none so far

# retrieving phylogenetic relationships among taxa in the form of a trimmed sub-tree
pop_tree <- rotl::tol_induced_subtree(ott_ids = ott_id(pop_taxa), label_format = "name")

# Compute branch lengths
set.seed(1) 
pop_tree <- ape::compute.brlen(pop_tree, method = "Grafen", power = 1)
pop_tree <- ape::multi2di(pop_tree, random = TRUE) # use a randomization approach to deal with polytomies

# Check tree is ultrametric
is.ultrametric(pop_tree) # TRUE

# Create correlation matrix for analysis
pop_phylo_cor <- vcv(pop_tree, cor = T)

# Fig production
pop_tree_tip_label <- pop_tree$tip.label # extract tree tip names
pop_species_list   <- levels(pop_data_clean$species_OTL) # extract species name from diet analysis

# Check if lengths match for both data and tree
length(unique(pop_tree_tip_label)) # 15
length(unique(pop_species_list)) # 15

pop_data_clean$species <- gsub(" ", "_", pop_data_clean$species_OTL)
pop_data_clean$species <- factor(pop_data_clean$species, levels = pop_tree_tip_label) # relevel order by tree
pop_species_trait_data     <- pop_data_clean %>% 
  dplyr::distinct(species) %>% 
  dplyr::mutate(thermal_strategy = pop_data_clean$thermal_strategy[match(species, pop_data_clean$species)],
                taxa             = pop_data_clean$taxa[match(species, pop_data_clean$species)]) %>% 
  tibble::column_to_rownames(var = 'species')

mycol <- viridis::viridis(4) # set 6 discrete colours
diversitree::trait.plot(pop_tree, pop_species_trait_data,
                        cols = list(thermal_strategy = c("#023FA5", "#8E063B"), taxa = mycol),
                        type = 'p', cex.lab = 0.7, w = 0.05)

# OVERALL MODEL #
pop_prior <- c(set_prior("cauchy(0, 1)", class = "sd")) # Cauchy on tau (random effect variance), normal on fixed effect

# Account for unbalanced sampling (square of inverse of effective sample size), and time-lag bias
pop_overall_model <- brms::brm(lnRR | se(v) ~ 1 + sqrt_inv_eff + year_centre + (1 | es_ID) + (1 | study_ID) + (1 | species_OTL) + (1|gr(species, cov = phylo)),
                               data    = pop_data_clean,
                               family  = gaussian,
                               data2   = list(phylo = pop_phylo_cor),
                               prior   = pop_prior,
                               iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                               control = list(adapt_delta = 0.999, max_treedepth = 20))

# Check convergence
brms::pp_check(pop_overall_model)

# Model summary
summary(pop_overall_model)
# Population-Level Effects = average effects
# Group-Level Effects = SD

# Heterogeneity
pop_overall_post <- posterior_samples(pop_overall_model) # extracting the posterior distributions from our models
parnames(pop_overall_post) # check parameter names

# WI = weight
Zr_WI <- na.omit(1 / pop_data_clean$v)

# s2I = measurement error variance = sigma2m
s2I_Zr <- sum(Zr_WI[is.finite(Zr_WI)] * (length(Zr_WI) - 1)) / (sum(Zr_WI[is.finite(Zr_WI)]) ^ 2 - sum(Zr_WI[is.finite(Zr_WI)] ^ 2))

# Total variance, including measurement error variance
total_var_Zr <- pop_overall_post$sd_es_ID__Intercept +
  pop_overall_post$sd_species__Intercept +
  pop_overall_post$sd_species_OTL__Intercept +
  pop_overall_post$sd_study_ID__Intercept +
  s2I_Zr

# Total heterogeneity I2
I2_total_Zr     <- (total_var_Zr - s2I_Zr) / total_var_Zr
I2_total_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_total_Zr),
                           bayestestR::hdi(I2_total_Zr, ci = 0.95)$CI_low,
                           bayestestR::hdi(I2_total_Zr, ci = 0.95)$CI_high), 3) * 100

# Observational level I2
I2_esID_Zr     <- pop_overall_post$sd_es_ID__Intercept / total_var_Zr
I2_esID_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_esID_Zr),
                          bayestestR::hdi(I2_esID_Zr, ci = 0.95)$CI_low,
                          bayestestR::hdi(I2_esID_Zr, ci = 0.95)$CI_high), 3) * 100

# Study ID I2
I2_studyID_Zr     <- pop_overall_post$sd_study_ID__Intercept / total_var_Zr
I2_studyID_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_studyID_Zr),
                             bayestestR::hdi(I2_studyID_Zr, ci = 0.95)$CI_low,
                             bayestestR::hdi(I2_studyID_Zr, ci = 0.95)$CI_high), 3) * 100

# Phylogeny I2: notice that s2I is substracted from this calculation as phylogenetic
# relatedness is a "fixed random effect"
I2_phylo_Zr     <- pop_overall_post$sd_species__Intercept / (total_var_Zr - s2I_Zr)
I2_phylo_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_phylo_Zr),
                           bayestestR::hdi(I2_phylo_Zr, ci = 0.95)$CI_low,
                           bayestestR::hdi(I2_phylo_Zr, ci = 0.95)$CI_high), 3) * 100

# Species ID I2
I2_species_Zr     <- pop_overall_post$sd_species_OTL__Intercept / total_var_Zr
I2_species_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_species_Zr),
                             bayestestR::hdi(I2_species_Zr, ci = 0.95)$CI_low,
                             bayestestR::hdi(I2_species_Zr, ci = 0.95)$CI_high), 3) * 100

# Funnel plot
metafor::funnel(x = pop_data_clean$lnRR, sei = pop_data_clean$sei, pch = 1,
                ylab = "Precision (1/SE)",
                xlab = "Effect size (lnRR)")

# TRAIT MODEL #
pop_trait_model <- brms::brm(lnRR | se(v) ~ -1 + trait + sex + sqrt_inv_eff + year_centre + (1 | es_ID) + (1 | study_ID) + (1 | species_OTL) + (1|gr(species, cov = phylo)),
                             data    = pop_data_clean %>%
                               dplyr::group_by(trait) %>% 
                               dplyr::filter(n() >= 5),
                             family  = gaussian,
                             data2   = list(phylo = pop_phylo_cor),
                             prior   = pop_prior,
                             iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                             control = list(adapt_delta = 0.999, max_treedepth = 20))

# Check convergence
brms::pp_check(pop_trait_model)

# Model summary
summary(pop_trait_model)

# Extract marginal effects from model
pop_trait_me <- marginal_effects(pop_trait_model, "trait")
pop_trait_me <- as.data.frame(pop_trait_me[[1]]) %>% 
  dplyr::mutate(estimate = estimate__, ci.lb = lower__, ci.ub = upper__)

# Plot model output
pop_trait_plot <- pop_trait_me %>% 
  ggplot(aes(x = trait, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(data = pop_data_clean, aes(x = trait, y = lnRR, size = sei_inv), position = position_jitter(0.1), colour = "grey", alpha = 0.5) +
  geom_point(aes(colour = trait), size = 4, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub, colour = trait), size = 0.8, width = 0.1, show.legend = FALSE) +
  colorspace::scale_colour_discrete_sequential(palette = "SunsetDark") +
  xlab(NULL) + ylab("Effect size (lnRR)") +
  coord_flip() +
  mytheme() + theme(legend.position = "bottom")


# TAXA MODEL - non-phylo #
pop_taxa_model <- brms::brm(lnRR | se(v) ~ -1 + taxa + sex + sqrt_inv_eff + year_centre + (1 | es_ID) + (1 | study_ID) + (1 | species_OTL),
                            data    = pop_data_clean %>%
                              dplyr::group_by(taxa) %>% 
                              dplyr::filter(n() >= 5),
                            family  = gaussian,
                            prior   = pop_prior,
                            iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                            control = list(adapt_delta = 0.999, max_treedepth = 20))

# Check convergence
brms::pp_check(pop_taxa_model)

# Model summary
summary(pop_taxa_model)

# Extract marginal effects from model
pop_taxa_me <- marginal_effects(pop_taxa_model, "taxa")
pop_taxa_me <- as.data.frame(pop_taxa_me[[1]]) %>% 
  dplyr::mutate(estimate = estimate__, ci.lb = lower__, ci.ub = upper__)

# Plot model output
pop_taxa_plot <- pop_taxa_me %>% 
  ggplot(aes(x = taxa, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(data = pop_data_clean, aes(x = taxa, y = lnRR, size = sei_inv), position = position_jitter(0.1), colour = "grey", alpha = 0.5) +
  geom_point(aes(colour = taxa), size = 4, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub, colour = taxa), size = 0.8, width = 0.1, show.legend = FALSE) +
  colorspace::scale_colour_discrete_sequential(palette = "Viridis") +
  xlab(NULL) + ylab("Effect size (lnRR)") +
  scale_x_discrete(position = "top") +
  coord_flip() +
  mytheme() + theme(legend.position = "bottom")

plot_grid(pop_trait_plot, pop_taxa_plot, ncol = 2, align = 'v', axis = 'lr')
