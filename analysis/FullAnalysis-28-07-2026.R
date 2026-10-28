###############################################
# Libraries
###############################################
library(dplyr)
library(tidyr)
library(tidyverse)
library(writexl)
library(ggplot2)
library(openxlsx)
library(stringr)
library(DBI)
library(RSQLite)
library(patchwork)
library(scales)

###############################################
# Load and Prepare Data
###############################################
# NTG Data
setwd("C:/Users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/InSilico/N_Drug/Data/")
#setwd("C:/Users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/InSilico/N_Drug/Data/")

NTG <- read.csv("BurnInFinal.csv")

# N target Data
NTreat <- read.csv("TreatmentFinal.csv")
NRelap <- read.csv("RelapseFinal.csv")
min(NRelap$Final_N)

# S target Data
setwd("C:/Users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/InSilico/S_Drug/Data/")
#setwd("C:/Users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/InSilico/S_Drug/Data/")

STreat <- read.csv("TreatmentFinal.csv")
SRelap <- read.csv("RelapseFinal.csv")
min(SRelap$Final_S)

###############################################
# Output directory
###############################################
out_dir <- "C:/Users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Plots/"
#out_dir <- "C:/Users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Plots/"

# Helper to save plots consistently
# Updated save_plot helper — adds title wrapping via theme
save_plot <- function(filename, plot, width = 174, height = 234, units = "mm", dpi = 300) {
  plot <- plot +
    theme(axis.text=element_text(size=12),plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 11))
  ggsave(paste0(out_dir, filename), plot = plot,
         width = width, height = height, units = units, dpi = dpi)
}

###############################################
# Define Pattern Labels (used throughout)
###############################################
pat_labels <- c(
  "<N>|I>|S|" = "P1",
  "<N>|I><S|" = "P2",
  "|N|<I>|S>" = "P3",
  "|N|<I>|S|" = "P4",
  "|N|<I><S>" = "P5"
)

PatInt           <- c("<N>|I>|S|", "<N>|I><S|", "|N|<I>|S>", "|N|<I>|S|", "|N|<I><S>")
pat_label_levels <- unname(pat_labels)
N_OptPatterns    <- c("<N>|I>|S|", "<N>|I><S|")
S_OptPatterns    <- c("|N|<I>|S>", "|N|<I>|S|", "|N|<I><S>")

# Shared colour palette
Icol  <- "#F8766D"
Ncol  <- "#7CAE00"
nocol <- "#00BFC4"
Scol  <- "#C77CFF"

###############################################
# Processing N data
###############################################
NTreat$Drug <- NTreat$Drug * 2
NRelap$Drug <- NRelap$Drug * 2

## Filter and Compute Relative Reduction
NTreat <- NTreat %>%
  mutate(N_Rel_reduction = (Initial_N - Final_N) / Initial_N)
NRelap <- NRelap %>%
  mutate(N_Rel_reduction = (Initial_N - Final_N) / Initial_N)

NTreat <- NTreat %>%
  mutate(K3Rate = case_when(
    K3NS==0 & K3NI==0 & K3IN==0 &
      K3IS==0 & K3SI==0 & K3SN==0 ~ "0",
    K3NS==0.1 | K3NI==0.1 | K3IN==0.1 |
      K3IS==0.1 | K3SI==0.1 | K3SN==0.1 ~ "0.1",
    K3NS==0.5 | K3NI==0.5 | K3IN==0.5 |
      K3IS==0.5 | K3SI==0.5 | K3SN==0.5 ~ "0.5",
    K3NS==0.9 | K3NI==0.9 | K3IN==0.9 |
      K3IS==0.9 | K3SI==0.9 | K3SN==0.9 ~ "0.9"
  )) %>%
  mutate(K3Rate = factor(K3Rate, levels = c("0","0.1","0.5","0.9")))

NRelap <- NRelap %>%
  mutate(K3Rate = case_when(
    K3NS==0 & K3NI==0 & K3IN==0 &
      K3IS==0 & K3SI==0 & K3SN==0 ~ "0",
    K3NS==0.1 | K3NI==0.1 | K3IN==0.1 |
      K3IS==0.1 | K3SI==0.1 | K3SN==0.1 ~ "0.1",
    K3NS==0.5 | K3NI==0.5 | K3IN==0.5 |
      K3IS==0.5 | K3SI==0.5 | K3SN==0.5 ~ "0.5",
    K3NS==0.9 | K3NI==0.9 | K3IN==0.9 |
      K3IS==0.9 | K3SI==0.9 | K3SN==0.9 ~ "0.9"
  )) %>%
  mutate(K3Rate = factor(K3Rate, levels = c("0","0.1","0.5","0.9")))

# Count occurrences of positive RR and score
score_table_NT <- NTreat %>%
  group_by(Interconversion) %>%
  summarise(Score = sum(N_Rel_reduction >= 0, na.rm = TRUE), .groups = "drop")

score_table_NR <- NRelap %>%
  group_by(Interconversion) %>%
  summarise(Score = sum(N_Rel_reduction >= 0, na.rm = TRUE), .groups = "drop")

###############################################
# Processing S data
###############################################
STreat$Drug <- STreat$Drug * 2
SRelap$Drug <- SRelap$Drug * 2

## Filter and Compute Relative Reduction
STreat <- STreat %>%
  mutate(S_Rel_reduction = (Initial_S - Final_S) / Initial_S)
SRelap <- SRelap %>%
  mutate(S_Rel_reduction = (Initial_S - Final_S) / Initial_S)

STreat <- STreat %>%
  mutate(K3Rate = case_when(
    K3NS==0 & K3NI==0 & K3IN==0 &
      K3IS==0 & K3SI==0 & K3SN==0 ~ "0",
    K3NS==0.1 | K3NI==0.1 | K3IN==0.1 |
      K3IS==0.1 | K3SI==0.1 | K3SN==0.1 ~ "0.1",
    K3NS==0.5 | K3NI==0.5 | K3IN==0.5 |
      K3IS==0.5 | K3SI==0.5 | K3SN==0.5 ~ "0.5",
    K3NS==0.9 | K3NI==0.9 | K3IN==0.9 |
      K3IS==0.9 | K3SI==0.9 | K3SN==0.9 ~ "0.9"
  )) %>%
  mutate(K3Rate = factor(K3Rate, levels = c("0","0.1","0.5","0.9")))

SRelap <- SRelap %>%
  mutate(K3Rate = case_when(
    K3NS==0 & K3NI==0 & K3IN==0 &
      K3IS==0 & K3SI==0 & K3SN==0 ~ "0",
    K3NS==0.1 | K3NI==0.1 | K3IN==0.1 |
      K3IS==0.1 | K3SI==0.1 | K3SN==0.1 ~ "0.1",
    K3NS==0.5 | K3NI==0.5 | K3IN==0.5 |
      K3IS==0.5 | K3SI==0.5 | K3SN==0.5 ~ "0.5",
    K3NS==0.9 | K3NI==0.9 | K3IN==0.9 |
      K3IS==0.9 | K3SI==0.9 | K3SN==0.9 ~ "0.9"
  )) %>%
  mutate(K3Rate = factor(K3Rate, levels = c("0","0.1","0.5","0.9")))

# Count occurrences of positive RR and score
score_table_ST <- STreat %>%
  group_by(Interconversion) %>%
  summarise(Score = sum(S_Rel_reduction >= 0, na.rm = TRUE), .groups = "drop")

score_table_SR <- SRelap %>%
  group_by(Interconversion) %>%
  summarise(Score = sum(S_Rel_reduction >= 0, na.rm = TRUE), .groups = "drop")

###############################################
# Processing NTG data
###############################################
## Calculate mean final cell populations for three cell types for each interconversion pattern

# 1. Classify K3Rate
NTG <- NTG %>%
  mutate(K3Rate = case_when(
    K3NS==0 & K3NI==0 & K3IN==0 &
      K3IS==0 & K3SI==0 & K3SN==0 ~ "0",
    K3NS==0.1 | K3NI==0.1 | K3IN==0.1 |
      K3IS==0.1 | K3SI==0.1 | K3SN==0.1 ~ "0.1",
    K3NS==0.5 | K3NI==0.5 | K3IN==0.5 |
      K3IS==0.5 | K3SI==0.5 | K3SN==0.5 ~ "0.5",
    K3NS==0.9 | K3NI==0.9 | K3IN==0.9 |
      K3IS==0.9 | K3SI==0.9 | K3SN==0.9 ~ "0.9"
  )) %>%
  mutate(K3Rate = factor(K3Rate, levels = c("0","0.1","0.5","0.9")))

# 2. Summaries
df_NTG <- NTG %>%
  group_by(Interconversion) %>%
  summarise(
    Mean_N = mean(Final_N, na.rm = TRUE),
    Mean_I = mean(Final_I, na.rm = TRUE),
    Mean_S = mean(Final_S, na.rm = TRUE)
  )

###############################################
# Ranking scores
###############################################
OptN <- df_NTG[which(df_NTG$Mean_N < 0.05),]
OptS <- df_NTG[which(df_NTG$Mean_S < 0.05),]

OptN <- OptN %>%
  left_join(score_table_ST %>% select(Interconversion, Score), by = "Interconversion") %>%
  rename(RRST = Score) %>%
  left_join(score_table_SR %>% select(Interconversion, Score), by = "Interconversion") %>%
  rename(RRSR = Score)

OptS <- OptS %>%
  left_join(score_table_NT %>% select(Interconversion, Score), by = "Interconversion") %>%
  rename(RRNT = Score) %>%
  left_join(score_table_NR %>% select(Interconversion, Score), by = "Interconversion") %>%
  rename(RRNR = Score)

###############################################
# Heatmaps of Cell pops (Patterns of Interest)
###############################################

# 1. Filter by pattern of interest
NTG_RR_PT <- NTG %>%
  filter(Interconversion %in% PatInt)

# 2. Summaries
df_NTG_long <- NTG_RR_PT %>%
  group_by(Interconversion, K3Rate) %>%
  summarise(
    Mean_N = mean(Final_N, na.rm = TRUE),
    Mean_I = mean(Final_I, na.rm = TRUE),
    Mean_S = mean(Final_S, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = starts_with("Mean_"), names_to = "State", values_to = "Value") %>%
  mutate(
    State = recode(State, "Mean_N" = "N", "Mean_I" = "I", "Mean_S" = "S"),
    State = factor(State, levels = c("N", "I", "S")),
    K3Rate = factor(K3Rate, levels = c("0", "0.1", "0.5", "0.9")),
    x_pos = as.numeric(K3Rate),
    Interconversion = recode(Interconversion, !!!pat_labels),
    Interconversion = factor(Interconversion, levels = pat_label_levels)
  )

# 3. Heatmap plot
p_ntg_heatmap <- ggplot(df_NTG_long, aes(x = x_pos, y = Interconversion)) +
  geom_tile(aes(fill = State, alpha = Value), width = 0.9, height = 0.9, colour = "white") +
  geom_text(aes(label = formatC(Value, digits = 2)), size = 5, fontface = "bold", colour = "black") +
  scale_fill_manual(values = c("N" = Ncol, "I" = Icol, "S" = Scol)) +
  scale_alpha(range = c(0.1, 1)) +
  scale_x_continuous(
    breaks = 1:length(levels(df_NTG_long$K3Rate)),
    labels = levels(df_NTG_long$K3Rate)
  ) +
  scale_y_discrete(limits = rev(pat_label_levels)) +
  facet_wrap(~ State, ncol = 1, strip.position = "top") +
  theme_classic(base_size = 14) +
  theme(axis.text=element_text(size=12),axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(
    x     = "Interconversion Rate (K3Rate)",
    y     = "Interconversion Pattern",
    fill  = "Cell Type",
    title = "Mean Cell Populations at the \nEnd of NTG (Faceted by Cell Type)",
    alpha = "Mean Cell Population"
  )

save_plot("NTG_heatmap.pdf", p_ntg_heatmap)

# 4. Bar chart plot
p_ntg_bar <- ggplot(df_NTG_long, aes(x = K3Rate, y = Value, fill = State)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("N" = Ncol, "I" = Icol, "S" = Scol)) +
  facet_wrap(~ Interconversion, ncol = 1, strip.position = "right") +
  theme_classic(base_size = 12) +
  theme(axis.text.y=element_text(size=7),axis.text.x = element_text(size=14,angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(
    x     = "Interconversion Rate",
    y     = "Mean Cell Population",
    fill  = "Cell Type",
    title = "Mean Cell Populations at the \nEnd of NTG by Interconversion Pattern"
  )

save_plot("NTG_barchart.pdf", p_ntg_bar)

###############################################
# Cross-inhibitor Score Heatmaps
###############################################

# Shared heatmap functions
create_heatmap_drug <- function(data, title, pattern_levels) {
  data_long <- data %>%
    pivot_longer(cols = -Interconversion, names_to = c("Phase", "Drug"),
                 names_sep = "_", values_to = "Score") %>%
    mutate(
      Phase = factor(Phase, levels = c("Treatment", "Relapse")),
      Interconversion = recode(Interconversion, !!!pat_labels),
      Interconversion = factor(Interconversion, levels = rev(pat_label_levels))
    )
  ggplot(data_long, aes(x = Drug, y = Interconversion, fill = Score)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "steelblue", mid = "white", high = "firebrick",
      midpoint = 6, limits = c(0, 12)
    ) +
    labs(title = title, x = "Drug", y = "Interconversion Pattern", fill = "Score") +
    geom_text(aes(label = Score), color = "black", size = 5) +
    facet_wrap(~ Phase, nrow = 1) +
    theme_classic() +
    theme(axis.text=element_text(size=12),axis.text.x = element_text(angle = 45, hjust = 1))
}

create_heatmap_k3 <- function(data, title, pattern_levels) {
  data_long <- data %>%
    pivot_longer(cols = -Interconversion, names_to = c("Phase", "K3Rate"),
                 names_sep = "_", values_to = "Score") %>%
    mutate(
      Phase = factor(Phase, levels = c("Treatment", "Relapse")),
      Interconversion = recode(Interconversion, !!!pat_labels),
      Interconversion = factor(Interconversion, levels = rev(pat_label_levels))
    )
  ggplot(data_long, aes(x = K3Rate, y = Interconversion, fill = Score)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "steelblue", mid = "white", high = "firebrick",
      midpoint = 6, limits = c(0, 12)
    ) +
    labs(title = title, x = "Interconversion Rate", y = "Interconversion Pattern", fill = "Score") +
    geom_text(aes(label = Score), color = "black", size = 5) +
    facet_wrap(~ Phase, nrow = 1) +
    theme_classic() +
    theme(axis.text=element_text(size=12),axis.text.x = element_text(angle = 45, hjust = 1))
}

#### N-Supressing patterns under S inhibitor ####
S_Treat_NOpt <- STreat %>% filter(Interconversion %in% N_OptPatterns)
S_Relap_NOpt <- SRelap %>% filter(Interconversion %in% N_OptPatterns)

score_S_NOpt_T <- S_Treat_NOpt %>%
  group_by(Drug, Interconversion) %>%
  summarise(Score = sum(S_Rel_reduction >= 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(Type = "Treatment")

score_S_NOpt_R <- S_Relap_NOpt %>%
  group_by(Drug, Interconversion) %>%
  summarise(Score = sum(S_Rel_reduction >= 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(Type = "Relapse")

result_S_NOpt_RR_Drug <- bind_rows(score_S_NOpt_T, score_S_NOpt_R) %>%
  pivot_wider(names_from = c(Type, Drug), values_from = Score)

score_S_NOpt_K3_T <- S_Treat_NOpt %>%
  group_by(K3Rate, Interconversion) %>%
  summarise(Score = sum(S_Rel_reduction >= 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(Type = "Treatment")

score_S_NOpt_K3_R <- S_Relap_NOpt %>%
  group_by(K3Rate, Interconversion) %>%
  summarise(Score = sum(S_Rel_reduction >= 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(Type = "Relapse")

result_S_NOpt_RR_K3 <- bind_rows(score_S_NOpt_K3_T, score_S_NOpt_K3_R) %>%
  pivot_wider(names_from = c(Type, K3Rate), values_from = Score)

p_NOpt_drug <- create_heatmap_drug(result_S_NOpt_RR_Drug,
                                   "N-Supressing Patterns under S Inhibitor: \nRR Score Across Drug Dose",
                                   N_OptPatterns)
p_NOpt_k3   <- create_heatmap_k3(result_S_NOpt_RR_K3,
                                 "N-Supressing Patterns under S Inhibitor: \nRR Score Across Interconversion Rate",
                                 N_OptPatterns)

save_plot("NOpt_SInhib_RR_Drug.pdf",  p_NOpt_drug)
save_plot("NOpt_SInhib_RR_K3.pdf",    p_NOpt_k3)

#### S-Supressing patterns under N inhibitor ####
N_Treat_SOpt <- NTreat %>% filter(Interconversion %in% S_OptPatterns)
N_Relap_SOpt <- NRelap %>% filter(Interconversion %in% S_OptPatterns)

score_N_SOpt_T <- N_Treat_SOpt %>%
  group_by(Drug, Interconversion) %>%
  summarise(Score = sum(N_Rel_reduction >= 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(Type = "Treatment")

score_N_SOpt_R <- N_Relap_SOpt %>%
  group_by(Drug, Interconversion) %>%
  summarise(Score = sum(N_Rel_reduction >= 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(Type = "Relapse")

result_N_SOpt_RR_Drug <- bind_rows(score_N_SOpt_T, score_N_SOpt_R) %>%
  pivot_wider(names_from = c(Type, Drug), values_from = Score)

score_N_SOpt_K3_T <- N_Treat_SOpt %>%
  group_by(K3Rate, Interconversion) %>%
  summarise(Score = sum(N_Rel_reduction >= 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(Type = "Treatment")

score_N_SOpt_K3_R <- N_Relap_SOpt %>%
  group_by(K3Rate, Interconversion) %>%
  summarise(Score = sum(N_Rel_reduction >= 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(Type = "Relapse")

result_N_SOpt_RR_K3 <- bind_rows(score_N_SOpt_K3_T, score_N_SOpt_K3_R) %>%
  pivot_wider(names_from = c(Type, K3Rate), values_from = Score)

p_SOpt_drug <- create_heatmap_drug(result_N_SOpt_RR_Drug,
                                   "S-Supressing Patterns under N Inhibitor: \nRR Score Across Drug Dose",
                                   S_OptPatterns)
p_SOpt_k3   <- create_heatmap_k3(result_N_SOpt_RR_K3,
                                 "S-Supressing Patterns under N Inhibitor: \nRR Score Across Interconversion Rate",
                                 S_OptPatterns)

save_plot("SOpt_NInhib_RR_Drug.pdf", p_SOpt_drug)
save_plot("SOpt_NInhib_RR_K3.pdf",   p_SOpt_k3)

###############################################
# Stability Analysis
###############################################



# --- LOAD DATA ---
#setwd("C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Stability/NoteBookOut")
setwd("C:/Users/Student/OneDrive/Documents/Year_3/Experiments/Stability/NoteBookOut")
NTGStab <- read.csv("stability_NTGFULL.csv")

# --- EXTRACT MAX REAL EIGENVALUE PER ROW ---
NTGStab <- NTGStab %>%
  mutate(
    Max_Re_Eig = pmax(eig_1, eig_2, eig_3)
  )

# --- EXTRACT K3 RATE ---
# All active k3 values are the same within a pattern so just take the max
#df <- df %>%
# mutate(
#  K3_rate = pmax(K3NS, K3NI, K3IN, K3IS, K3SI, K3SN)
#)

# --- KEEP ONE ROW PER PATTERN + K3 RATE ---
heatmap_NTGStab <- NTGStab %>%
  select(pattern, K3, Max_Re_Eig) %>%
  distinct()

# --- SET PATTERN ORDER ---

pattern_order <- c( "pattern_1", "pattern_2", "pattern_3", "pattern_4", "pattern_5")

# Define label mapping
pat_labels <- c(
  "pattern_1" = "P1",
  "pattern_2" = "P2",
  "pattern_3" = "P3",
  "pattern_4" = "P4",
  "pattern_5" = "P5"
)

PatInt <- c("pattern_1", "pattern_2", "pattern_3", "pattern_4", "pattern_5")
pat_label_levels <- unname(pat_labels)

heatmap_df <- heatmap_NTGStab %>%
  mutate(
    pattern = recode(pattern, !!!pat_labels),
    pattern = factor(pattern, levels = pat_label_levels)
  )
# --- HEATMAP ---
ggplot(heatmap_df, aes(x = factor(K3), y = pattern, fill = Max_Re_Eig)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = formatC(Max_Re_Eig, digits = 2)), size = 3.5, colour = "black")+
  scale_fill_gradient2(
    low      = "steelblue",
    mid      = "white",
    high     = "firebrick",
    midpoint = 0,
    name     = "Max Re(λ)"
  ) +
  scale_y_discrete(limits = rev(pat_label_levels)) +
  labs(
    title    = "Stability Heatmap — Maximum Real Eigenvalue at Origin During NTG",
    subtitle = "k1_n=0.75, k1_s=0.25, k1_i=0.5 | k2=0, D=0, K=1",
    x        = "Interconversion Rate (K3)",
    y        = "Pattern"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text=element_text(size=12),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold"),
    legend.position = "right"
  )

########## N pattern S inhibitor Origin ###############
# --- LABEL MAPPING ---
pat_labels <- c(
  "pattern_1" = "P1",
  "pattern_2" = "P2",
  "pattern_3" = "P3",
  "pattern_4" = "P4",
  "pattern_5" = "P5"
)

########## N pattern S inhibitor Origin ###############
# --- LOAD DATA ---
SDrugNPat_Origin <- read.csv("stability_NpatSinhibFULL.csv")

extract_real <- function(x) {
  as.numeric(gsub("\\s*[+\\-]\\s*[0-9eE.]+\\*I.*", "", as.character(x)))
}

SDrugNPat_Origin <- SDrugNPat_Origin %>%
  mutate(
    eig1_real  = extract_real(eig_1),
    eig2_real  = extract_real(eig_2),
    eig3_real  = extract_real(eig_3),
    Max_Re_Eig = pmax(eig1_real, eig2_real, eig3_real, na.rm = TRUE)
  )

range(SDrugNPat_Origin$Max_Re_Eig, na.rm = TRUE)

pattern_order <- c('pattern_1', 'pattern_2')
SDrugNPat_Origin <- SDrugNPat_Origin %>%
  mutate(
    Pattern    = recode(pattern, !!!pat_labels),
    Pattern    = factor(Pattern, levels = c("P1", "P2")),
    D          = factor(C),
    K3_rate    = factor(K3),
    Max_Re_Eig = as.numeric(Max_Re_Eig)
  )

k2_combos <- SDrugNPat_Origin %>%
  select(K2I, K2N) %>%
  distinct()

for (i in 1:nrow(k2_combos)) {
  
  k2i_val <- k2_combos$K2I[i]
  k2n_val <- k2_combos$K2N[i]
  
  plot_data <- SDrugNPat_Origin %>%
    filter(K2I == k2i_val, K2N == k2n_val)
  
  p <- ggplot(plot_data, aes(x = K3_rate, y = D, fill = Max_Re_Eig)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = signif(Max_Re_Eig, 2)), size = 3, colour = "black") +
    facet_wrap(~ Pattern, ncol = 1) +
    scale_fill_gradient2(
      low      = "steelblue",
      mid      = "white",
      high     = "firebrick",
      midpoint = 0,
      name     = "Max Re(λ)"
    ) +
    labs(
      title    = "N-optimal Patterns, S Inhibitor — Origin",
      subtitle = paste0("K2_N=", k2n_val, ", K2_I=", k2i_val, " | K2_S=1"),
      x        = "Interconversion Rate (K3)",
      y        = "Drug Concentration (D)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text=element_text(size=12),
      panel.grid       = element_blank(),
      strip.text       = element_text(face = "bold", size = 11),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      plot.title       = element_text(face = "bold"),
      legend.position  = "right",
      plot.background  = element_rect(fill = "white", colour = "white"),
      panel.background = element_rect(fill = "white", colour = "white")
    )
  
  filename <- paste0('NpatSinhibHM_K2I', k2i_val, '_K2N', k2n_val, '.png')
  ggsave(filename, plot = p, width = 16, height = 12)
  print(p)
  print(paste('Saved:', filename))
  save_plot(paste0("NpatSinhibHM_K2I", k2i_val, "_K2N", k2n_val, ".pdf"), p)
  
}


########## S pattern N inhibitor Origin ###############
# --- LOAD DATA ---
NDrugSPat_Origin <- read.csv("stability_SpatNinhibFULL.csv")

NDrugSPat_Origin <- NDrugSPat_Origin %>%
  mutate(
    eig1_real  = extract_real(eig_1),
    eig2_real  = extract_real(eig_2),
    eig3_real  = extract_real(eig_3),
    Max_Re_Eig = pmax(eig1_real, eig2_real, eig3_real, na.rm = TRUE)
  )

range(NDrugSPat_Origin$Max_Re_Eig, na.rm = TRUE)

pattern_order <- c('pattern_3', 'pattern_4', 'pattern_5')
NDrugSPat_Origin <- NDrugSPat_Origin %>%
  mutate(
    Pattern    = recode(pattern, !!!pat_labels),
    Pattern    = factor(Pattern, levels = c("P3", "P4", "P5")),
    D          = factor(C),
    K3_rate    = factor(K3),
    Max_Re_Eig = as.numeric(Max_Re_Eig)
  )

k2_combos <- NDrugSPat_Origin %>%
  select(K2I, K2S) %>%
  distinct()

for (i in 1:nrow(k2_combos)) {
  
  k2i_val <- k2_combos$K2I[i]
  k2s_val <- k2_combos$K2S[i]
  
  plot_data <- NDrugSPat_Origin %>%
    filter(K2I == k2i_val, K2S == k2s_val)
  
  p <- ggplot(plot_data, aes(x = K3_rate, y = D, fill = Max_Re_Eig)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = signif(Max_Re_Eig, 2)), size = 3, colour = "black") +
    facet_wrap(~ Pattern, ncol = 1) +
    scale_fill_gradient2(
      low      = "steelblue",
      mid      = "white",
      high     = "firebrick",
      midpoint = 0,
      name     = "Max Re(λ)"
    ) +
    labs(
      title    = "S-optimal Patterns, N Inhibitor — Origin",
      subtitle = paste0("K2_S=", k2s_val, ", K2_I=", k2i_val, " | K2_N=1"),
      x        = "Interconversion Rate (K3)",
      y        = "Drug Concentration (D)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text=element_text(size=12),
      panel.grid       = element_blank(),
      strip.text       = element_text(face = "bold", size = 11),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      plot.title       = element_text(face = "bold"),
      legend.position  = "right",
      plot.background  = element_rect(fill = "white", colour = "white"),
      panel.background = element_rect(fill = "white", colour = "white")
    )
  
  filename <- paste0('SpatNinhibHM_K2I', k2i_val, '_K2S', k2s_val, '.png')
  ggsave(filename, plot = p, width = 16, height = 12)
  print(p)
  save_plot(paste0("SpatNinhibHM_K2I", k2i_val, "_K2S", k2s_val, ".pdf"), p)
  print(paste('Saved:', filename))
}


###############################################
# Stability Analysis
###############################################



# --- LOAD DATA ---
#setwd("C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Stability/NoteBookOut")
setwd("C:/Users/Student/OneDrive/Documents/Year_3/Experiments/Stability/NoteBookOut")
NTGStab <- read.csv("stability_NTGFULL.csv")

# --- EXTRACT MAX REAL EIGENVALUE PER ROW ---
NTGStab <- NTGStab %>%
  mutate(
    Max_Re_Eig = pmax(eig_1, eig_2, eig_3)
  )

# --- EXTRACT K3 RATE ---
# All active k3 values are the same within a pattern so just take the max
#df <- df %>%
# mutate(
#  K3_rate = pmax(K3NS, K3NI, K3IN, K3IS, K3SI, K3SN)
#)

# --- KEEP ONE ROW PER PATTERN + K3 RATE ---
heatmap_NTGStab <- NTGStab %>%
  select(pattern, K3, Max_Re_Eig) %>%
  distinct()

# --- SET PATTERN ORDER ---

pattern_order <- c( "pattern_1", "pattern_2", "pattern_3", "pattern_4", "pattern_5")

# Define label mapping
pat_labels <- c(
  "pattern_1" = "P1",
  "pattern_2" = "P2",
  "pattern_3" = "P3",
  "pattern_4" = "P4",
  "pattern_5" = "P5"
)

PatInt <- c("pattern_1", "pattern_2", "pattern_3", "pattern_4", "pattern_5")
pat_label_levels <- unname(pat_labels)

heatmap_df <- heatmap_NTGStab %>%
  mutate(
    pattern = recode(pattern, !!!pat_labels),
    pattern = factor(pattern, levels = pat_label_levels)
  )
# --- HEATMAP ---
ggplot(heatmap_df, aes(x = factor(K3), y = pattern, fill = Max_Re_Eig)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = formatC(Max_Re_Eig, digits = 2)), size = 3.5, colour = "black")+
  scale_fill_gradient2(
    low      = "steelblue",
    mid      = "white",
    high     = "firebrick",
    midpoint = 0,
    name     = "Max Re(λ)"
  ) +
  scale_y_discrete(limits = rev(pat_label_levels)) +
  labs(
    title    = "Stability Heatmap — Maximum Real Eigenvalue at Origin During NTG",
    subtitle = "k1_n=0.75, k1_s=0.25, k1_i=0.5 | k2=0, D=0, K=1",
    x        = "Interconversion Rate (K3)",
    y        = "Pattern"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text=element_text(size=12),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold"),
    legend.position = "right"
  )


###############################################
# Validation Analysis
###############################################
###############################################
# Validation
###############################################
# -------------------------- Load Functions --------------------------

load_ts_avg <- function(db_path, table_name, schedule_label, cycle_label) {
  con <- dbConnect(SQLite(), db_path)
  df <- dbGetQuery(con, paste0("
    SELECT dominant,pattern_id, cycle, phase, time,
           AVG(N) as mean_N,
           AVG(I) as mean_I,
           AVG(S) as mean_S,
           AVG(N + I + S) as mean_total
    FROM ", table_name, "
    GROUP BY dominant, cycle, phase, time
    ORDER BY dominant, cycle, phase, time
  "))
  dbDisconnect(con)
  df$schedule <- schedule_label
  df$n_cycles <- cycle_label
  return(df)
}

load_ts_by_pattern <- function(db_path, table_name, schedule_label, cycle_label) {
  con <- dbConnect(SQLite(), db_path)
  df <- dbGetQuery(con, paste0("
    SELECT dominant, pattern_id, cycle, phase, time,
           AVG(N) as mean_N,
           AVG(I) as mean_I,
           AVG(S) as mean_S,
           AVG(N + I + S) as mean_total
    FROM ", table_name, "
    GROUP BY dominant, pattern_id, cycle, phase, time
    ORDER BY dominant, pattern_id, cycle, phase, time
  "))
  dbDisconnect(con)
  df$schedule <- schedule_label
  df$n_cycles <- cycle_label
  return(df)
}

# -------------------------- Load — averaged --------------------------
#Three phase
ts_avg_3N <- load_ts_avg(
  "C:/users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/MultiSimpleN-Inhib/timeseries.db",
  #"C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/MultiSimpleN-Inhib/timeseries.db",
  "Simple_TimeSeries", "Three Phase", "N-Inhibitor"
)
#ts_avg_3N<-subset(ts_avg_3N, ts_avg_3N$pattern_id != "P5")


ts_avg_3S <- load_ts_avg(
  "C:/users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/MultiSimpleS-Inhib/timeseries.db",
  #"C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/MultiSimpleS-Inhib/timeseries.db",
  "Simple_TimeSeries", "Three Phase", "S-Inhibitor"
)
#ts_avg_3S<-subset(ts_avg_3S, ts_avg_3S$pattern_id != "P5")

#Four phase
ts_avg_4c <- load_ts_avg(
  "C:/users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/Multi ComplexTS/timeseries.db",
  #"C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/Multi ComplexTS/timeseries.db",
  "Complex_TimeSeries", "Four Phase", "4-Cycle"
)
ts_avg_4N <- load_ts_avg(
  "C:/users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/Multi ComplexTS/timeseries.db",
  #"C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/Multi ComplexTS/timeseries.db",
  "Complex_TimeSeries", "Four Phase", "N-Inhibitor"
)
#ts_avg_4N<-subset(ts_avg_4N, ts_avg_4N$pattern_id != "P5")

ts_avg_4S <- load_ts_avg(
    "C:/users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/Multi ComplexTS/timeseries.db",
  #"C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/Multi ComplexTS/timeseries.db",
  "Complex_TimeSeries", "Four Phase", "S-Inhibitor"
)
#ts_avg_4S<-subset(ts_avg_4S, ts_avg_4S$pattern_id != "P5")

# -------------------------- Load — by pattern --------------------------
#Two phase
ts_pat_2phase <- load_ts_by_pattern(
    "C:/users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/P1End/timeseries.db",
  #"C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/P1End/timeseries.db",
  "Simple_TimeSeries", "Two Phase", "No-Inhibitor"
)
ts_pat_2phase<-subset(ts_pat_2phase, ts_pat_2phase$pattern_id != "P5")
#Three phase
ts_pat_3N <- load_ts_by_pattern(
    "C:/users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/MultiSimpleN-Inhib/timeseries.db",
  #"C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/MultiSimpleN-Inhib/timeseries.db",
  "Simple_TimeSeries", "Three Phase", "N-Inhibitor"
)
ts_pat_3N<-subset(ts_pat_3N, ts_pat_3N$pattern_id != "P5")

ts_pat_3S <- load_ts_by_pattern(
  "C:/users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/MultiSimpleS-Inhib/timeseries.db",
  #"C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/MultiSimpleS-Inhib/timeseries.db",
  "Simple_TimeSeries", "Three Phase", "S-Inhibitor"
)
ts_pat_3S<-subset(ts_pat_3S, ts_pat_3S$pattern_id != "P5")
#Four phase
ts_pat_4N <- load_ts_by_pattern(
    "C:/users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/Multi ComplexTS/timeseries.db",
  #"C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/Multi ComplexTS/timeseries.db",
  "Complex_TimeSeries", "Four Phase", "N-Inhibitor"
)
ts_pat_4N<-subset(ts_pat_4N, ts_pat_4N$pattern_id != "P5")

ts_pat_4S <- load_ts_by_pattern(
  "C:/users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/Multi ComplexTS/timeseries.db",
  #"C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/Multi ComplexTS/timeseries.db",
  "Complex_TimeSeries", "Four Phase", "S-Inhibitor"
)
ts_pat_4S<-subset(ts_pat_4S, ts_pat_4S$pattern_id != "P5")

ts_pat_4c <- load_ts_by_pattern(
    "C:/users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/Multi ComplexTS/timeseries.db",
  #"C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/Multi ComplexTS/timeseries.db",
  "Complex_TimeSeries", "Four Phase", "4-Cycle"
)

# -------------------------- Combine --------------------------

ts_avg <- bind_rows(ts_avg_3N, ts_avg_3S, ts_avg_4N, ts_avg_4S) %>%#ts_avg_4c
  mutate(
    schedule = factor(schedule, levels = c("Three Phase", "Four Phase")),
    n_cycles = factor(n_cycles, levels = c("N-Inhibitor", "S-Inhibitor", "4-Cycle")),
    phase    = factor(phase, levels = c(1, 2, 3, 4),
                      labels = c("Phase 1", "Phase 2", "Phase 3", "Phase 4"))
  )

ts_pat <- bind_rows(ts_pat_2phase, ts_pat_3N, ts_pat_3S, ts_pat_4N, ts_pat_4S) %>% #ts_pat_4c
  mutate(
    schedule   = factor(schedule,   levels = c("Two Phase","Three Phase", "Four Phase")),
    n_cycles   = factor(n_cycles,   levels = c("No-Inhibitor","N-Inhibitor", "S-Inhibitor", "4-Cycle")),
    pattern_id = factor(pattern_id, levels = c("P3", "P4")),
    phase      = factor(phase, levels = c(1, 2, 3, 4),
                        labels = c("Phase 1", "Phase 2", "Phase 3", "Phase 4"))
  )

# -------------------------- Helper: pivot to long --------------------------

pivot_ts <- function(df) {
  df %>%
    pivot_longer(cols = c(mean_N, mean_I, mean_S),
                 names_to  = "cell_type",
                 values_to = "mean_pop") %>%
    mutate(
      cell_type = recode(cell_type,
                         mean_N = "N Cell",
                         mean_I = "I Cell",
                         mean_S = "S Cell"),
      cell_type = factor(cell_type, levels = c("N Cell", "I Cell", "S Cell")),
      phase_fct = factor(paste0("Phase ", phase),
                         levels = paste0("Phase ", sort(unique(phase)))),
      cycle_fct = factor(paste0("Cycle ", cycle),
                         levels = paste0("Cycle ", sort(unique(cycle))))
    )
}

# -------------------------- Shared theme --------------------------

ts_theme <- list(
  scale_colour_manual(
    name   = "Initial Dominant",
    values = c("I Dominant" = Icol, "N Dominant" = Ncol,
               "No Dominant" = nocol, "S Dominant" = Scol)
  ),
  scale_linetype_manual(
    name   = "Cell Type",
    values = c("N Cell" = "solid", "I Cell" = "dashed", "S Cell" = "dotted")
  ),
  theme_classic(),
  theme(axis.text=element_text(size=12),legend.position = "bottom",
        strip.text      = element_text(size = 9, face = "bold"),
        panel.spacing   = unit(0.5, "lines")),
  guides(colour   = guide_legend(ncol = 1),
         linetype = guide_legend(ncol = 1))
)

# -------------------------- Plot functions --------------------------

make_ts_plot_all <- function(df, schedule_label, cycle_label) {
  df %>%
    filter(schedule == schedule_label, n_cycles == cycle_label) %>%
    pivot_ts() %>%
    ggplot(aes(x = time, y = mean_pop, colour = dominant, linetype = cell_type)) +
    geom_line(linewidth = 0.8) +
    facet_grid(cycle_fct ~ phase_fct, scales = "free_x") +
    ts_theme +
    labs(title    = paste(cycle_label, "|", schedule_label, "| All Patterns"),
         x        = "Timestep",
         y        = "Mean Population",
         colour   = "Initial Dominant",
         linetype = "Cell Type")
}

make_ts_plot_P4 <- function(df, schedule_label, cycle_label) {
  df %>%
    filter(schedule == schedule_label, n_cycles == cycle_label, pattern_id == "P4") %>%
    pivot_ts() %>%
    ggplot(aes(x = time, y = mean_pop, colour = dominant, linetype = cell_type)) +
    geom_line(linewidth = 0.8) +
    facet_grid(cycle_fct ~ phase_fct, scales = "free_x") +
    ts_theme +
    labs(title    = paste(cycle_label, "|", schedule_label, "| P4"),
         x        = "Timestep",
         y        = "Mean Population",
         colour   = "Initial Dominant",
         linetype = "Cell Type")
}
make_ts_plot_dom_inhib <- function(df, dominant_label, inhib_label) {
  dom_col <- case_when(
    dominant_label == "I Dominant"  ~ Icol,
    dominant_label == "N Dominant"  ~ Ncol,
    dominant_label == "No Dominant" ~ nocol,
    dominant_label == "S Dominant"  ~ Scol
  )
  df %>%
    filter(dominant == dominant_label) %>%
    group_by(dominant, schedule, n_cycles, cycle, phase, time) %>%
    summarise(
      mean_N     = mean(mean_N),
      mean_I     = mean(mean_I),
      mean_S     = mean(mean_S),
      mean_total = mean(mean_total),
      .groups = "drop"
    ) %>%
    pivot_ts() %>%
    ggplot(aes(x = time, y = mean_pop, linetype = cell_type)) +
    geom_line(linewidth = 0.8, colour = dom_col) +
    scale_linetype_manual(
      name   = "Cell Type",
      values = c("N Cell" = "solid", "I Cell" = "dashed", "S Cell" = "dotted")
    ) +
    facet_grid(cycle_fct ~ phase_fct, scales = "free_x") +
    theme_classic() +
    theme(axis.text=element_text(size=12),
          legend.position = "bottom",
          strip.text      = element_text(size = 9, face = "bold"),
          panel.spacing   = unit(0.4, "lines")) +
    guides(linetype = guide_legend(nrow = 1)) +
    labs(title    = paste(dominant_label, "|", inhib_label, "| All Patterns"),
         x        = "Timestep",
         y        = "Mean Population",
         linetype = "Cell Type")
}

# -------------------------- Build plots --------------------------

all_3N <- make_ts_plot_all(ts_avg_3N, "Three Phase", "N-Inhibitor")
all_3S <- make_ts_plot_all(ts_avg_3S, "Three Phase", "S-Inhibitor")
all_4N <- make_ts_plot_all(ts_avg_4N, "Four Phase",  "N-Inhibitor")
all_4S <- make_ts_plot_all(ts_avg_4S, "Four Phase",  "S-Inhibitor")
all_4c <- make_ts_plot_all(ts_avg_4c, "Four Phase",  "4-Cycle")

p4_3N <- make_ts_plot_P4(ts_pat_3N, "Three Phase", "N-Inhibitor")
p4_3S <- make_ts_plot_P4(ts_pat_3S, "Three Phase", "S-Inhibitor")
p4_4N <- make_ts_plot_P4(ts_pat_4N, "Four Phase",  "N-Inhibitor")
p4_4S <- make_ts_plot_P4(ts_pat_4S, "Four Phase",  "S-Inhibitor")
p4_4c <- make_ts_plot_P4(ts_pat_4c, "Four Phase",  "4-Cycle")

dom2_I  <- make_ts_plot_dom_inhib(ts_pat_2phase, "I Dominant",  "No-Inhibitor")
dom2_N  <- make_ts_plot_dom_inhib(ts_pat_2phase, "N Dominant",  "No-Inhibitor")
dom2_S  <- make_ts_plot_dom_inhib(ts_pat_2phase, "S Dominant",  "No-Inhibitor")
dom2_No  <- make_ts_plot_dom_inhib(ts_pat_2phase, "No Dominant",  "No-Inhibitor")

dom3_I_N  <- make_ts_plot_dom_inhib(ts_pat_3N, "I Dominant",  "N-Inhibitor")
dom3_N_N  <- make_ts_plot_dom_inhib(ts_pat_3N, "N Dominant",  "N-Inhibitor")
dom3_S_N  <- make_ts_plot_dom_inhib(ts_pat_3N, "S Dominant",  "N-Inhibitor")
dom3_No_N <- make_ts_plot_dom_inhib(ts_pat_3N, "No Dominant", "N-Inhibitor")

dom4_I_N  <- make_ts_plot_dom_inhib(ts_pat_4N, "I Dominant",  "N-Inhibitor and S-Inhibitor")
dom4_N_N  <- make_ts_plot_dom_inhib(ts_pat_4N, "N Dominant",  "N-Inhibitor and S-Inhibitor")
dom4_S_N  <- make_ts_plot_dom_inhib(ts_pat_4N, "S Dominant",  "N-Inhibitor and S-Inhibitor")
dom4_No_N <- make_ts_plot_dom_inhib(ts_pat_4N, "No Dominant", "N-Inhibitor and S-Inhibitor")

dom3_I_S  <- make_ts_plot_dom_inhib(ts_pat_3S, "I Dominant",  "S-Inhibitor")
dom3_N_S  <- make_ts_plot_dom_inhib(ts_pat_3S, "N Dominant",  "S-Inhibitor")
dom3_S_S  <- make_ts_plot_dom_inhib(ts_pat_3S, "S Dominant",  "S-Inhibitor")
dom3_No_S <- make_ts_plot_dom_inhib(ts_pat_3S, "No Dominant", "S-Inhibitor")

dom4_I_S  <- make_ts_plot_dom_inhib(ts_pat_4S, "I Dominant",  "S-Inhibitor and N-Inhibitor")
dom4_N_S  <- make_ts_plot_dom_inhib(ts_pat_4S, "N Dominant",  "S-Inhibitor and N-Inhibitor")
dom4_S_S  <- make_ts_plot_dom_inhib(ts_pat_4S, "S Dominant",  "S-Inhibitor and N-Inhibitor")
dom4_No_S <- make_ts_plot_dom_inhib(ts_pat_4S, "No Dominant", "S-Inhibitor and N-Inhibitor")
# -------------------------- Patchwork --------------------------

ThreePhase_ts_all <-
  (all_3N / all_3S) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase Schedule: Mean Time Series — All Patterns",
    subtitle = "Top = N Inhibitor | Bottom = S Inhibitor\nSolid = N Cell | Dashed = I Cell | Dotted = S Cell",
    theme    = theme(axis.text=element_text(size=12),plot.title = element_text(size = 14, face = "bold"))
  ) & theme(legend.position = "bottom")

FourPhase_ts_all <-
  (all_4N / all_4S) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Four Phase Schedule: Mean Time Series — All Patterns",
    subtitle = "Top = N Inhibitor | Bottom = S Inhibitor\nSolid = N Cell | Dashed = I Cell | Dotted = S Cell",
    theme    = theme(axis.text=element_text(size=12),plot.title = element_text(size = 14, face = "bold"))
  ) & theme(legend.position = "bottom")

ThreePhase_ts_P4 <-
  (p4_3N / p4_3S) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase Schedule: Mean Time Series — P4",
    subtitle = "Top = N Inhibitor | Bottom = S Inhibitor\nSolid = N Cell | Dashed = I Cell | Dotted = S Cell",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) & theme(legend.position = "bottom")

FourPhase_ts_P4 <-
  p4_4c +
  plot_annotation(
    title    = "Four Phase Schedule: Mean Time Series — P4",
    subtitle = "Solid = N Cell | Dashed = I Cell | Dotted = S Cell",
    theme    = theme(axis.text=element_text(size=12),plot.title = element_text(size = 14, face = "bold"))
  ) & theme(legend.position = "bottom")

all_dom2 <-
  (dom2_No / dom2_N / dom2_I / dom2_S) +
  plot_annotation(
    title    = "Two Phase Schedule — No Inhibitor: Mean Time Series by Initial Dominant",
    subtitle = "Solid = N Cell | Dashed = I Cell | Dotted = S Cell",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) & theme(legend.position = "bottom")

all_dom3_N <-
  (dom3_No_N / dom3_N_N / dom3_I_N / dom3_S_N) +
  plot_annotation(
    title    = "Three Phase Schedule — N Inhibitor: Mean Time Series by Initial Dominant",
    subtitle = "Solid = N Cell | Dashed = I Cell | Dotted = S Cell",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) & theme(legend.position = "bottom")

all_dom4_N <-
  (dom4_No_N / dom4_N_N / dom4_I_N / dom4_S_N) +
  plot_annotation(
    title    = "Four Phase Schedule — Two Inhibitor: Mean Time Series by Initial Dominant",
    subtitle = "Solid = N Cell | Dashed = I Cell | Dotted = S Cell",
    theme    = theme(axis.text=element_text(size=12),plot.title = element_text(size = 14, face = "bold"))
  ) & theme(legend.position = "bottom")

all_dom3_S <-
  (dom3_No_S / dom3_N_S / dom3_I_S / dom3_S_S) +
  plot_annotation(
    title    = "Three Phase Schedule — S Inhibitor: Mean Time Series by Initial Dominant",
    subtitle = "Solid = N Cell | Dashed = I Cell | Dotted = S Cell",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) & theme(legend.position = "bottom")

all_dom4_S <-
  (dom4_No_S / dom4_N_S / dom4_I_S / dom4_S_S) +
  plot_annotation(
    title    = "Four Phase Schedule — S Inhibitor: Mean Time Series by Initial Dominant",
    subtitle = "Solid = N Cell | Dashed = I Cell | Dotted = S Cell",
    theme    = theme(axis.text=element_text(size=12),plot.title = element_text(size = 14, face = "bold"))
  ) & theme(legend.position = "bottom")

# -------------------------- Save all plots --------------------------
save_plot("FourPhase_TimeToRemission.pdf",complex_drug_phases_4C)
save_plot("FourPhase_TimeToRelapse.pdf",complex_relapse_phases_4C)
save_plot("ThreePhase_TimeToRemission.pdf", simple_drug_phases)
save_plot("ThreePhase_TimeToRelapse.pdf", simple_relapse_phases)
save_plot("ts_ThreePhase_all.pdf",       ThreePhase_ts_all)
save_plot("ts_FourPhase_all.pdf",        FourPhase_ts_all)
save_plot("ts_ThreePhase_P4.pdf",        ThreePhase_ts_P4)
save_plot("ts_FourPhase_P4.pdf",         FourPhase_ts_P4)
save_plot("ts_ThreePhase_all_N.pdf",     all_3N)
save_plot("ts_ThreePhase_all_S.pdf",     all_3S)

save_plot("ts_FourPhase_all_N.pdf",      all_4N)
save_plot("ts_FourPhase_all_S.pdf",      all_4S)
save_plot("ts_FourPhase_all_4c.pdf",     all_4c)
save_plot("ts_ThreePhase_P4_N.pdf",      p4_3N)
save_plot("ts_ThreePhase_P4_S.pdf",      p4_3S)
save_plot("ts_FourPhase_P4_N.pdf",       p4_4N)
save_plot("ts_FourPhase_P4_S.pdf",       p4_4S)
save_plot("ts_FourPhase_P4_4c.pdf",      p4_4c)


save_plot("ts_dom2_all.pdf",           all_dom2)
save_plot("ts_dom3_N_all.pdf",           all_dom3_N)
save_plot("ts_dom4_N_all.pdf",           all_dom4_N)
save_plot("ts_dom3_S_all.pdf",           all_dom3_S)
save_plot("ts_dom4_S_all.pdf",           all_dom4_S)

save_plot("ts_2_Idom.pdf",              dom2_I)
save_plot("ts_2_Ndom.pdf",              dom2_N)
save_plot("ts_2_Sdom.pdf",              dom2_S)
save_plot("ts_2_Nodom.pdf",             dom2_No)

save_plot("ts_3N_Idom.pdf",              dom3_I_N)
save_plot("ts_3N_Ndom.pdf",              dom3_N_N)
save_plot("ts_3N_Sdom.pdf",              dom3_S_N)
save_plot("ts_3N_Nodom.pdf",             dom3_No_N)
save_plot("ts_3S_Idom.pdf",              dom3_I_S)
save_plot("ts_3S_Ndom.pdf",              dom3_N_S)
save_plot("ts_3S_Sdom.pdf",              dom3_S_S)
save_plot("ts_3S_Nodom.pdf",             dom3_No_S)

save_plot("ts_4N_Idom.pdf",              dom4_I_N)
save_plot("ts_4N_Ndom.pdf",              dom4_N_N)
save_plot("ts_4N_Sdom.pdf",              dom4_S_N)
save_plot("ts_4N_Nodom.pdf",             dom4_No_N)

save_plot("ts_4_Idom.pdf",              dom4_I_S)
save_plot("ts_4S_Ndom.pdf",              dom4_N_S)
save_plot("ts_4S_Sdom.pdf",              dom4_S_S)
save_plot("ts_4S_Nodom.pdf",             dom4_No_S)

