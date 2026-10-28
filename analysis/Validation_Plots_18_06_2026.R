# -------------------------- Validation Plots --------------------------
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)

# -------------------------- Load Data --------------------------

#setwd("C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/3 cycle/")
#Multi3Simple  <- read.csv("Multi SimpleTS/SimpleResults.csv")
#Multi3Complex <- read.csv("Multi ComplexTS/ComplexResults.csv")

setwd("C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle/")
#setwd("C:/Users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/4 cycle")
Multi4SimpleNihib  <- read.csv("MultiSimpleN-Inhib/SimpleResults.csv")
Multi4SimpleSihib  <- read.csv("MultiSimpleS-Inhib/SimpleResults.csv")
Multi4Complex <- read.csv("Multi ComplexTS/ComplexResults.csv")

# -------------------------- Output Directory --------------------------

out_dir <- "C:/users/892215/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Validity/Plots/Figures/"
dir.create(out_dir, showWarnings = FALSE)

# -------------------------- Prep Functions --------------------------

prep_simple <- function(df) {
  df %>%
    mutate(
      Phase1_Relapse      = Phase1_Relapse == "True",
      Phase2_Reached_Zero = Phase2_Reached_Zero == "True",
      Phase3_Relapse      = Phase3_Relapse == "True",
      Phase1_Recovery_Time = if_else(Phase1_Relapse & is.na(Phase1_Recovery_Time), 30,  Phase1_Recovery_Time),
      Phase2_Zero_Time     = if_else(Phase2_Reached_Zero & is.na(Phase2_Zero_Time), 100, Phase2_Zero_Time),
      Phase3_Recovery_Time = if_else(Phase3_Relapse & is.na(Phase3_Recovery_Time), 30,  Phase3_Recovery_Time)
    )
}

prep_complex <- function(df) {
  df %>%
    mutate(
      Phase1_Relapse      = Phase1_Relapse == "True",
      Phase2_Reached_Zero = Phase2_Reached_Zero == "True",
      Phase3_Relapse      = Phase3_Relapse == "True",
      Phase4_Reached_Zero = Phase4_Reached_Zero == "True",
      Phase1_Recovery_Time = if_else(Phase1_Relapse & is.na(Phase1_Recovery_Time), 30,  Phase1_Recovery_Time),
      Phase2_Zero_Time     = if_else(Phase2_Reached_Zero & is.na(Phase2_Zero_Time), 100, Phase2_Zero_Time),
      Phase3_Recovery_Time = if_else(Phase3_Relapse & is.na(Phase3_Recovery_Time), 30,  Phase3_Recovery_Time),
      Phase4_Zero_Time     = if_else(Phase4_Reached_Zero & is.na(Phase4_Zero_Time), 100, Phase4_Zero_Time)
    )
}

# -------------------------- Apply Prep --------------------------

#Multi3Simple  <- prep_simple(Multi3Simple)
Multi4SimpleNihib  <- prep_simple(Multi4SimpleNihib)
Multi4SimpleSihib  <- prep_simple(Multi4SimpleSihib)
#Multi3Complex <- prep_complex(Multi3Complex)
Multi4Complex <- prep_complex(Multi4Complex)

# -------------------------- Crossing Time Helpers --------------------------

make_simple_crossing_plots <- function(df, label) {
  
  p1_data <- df %>% filter(!is.na(Phase1_Recovery_Time)) %>% mutate(Cycle = factor(Cycle))
  p2_data <- df %>% filter(!is.na(Phase2_Zero_Time))     %>% mutate(Cycle = factor(Cycle))
  p3_data <- df %>% filter(!is.na(Phase3_Recovery_Time)) %>% mutate(Cycle = factor(Cycle))
  
  base_theme <- list(
    geom_boxplot(alpha = 0.3, outlier.size = 1),
    facet_grid(Pattern.ID ~ .),
    theme_bw(),
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  )
  
  p1 <- ggplot(p1_data, aes(x = Cycle, y = Phase1_Recovery_Time, colour = Dominant)) +
    base_theme +
    labs(title    = paste(label, "| Phase 1: Rise Above 1e-3"),
         subtitle = paste0("n = ", nrow(p1_data), " runs"),
         x = "Cycle", y = "Timestep")
  
  p2 <- ggplot(p2_data, aes(x = Cycle, y = Phase2_Zero_Time, colour = Dominant)) +
    base_theme +
    labs(title    = paste(label, "| Phase 2: Drop Below 1e-3"),
         subtitle = paste0("n = ", nrow(p2_data), " runs"),
         x = "Cycle", y = "Timestep")
  
  p3 <- ggplot(p3_data, aes(x = Cycle, y = Phase3_Recovery_Time, colour = Dominant)) +
    base_theme +
    labs(title    = paste(label, "| Phase 3: Rise Above 1e-3"),
         subtitle = paste0("n = ", nrow(p3_data), " runs"),
         x = "Cycle", y = "Timestep") +
    theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))
  
  list(p1 = p1, p2 = p2, p3 = p3)
}

make_complex_crossing_plots <- function(df, label) {
  
  p1_data <- df %>% filter(!is.na(Phase1_Recovery_Time)) %>% mutate(Cycle = factor(Cycle))
  p2_data <- df %>% filter(!is.na(Phase2_Zero_Time))     %>% mutate(Cycle = factor(Cycle))
  p3_data <- df %>% filter(!is.na(Phase3_Recovery_Time)) %>% mutate(Cycle = factor(Cycle))
  p4_data <- df %>% filter(!is.na(Phase4_Zero_Time))     %>% mutate(Cycle = factor(Cycle))
  
  base_theme <- list(
    geom_boxplot(alpha = 0.3, outlier.size = 1),
    facet_grid(Pattern.ID ~ .),
    theme_bw(),
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  )
  
  p1 <- ggplot(p1_data, aes(x = Cycle, y = Phase1_Recovery_Time, colour = Dominant)) +
    base_theme +
    labs(title    = paste(label, "| Phase 1: Rise Above 1e-3"),
         subtitle = paste0("n = ", nrow(p1_data), " runs"),
         x = "Cycle", y = "Timestep")
  
  p2 <- ggplot(p2_data, aes(x = Cycle, y = Phase2_Zero_Time, colour = Dominant)) +
    base_theme +
    labs(title    = paste(label, "| Phase 2: Drop Below 1e-3"),
         subtitle = paste0("n = ", nrow(p2_data), " runs"),
         x = "Cycle", y = "Timestep")
  
  p3 <- ggplot(p3_data, aes(x = Cycle, y = Phase3_Recovery_Time, colour = Dominant)) +
    base_theme +
    labs(title    = paste(label, "| Phase 3: Rise Above 1e-3"),
         subtitle = paste0("n = ", nrow(p3_data), " runs"),
         x = "Cycle", y = "Timestep")
  
  p4 <- ggplot(p4_data, aes(x = Cycle, y = Phase4_Zero_Time, colour = Dominant)) +
    base_theme +
    labs(title    = paste(label, "| Phase 4: Drop Below 1e-3"),
         subtitle = paste0("n = ", nrow(p4_data), " runs"),
         x = "Cycle", y = "Timestep") +
    theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))
  
  
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
}

# -------------------------- Run Count Helpers --------------------------

make_simple_count_plot <- function(df, label) {
  df %>%
    group_by(Pattern.ID, Dominant, Cycle) %>%
    summarise(
      Phase1_Relapsed = sum(Phase1_Relapse      == TRUE, na.rm = TRUE),
      Phase2_Crossed  = sum(Phase2_Reached_Zero == TRUE, na.rm = TRUE),
      Phase3_Relapsed = sum(Phase3_Relapse      == TRUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(Phase1_Relapsed, Phase2_Crossed, Phase3_Relapsed),
                 names_to = "Event", values_to = "Count") %>%
    mutate(Event = case_when(
      Event == "Phase1_Relapsed" ~ "Phase 1 Above 1e-3",
      Event == "Phase2_Crossed"  ~ "Phase 2 Below 1e-3",
      Event == "Phase3_Relapsed" ~ "Phase 3 Above 1e-3"
    ),
    Event = factor(Event, levels = c(
      "Phase 1 Above 1e-3", "Phase 2 Below 1e-3", "Phase 3 Above 1e-3"
    ))) %>%
    ggplot(aes(x = factor(Cycle), y = Count, fill = Dominant)) +
    geom_col(position = "dodge") +
    facet_grid(Pattern.ID ~ Event) +
    labs(title = paste(label, "| Run Counts"),
         x = "Cycle", y = "Number of runs") +
    theme_bw() +
    theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))
  
}

make_complex_count_plot <- function(df, label) {
  df %>%
    group_by(Pattern.ID, Dominant, Cycle) %>%
    summarise(
      Phase1_Relapsed = sum(Phase1_Relapse      == TRUE, na.rm = TRUE),
      Phase2_Crossed  = sum(Phase2_Reached_Zero == TRUE, na.rm = TRUE),
      Phase3_Relapsed = sum(Phase3_Relapse      == TRUE, na.rm = TRUE),
      Phase4_Crossed  = sum(Phase4_Reached_Zero == TRUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(Phase1_Relapsed, Phase2_Crossed,
                          Phase3_Relapsed, Phase4_Crossed),
                 names_to = "Event", values_to = "Count") %>%
    mutate(Event = case_when(
      Event == "Phase1_Relapsed" ~ "Phase 1 Above 1e-3",
      Event == "Phase2_Crossed"  ~ "Phase 2 Below 1e-3",
      Event == "Phase3_Relapsed" ~ "Phase 3 Above 1e-3",
      Event == "Phase4_Crossed"  ~ "Phase 4 Below 1e-3"
    ),
    Event = factor(Event, levels = c(
      "Phase 1 Above 1e-3", "Phase 2 Below 1e-3",
      "Phase 3 Above 1e-3", "Phase 4 Below 1e-3"
    ))) %>%
    ggplot(aes(x = factor(Cycle), y = Count, fill = Dominant)) +
    geom_col(position = "dodge") +
    facet_grid(Pattern.ID ~ Event) +
    labs(title = paste(label, "| Run Counts"),
         x = "Cycle", y = "Number of runs") +
    theme_bw() +
    theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))
  
}

# -------------------------- Appendix Heatmap Helpers (A and B) --------------------------

make_heatmap_A <- function(df, label) {
  df %>%
    group_by(Dominant, Pattern.ID, Cycle) %>%
    summarise(
      Phase1 = mean(Phase1_Total, na.rm = TRUE),
      Phase2 = mean(Phase2_Total, na.rm = TRUE),
      Phase3 = mean(Phase3_Total, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(Phase1, Phase2, Phase3),
                 names_to = "Phase", values_to = "Mean_Pop") %>%
    mutate(
      Phase      = factor(Phase, levels = c("Phase1", "Phase2", "Phase3"),
                          labels = c("Phase 1", "Phase 2", "Phase 3")),
      Cycle      = factor(Cycle),
      label_text = ifelse(Mean_Pop < 0.001, "~0",
                          formatC(Mean_Pop, format = "f", digits = 3))
    ) %>%
    ggplot(aes(x = Cycle, y = Dominant, fill = Mean_Pop)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = label_text), size = 3) +
    facet_grid(Pattern.ID ~ Phase) +
    scale_fill_gradient(low = "white", high = "red",
                        name = "Mean/nPopulation", limits = c(0, 1)) +
    labs(title = paste(label, "| Mean End-of-Phase Total Population"),
         x = "Cycle", y = "Initial Dominant") +
    theme_bw() +
    theme(legend.position = "right",
          strip.text      = element_text(size = 9, face = "bold"),
          axis.text       = element_text(size = 9),
          panel.spacing   = unit(0.5, "lines"))
}

make_heatmap_B <- function(df, label) {
  df %>%
    group_by(Dominant, Pattern.ID, Cycle) %>%
    summarise(
      Phase1 = mean(Phase1_Total >  1e-3, na.rm = TRUE),
      Phase2 = mean(Phase2_Total <= 1e-3, na.rm = TRUE),
      Phase3 = mean(Phase3_Total >  1e-3, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(Phase1, Phase2, Phase3),
                 names_to = "Phase", values_to = "Prop") %>%
    mutate(
      Phase      = factor(Phase,
                          levels = c("Phase1", "Phase2", "Phase3"),
                          labels = c("Phase 1/n(% above 1e-3)",
                                     "Phase 2/n(% below 1e-3)",
                                     "Phase 3/n(% above 1e-3)")),
      Cycle      = factor(Cycle),
      label_text = paste0(round(Prop * 100), "%")
    ) %>%
    ggplot(aes(x = Cycle, y = Dominant, fill = Prop)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = label_text), size = 3) +
    facet_grid(Pattern.ID ~ Phase) +
    scale_fill_gradient(low = "white", high = "red",
                        name = "Proportion",
                        limits = c(0, 1),
                        labels = scales::percent) +
    labs(title    = paste(label, "| Proportion of Runs Crossing Threshold at End of Phase"),
         subtitle = "Phase 1 & 3: % above 1e-3 | Phase 2: % below 1e-3",
         x = "Cycle", y = "Initial Dominant") +
    theme_bw() +
    theme(legend.position = "right",
          strip.text      = element_text(size = 9, face = "bold"),
          axis.text       = element_text(size = 9),
          panel.spacing   = unit(0.5, "lines"))
}
# -------------------------- Main Text Heatmap Helper (D) --------------------------


make_heatmap_D <- function(df, label, include_phase4 = FALSE) {
  
  base <- df %>%
    filter(Pattern.ID == "P4") %>%
    group_by(Dominant, Cycle) %>%
    summarise(
      Phase1_mean = mean(Phase1_Total,        na.rm = TRUE),
      Phase2_mean = mean(Phase2_Total,        na.rm = TRUE),
      Phase3_mean = mean(Phase3_Total,        na.rm = TRUE),
      Phase1_prop = mean(Phase1_Total >  1e-3, na.rm = TRUE),
      Phase2_prop = mean(Phase2_Total <= 1e-3, na.rm = TRUE),
      Phase3_prop = mean(Phase3_Total >  1e-3, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (include_phase4 && "Phase4_Total" %in% names(df)) {
    base <- df %>%
      filter(Pattern.ID == "P4") %>%
      group_by(Dominant, Cycle) %>%
      summarise(
        Phase1_mean = mean(Phase1_Total,        na.rm = TRUE),
        Phase2_mean = mean(Phase2_Total,        na.rm = TRUE),
        Phase3_mean = mean(Phase3_Total,        na.rm = TRUE),
        Phase4_mean = mean(Phase4_Total,        na.rm = TRUE),
        Phase1_prop = mean(Phase1_Total >  1e-3, na.rm = TRUE),
        Phase2_prop = mean(Phase2_Total <= 1e-3, na.rm = TRUE),
        Phase3_prop = mean(Phase3_Total >  1e-3, na.rm = TRUE),
        Phase4_prop = mean(Phase4_Total <= 1e-3, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  base %>%
    pivot_longer(cols = -c(Dominant, Cycle),
                 names_to = "variable", values_to = "Value") %>%
    mutate(
      Phase  = case_when(
        grepl("Phase1", variable) ~ "Phase 1",
        grepl("Phase2", variable) ~ "Phase 2",
        grepl("Phase3", variable) ~ "Phase 3",
        grepl("Phase4", variable) ~ "Phase 4"
      ),
      Metric = case_when(
        grepl("mean", variable) ~ "Mean Population",
        grepl("prop", variable) ~ "% Crossing Threshold"
      ),
      Phase      = factor(Phase,  levels = c("Phase 1", "Phase 2",
                                             "Phase 3", "Phase 4")),
      Metric     = factor(Metric, levels = c("Mean Population",
                                             "% Crossing Threshold")),
      Cycle      = factor(Cycle),
      label_text = case_when(
        Metric == "Mean Population" & Value < 0.001 ~ "~0",
        Metric == "Mean Population"                 ~ formatC(Value, format = "f", digits = 3),
        Metric == "% Crossing Threshold"            ~ paste0(round(Value * 100), "%")
      )
    ) %>%
    ggplot(aes(x = Cycle, y = Dominant, fill = Value)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = label_text), size = 3.5) +
    facet_grid(Metric ~ Phase, scales = "free") +
    scale_fill_gradient(low = "white", high = "red",
                        name = "Value", limits = c(0, 1)) +
    labs(title    = paste(label, "| P4: Dominant Type x Cycle"),
         subtitle = "Top = mean end-of-phase population | Bottom = % crossing threshold/n(Phase 1 & 3: % above 1e-3 | Phase 2 & 4: % below 1e-3)",
         x = "Cycle", y = "Initial Dominant") +
    theme_bw() +
    theme(legend.position = "right",
          strip.text      = element_text(size = 9, face = "bold"),
          axis.text       = element_text(size = 9),
          panel.spacing   = unit(0.6, "lines"))
}

# -------------------------- Build All Plots --------------------------

#s3 <- make_simple_crossing_plots(Multi3Simple,   "3-Cycle")
s4N <- make_simple_crossing_plots(Multi4SimpleNihib,   "4-Cycle")
s4S <- make_simple_crossing_plots(Multi4SimpleSihib,   "4-Cycle")
#c3 <- make_complex_crossing_plots(Multi3Complex, "3-Cycle")
c4 <- make_complex_crossing_plots(Multi4Complex, "4-Cycle")

#sc3 <- make_simple_count_plot(Multi3Simple,   "3-Cycle")
sc4N <- make_simple_count_plot(Multi4SimpleNihib,   "4-Cycle")
sc4S <- make_simple_count_plot(Multi4SimpleSihib,   "4-Cycle")
#cc3 <- make_complex_count_plot(Multi3Complex, "3-Cycle")
cc4 <- make_complex_count_plot(Multi4Complex, "4-Cycle")

# Appendix heatmaps — simple
#hA_s3 <- make_heatmap_A(Multi3Simple, "3-Cycle Simple")
hA_s4N <- make_heatmap_A(Multi4SimpleNihib, "4-Cycle Simple")
hA_s4S <- make_heatmap_A(Multi4SimpleSihib, "4-Cycle Simple")
#hB_s3 <- make_heatmap_B(Multi3Simple, "3-Cycle Simple")
hB_s4 <- make_heatmap_B(Multi4SimpleNihib, "4-Cycle Simple")

# Appendix heatmaps — complex
#hA_c3 <- make_heatmap_A(Multi3Complex, "3-Cycle Complex")
hA_c4 <- make_heatmap_A(Multi4Complex, "4-Cycle Complex")
#hB_c3 <- make_heatmap_B(Multi3Complex, "3-Cycle Complex")
hB_c4 <- make_heatmap_B(Multi4Complex, "4-Cycle Complex")

# Main text heatmaps — D
#hD_s3 <- make_heatmap_D(Multi3Simple,  "3-Cycle Simple",  include_phase4 = FALSE)
hD_s4N <- make_heatmap_D(Multi4SimpleNihib,  "4-Cycle Simple",  include_phase4 = FALSE)
hD_s4S <- make_heatmap_D(Multi4SimpleSihib,  "4-Cycle Simple",  include_phase4 = FALSE)
#hD_c3 <- make_heatmap_D(Multi3Complex, "3-Cycle Complex", include_phase4 = TRUE)
hD_c4 <- make_heatmap_D(Multi4Complex, "4-Cycle Complex", include_phase4 = TRUE)

# -------------------------- Patchwork: Crossing Times --------------------------

simple_drug_phases <-
  (s4N$p2 | s4S$p2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase Schedule: Time to Remission in Drug Phase",
    subtitle = "Left = N Inhibitor | Right = S Inhibitor",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))

simple_drug_phases_4CS <-
  (s4S$p2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase S Inhibitor Schedule: Time to Remission in Drug Phase",
    subtitle = "4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))

simple_drug_phases_4CN <-
 (s4N$p2) +
plot_layout(guides = "collect") +
  plot_annotation(
   title    = "Three Phase N Inhibitor Schedule: Time to Remission in Drug Phase",
  subtitle = "4-cycle",
 theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
 theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))


simple_relapse_phases <-
  (s4N$p1 | s4S$p1) /
  (s4N$p3 | s4S$p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase Schedule: Time to Relapse in Off-Drug Phases",
    subtitle = "Left = N Inhibitor | Right = S Inhibitor",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))

simple_relapse_phases_4CN <-
  ( s4N$p1) /
  ( s4N$p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase N Inhibitor Schedule: Time to Relapse in Off-Drug Phases",
    subtitle = "4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))

simple_relapse_phases_4CS <-
  ( s4S$p1) /
  ( s4S$p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase S Inhibitor Schedule: Time to Relapse in Off-Drug Phases",
    subtitle = "4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))

#complex_drug_phases <-
 # (c3$p2 | c4$p2) /
  #(c3$p4 | c4$p4) +
  #plot_layout(guides = "collect") +
  #plot_annotation(
#    title    = "Four Phase Schedule: Time to Remission in Drug Phases",
 #   subtitle = "Left = 3-cycle | Right = 4-cycle",
  #  theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  #) &
  #theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))


#complex_relapse_phases <-
 # (c3$p1 | c4$p1) /
  #(c3$p3 | c4$p3) +
#  plot_layout(guides = "collect") +
 # plot_annotation(
  #  title    = "Four Phase Schedule: Time to Relapse in Off-Drug Phases",
   # subtitle = "Left = 3-cycle | Right = 4-cycle",
    #theme    = theme(plot.title = element_text(size = 14, face = "bold"))
#  ) &
 # theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))

complex_drug_phases_4C <-
  ( c4$p2) /
  (c4$p4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Four Phase Schedule: Time to Remission in Drug Phases",
    subtitle = "4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))


complex_relapse_phases_4C <-
  (c4$p1) /
  (c4$p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Four Phase Schedule: Time to Relapse in Off-Drug Phases",
    subtitle = " 4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))



# -------------------------- Patchwork: Run Counts --------------------------

simple_counts <-
  (sc4N | sc4S) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase Schedule: Number of Runs Crossing Threshold",
    subtitle = "Left = N Inhibitor | Right = S Inhibitor",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))


complex_counts <-
  (cc3 | cc4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Four Phase Schedule: Number of Runs Crossing Threshold",
    subtitle = "Left = 3-cycle | Right = 4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))


simple_counts_4CN <-
  (sc4N) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase N Inhibitor Schedule: Number of Runs Crossing Threshold",
    subtitle = "4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))

simple_counts_4CS <-
  (sc4S) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase S Inhibitor Schedule: Number of Runs Crossing Threshold",
    subtitle = "4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))


complex_counts_4C <-
  (cc4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Four Phase Schedule: Number of Runs Crossing Threshold",
    subtitle = " 4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "bottom", axis.text=element_text(size=12), legend.title = element_text(size = 15), legend.text = element_text(size=14))


# -------------------------- Patchwork: Appendix Heatmaps (A and B) --------------------------

simple_heatmaps_A <-
  (hA_s3 / hA_s4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase Schedule: Mean End-of-Phase Population",
    subtitle = "Top = 3-cycle | Bottom = 4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "right")

simple_heatmaps_B <-
  (hB_s3 / hB_s4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase Schedule: Proportion of Runs Below 1e-3",
    subtitle = "Top = 3-cycle | Bottom = 4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "right")

complex_heatmaps_A <-
  (hA_c3 / hA_c4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Four Phase Schedule: Mean End-of-Phase Population",
    subtitle = "Top = 3-cycle | Bottom = 4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "right")

complex_heatmaps_B <-
  (hB_c3 / hB_c4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Four Phase Schedule: Proportion of Runs Below 1e-3",
    subtitle = "Top = 3-cycle | Bottom = 4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "right")

# -------------------------- Patchwork: Main Text Heatmaps (D) --------------------------

simple_heatmaps_D <-
  (hD_s4N / hD_s4S) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase Schedule: P4 Dominant Type x Cycle",
    subtitle = "Left = N Inhibitor | Right = S Inhibitor",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "right")

complex_heatmaps_D <-
  (hD_c3 / hD_c4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Four Phase Schedule: P4 Dominant Type x Cycle",
    subtitle = "Top = 3-cycle | Bottom = 4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "right")


simple_heatmaps_D_4CN <-
  (hD_s4N) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase N Inhibitor Schedule: P4 Dominant Type x Cycle",
    subtitle = "4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "right")

simple_heatmaps_D_4CS <-
  (hD_s4S) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three Phase S Inhibitor Schedule: P4 Dominant Type x Cycle",
    subtitle = "4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "right")

complex_heatmaps_D_4C <-
  (hD_c4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Four Phase Schedule: P4 Dominant Type x Cycle",
    subtitle = "4-cycle",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "right")


# -------------------------- Save All --------------------------

# Crossing times
ggsave(paste0(out_dir, "simple_drug_phases.png"),     simple_drug_phases,     width = 14, height = 6,  dpi = 300)
ggsave(paste0(out_dir, "simple_relapse_phases.png"),  simple_relapse_phases,  width = 14, height = 10, dpi = 300)
ggsave(paste0(out_dir, "complex_drug_phases.png"),    complex_drug_phases,    width = 14, height = 10, dpi = 300)
ggsave(paste0(out_dir, "complex_relapse_phases.png"), complex_relapse_phases, width = 14, height = 10, dpi = 300)

# Run counts
ggsave(paste0(out_dir, "simple_counts.png"),          simple_counts,          width = 14, height = 8,  dpi = 300)
ggsave(paste0(out_dir, "complex_counts.png"),         complex_counts,         width = 14, height = 8,  dpi = 300)

# Appendix heatmaps
ggsave(paste0(out_dir, "appendix_simple_heatmap_A.png"),  simple_heatmaps_A,  width = 12, height = 10, dpi = 300)
ggsave(paste0(out_dir, "appendix_simple_heatmap_B.png"),  simple_heatmaps_B,  width = 12, height = 10, dpi = 300)
ggsave(paste0(out_dir, "appendix_complex_heatmap_A.png"), complex_heatmaps_A, width = 12, height = 10, dpi = 300)
ggsave(paste0(out_dir, "appendix_complex_heatmap_B.png"), complex_heatmaps_B, width = 12, height = 10, dpi = 300)

# Main text heatmaps
ggsave(paste0(out_dir, "main_simple_heatmap_D.png"),  simple_heatmaps_D,  width = 10, height = 8,  dpi = 300)
ggsave(paste0(out_dir, "main_complex_heatmap_D.png"), complex_heatmaps_D, width = 10, height = 8,  dpi = 300)

# Main text 4 cycle
ggsave(paste0(out_dir, "simple_drug_time_4C.png"),  simple_drug_phases_4C,  width = 10, height = 8,  dpi = 300)
ggsave(paste0(out_dir, "simple_drug_time_4CN.png"),  simple_drug_phases_4CN,  width = 10, height = 8,  dpi = 300)
ggsave(paste0(out_dir, "simple_drug_time_4CS.png"),  simple_drug_phases_4CS,  width = 10, height = 8,  dpi = 300)

ggsave(paste0(out_dir, "simple_relapse_time_4C.png"),  simple_relapse_phases_4C,  width = 10, height = 8,  dpi = 300)
ggsave(paste0(out_dir, "simple_relapse_time_4CS.png"),  simple_relapse_phases_4CS,  width = 10, height = 8,  dpi = 300)
ggsave(paste0(out_dir, "simple_relapse_time_4CN.png"),  simple_relapse_phases_4CN,  width = 10, height = 8,  dpi = 300)

ggsave(paste0(out_dir, "complex_drug_time_4C.png"), complex_drug_phases_4C, width = 10, height = 8,  dpi = 300)
ggsave(paste0(out_dir, "complex_relapse_time_4C.png"), complex_relapse_phases_4C, width = 10, height = 8,  dpi = 300)

ggsave(paste0(out_dir, "simple_counts_4C.png"),          simple_counts_4C,          width = 14, height = 8,  dpi = 300)
ggsave(paste0(out_dir, "simple_counts_4CN.png"),          simple_counts_4CN,          width = 14, height = 8,  dpi = 300)
ggsave(paste0(out_dir, "simple_counts_4CS.png"),          simple_counts_4CS,          width = 14, height = 8,  dpi = 300)

ggsave(paste0(out_dir, "complex_counts_4C.png"),         complex_counts_4C,         width = 14, height = 8,  dpi = 300)

ggsave(paste0(out_dir, "simple_heatmap_P4_4C.png"),  simple_heatmaps_D_4C,  width = 10, height = 8,  dpi = 300)
ggsave(paste0(out_dir, "simple_heatmap_P4_4CN.png"),  simple_heatmaps_D_4CN,  width = 10, height = 8,  dpi = 300)
ggsave(paste0(out_dir, "simple_heatmap_P4_4CS.png"),  simple_heatmaps_D_4CS,  width = 10, height = 8,  dpi = 300)

ggsave(paste0(out_dir, "complex_heatmap_P4_4C.png"), complex_heatmaps_D_4C, width = 10, height = 8,  dpi = 300)
