library(DBI)
library(RSQLite)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
out_dir <- "C:/Users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/Plots/"
save_plot <- function(filename, plot, width = 180, height = 234, units = "mm", dpi = 300) {
  plot <- plot +
    theme(axis.text=element_text(size=12),plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 11))
  ggsave(paste0(out_dir, filename), plot = plot,
         width = width, height = height, units = units, dpi = dpi)
}

# --- Connect ---
db_path <- "C:/Users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/TS/N_Drug/Data/Time_Series.db"
con <- dbConnect(RSQLite::SQLite(), db_path)

target_run <- "Combo17_mod_13.pkl"

# --- Pull each phase for this RunID ---
burnin    <- dbGetQuery(con, "SELECT * FROM BurnIn_TimeSeries    WHERE RunID = ?", params = list(target_run))
treatment <- dbGetQuery(con, "SELECT * FROM Treatment_TimeSeries WHERE RunID = ?", params = list(target_run))
relapse   <- dbGetQuery(con, "SELECT * FROM Relapse_TimeSeries   WHERE RunID = ?", params = list(target_run))

# Sanity check — if any of these are 0 rows, the RunID string doesn't match what's stored
cat("BurnIn rows:", nrow(burnin), "| Treatment rows:", nrow(treatment), "| Relapse rows:", nrow(relapse), "\n")

if (nrow(burnin) == 0) {
  # Help find the correct RunID if the guess above was wrong
  print(dbGetQuery(con, "SELECT DISTINCT RunID FROM BurnIn_TimeSeries LIMIT 20"))
}

dbDisconnect(con)

# --- Tag each phase and stack into one long dataframe ---
burnin$Phase    <- "Natural Tumour Growth"
treatment$Phase <- "Treatment"
relapse$Phase   <- "Relapse"

all_data <- bind_rows(burnin, treatment, relapse) %>%
  mutate(Phase = factor(Phase, levels = c("Natural Tumour Growth", "Treatment", "Relapse")))

# Long format for N/I/S cell populations
pop_long <- all_data %>%
  select(Phase, time, N, I, S) %>%
  pivot_longer(cols = c(N, I, S), names_to = "CellType", values_to = "Population") %>%
  mutate(CellType = factor(CellType, levels = c("N", "I", "S"),
                           labels = c("N Cell", "I Cell", "S Cell")))

# Drug concentration only exists in the Treatment table (column C)
drug_line <- treatment %>%
  select(Phase, time, C) %>%
  rename(Drug = C)

phase_levels <- c("Natural Tumour Growth", "Treatment", "Relapse")

pop_long <- pop_long %>%
  mutate(Phase = factor(Phase, levels = phase_levels))

drug_line <- drug_line %>%
  mutate(Phase = factor(Phase, levels = phase_levels))

ann_df <- ann_df %>%
  mutate(Phase = factor(Phase, levels = phase_levels))
# --- Build annotation text from this run's own metadata ---
pattern_str <- unique(treatment$Interconversion)[1]
drug_val    <- unique(treatment$kw)[1]
ann_text <- paste0(
  "Pattern:              ", pattern_str, "\n",
  "Interconversion Rate: 0.5\n",   # adjust if you want to pull this from your summary CSV instead
  "Drug Concentration:   0.9\n",
  "Drug Target:         N cells"
)

ann_df <- data.frame(
  Phase = factor("Relapse", levels = levels(pop_long$Phase)),
  x = max(relapse$time, na.rm = TRUE) * 0.5,
  y = 0.85,
  label = ann_text
)

# --- Plot ---
NDrugTS<-ggplot() +
  geom_line(data = pop_long, aes(x = time, y = Population, color = CellType), linewidth = 1) +
  geom_line(data = drug_line, aes(x = time, y = 0.5), color = "red",
            linetype = "dashed", linewidth = 1.2) +
  geom_label(data = ann_df, aes(x = x, y = y, label = label),
             size = 3, hjust = 0.5, fill = "lightyellow", label.size = 0.3) +
  facet_wrap(~Phase, scales = "free_x") +
  # Shared colour palette
  scale_color_manual(values = c("N Cell" =  "#7CAE00", "I Cell" = "#F8766D", "S Cell" = "#C77CFF")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Time", y = "Population", color = NULL) +
  theme_classic() +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

save_plot("NDrugTS.pdf",              NDrugTS)




# --- Connect ---
db_path <- "C:/Users/Student/OneDrive - hull.ac.uk/Desktop/Year3/Experiment_3/TS/S_Drug/Data/Time_Series.db"
con <- dbConnect(RSQLite::SQLite(), db_path)

target_run <- "Combo17_mod_52.pkl"

# --- Pull each phase for this RunID ---
burnin    <- dbGetQuery(con, "SELECT * FROM BurnIn_TimeSeries    WHERE RunID = ?", params = list(target_run))
treatment <- dbGetQuery(con, "SELECT * FROM Treatment_TimeSeries WHERE RunID = ?", params = list(target_run))
relapse   <- dbGetQuery(con, "SELECT * FROM Relapse_TimeSeries   WHERE RunID = ?", params = list(target_run))

# Sanity check — if any of these are 0 rows, the RunID string doesn't match what's stored
cat("BurnIn rows:", nrow(burnin), "| Treatment rows:", nrow(treatment), "| Relapse rows:", nrow(relapse), "\n")

if (nrow(burnin) == 0) {
  # Help find the correct RunID if the guess above was wrong
  print(dbGetQuery(con, "SELECT DISTINCT RunID FROM BurnIn_TimeSeries LIMIT 20"))
}

dbDisconnect(con)

# --- Tag each phase and stack into one long dataframe ---
burnin$Phase    <- "Natural Tumour Growth"
treatment$Phase <- "Treatment"
relapse$Phase   <- "Relapse"

all_data <- bind_rows(burnin, treatment, relapse) %>%
  mutate(Phase = factor(Phase, levels = c("Natural Tumour Growth", "Treatment", "Relapse")))

# Long format for N/I/S cell populations
pop_long <- all_data %>%
  select(Phase, time, N, I, S) %>%
  pivot_longer(cols = c(N, I, S), names_to = "CellType", values_to = "Population") %>%
  mutate(CellType = factor(CellType, levels = c("N", "I", "S"),
                           labels = c("N Cell", "I Cell", "S Cell")))

# Drug concentration only exists in the Treatment table (column C)
drug_line <- treatment %>%
  select(Phase, time, C) %>%
  rename(Drug = C)

phase_levels <- c("Natural Tumour Growth", "Treatment", "Relapse")

pop_long <- pop_long %>%
  mutate(Phase = factor(Phase, levels = phase_levels))

drug_line <- drug_line %>%
  mutate(Phase = factor(Phase, levels = phase_levels))

ann_df <- ann_df %>%
  mutate(Phase = factor(Phase, levels = phase_levels))
# --- Build annotation text from this run's own metadata ---
pattern_str <- unique(treatment$Interconversion)[1]
drug_val    <- unique(treatment$kw)[1]
ann_text <- paste0(
  "Pattern:              ", pattern_str, "\n",
  "Interconversion Rate: 0.5\n",   # adjust if you want to pull this from your summary CSV instead
  "Drug Concentration:   0.9\n",
  "Drug Target:         S cells"
)

ann_df <- data.frame(
  Phase = factor("Relapse", levels = levels(pop_long$Phase)),
  x = max(relapse$time, na.rm = TRUE) * 0.5,
  y = 0.85,
  label = ann_text
)

# --- Plot ---
SDrugTS<-ggplot() +
  geom_line(data = pop_long, aes(x = time, y = Population, color = CellType), linewidth = 1) +
  geom_line(data = drug_line, aes(x = time, y = 0.5), color = "red",
            linetype = "dashed", linewidth = 1.2) +
  geom_label(data = ann_df, aes(x = x, y = y, label = label),
             size = 3, hjust = 0.5, fill = "lightyellow", label.size = 0.3) +
  facet_wrap(~Phase, scales = "free_x") +
  # Shared colour palette
  scale_color_manual(values = c("N Cell" =  "#7CAE00", "I Cell" = "#F8766D", "S Cell" = "#C77CFF")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Time", y = "Population", color = NULL) +
  theme_classic() +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

save_plot("SDrugTS.pdf",              SDrugTS)
