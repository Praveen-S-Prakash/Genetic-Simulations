
workingDir <- "/Users/praveenp/Desktop/Slim_simulations_100_replicate/"
setwd(workingDir)

library(ggplot2)


# Custom colours for selfing levels

custom_colors <- c(
  "0" = "#004D40",
  "25" = "#A5D6A7",
  "50" = "#FFCC80",
  "75" = "#F57F17",
  "100" = "#BF360C"
)


###### FST: invasive vs invasive ######

FST_p2_p3 <- read.csv(
  "fst_invaive_invasive_slim_final_ggpl.csv",
  header = TRUE
)

FST_p2_p3$Selfing <- as.factor(FST_p2_p3$Selfing)

FST_p2_p3$Scenario <- factor(
  FST_p2_p3$Scenario,
  levels = unique(FST_p2_p3$Scenario)
)

ggplot(
  FST_p2_p3,
  aes(x = Scenario, y = Value, fill = Selfing)
) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Differentiation between two invasive populations",
    x = "Scenario",
    y = "FST"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.5
    )
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = custom_colors)


###### FST: native vs invasive ######

FST_p1_p2 <- read.csv(
  "fst_native_invasive_final_ggpl.csv",
  header = TRUE
)

FST_p1_p2$Selfing <- as.factor(FST_p1_p2$Selfing)

FST_p1_p2$Scenario <- factor(
  FST_p1_p2$Scenario,
  levels = unique(FST_p1_p2$Scenario)
)

ggplot(
  FST_p1_p2,
  aes(x = Scenario, y = Value, fill = Selfing)
) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Differentiation between native and invasive population",
    x = "Scenario",
    y = "FST"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.5
    )
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = custom_colors)


###### Heterozygosity (He) ######

het_p1_p2_box <- read.csv(
  "He_smil_final_ggpl.csv",
  header = TRUE
)

het_p1_p2_box$Selfing <- as.factor(het_p1_p2_box$Selfing)

het_p1_p2_box$Scenario <- factor(
  het_p1_p2_box$Scenario,
  levels = unique(het_p1_p2_box$Scenario)
)

ggplot(
  het_p1_p2_box,
  aes(x = Scenario, y = Value, fill = Selfing)
) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Heterozygocity",
    x = "Scenario",
    y = "He"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.5
    )
  ) +
  scale_fill_manual(values = custom_colors)


###### Nucleotide diversity (Theta) ######

theta_p1_p2_box <- read.csv(
  "theta_slim_final_ggpl.csv",
  header = TRUE
)

theta_p1_p2_box$Selfing <- as.factor(theta_p1_p2_box$Selfing)

theta_p1_p2_box$Scenario <- factor(
  theta_p1_p2_box$Scenario,
  levels = unique(theta_p1_p2_box$Scenario)
)

ggplot(
  theta_p1_p2_box,
  aes(x = Scenario, y = Value, fill = Selfing)
) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Nucleotide Diversity",
    x = "Scenario",
    y = "Theta"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.5
    )
  )
