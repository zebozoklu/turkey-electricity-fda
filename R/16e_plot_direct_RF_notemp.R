library(dplyr)
library(ggplot2)
library(readr)

by_hour <- read_csv("output/tables/direct_RF_by_hour.csv",
                    show_col_types = FALSE)

keep <- c(
  "Direct RF (no temperature)",
  "FPCA VAR(1,7) + holiday/load/derivative regime",
  "Official forecast"
)

labels <- c(
  "Direct RF (no temperature)" = "Direct RF",
  "FPCA VAR(1,7) + holiday/load/derivative regime" = "FPCA-OLS",
  "Official forecast" = "Official forecast"
)

p <- by_hour |>
  filter(model %in% keep) |>
  mutate(model = labels[model],
         model = factor(model, levels = c("Direct RF", "FPCA-OLS", "Official forecast"))) |>
  ggplot(aes(hour, mae, colour = model, linetype = model)) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(breaks = seq(0, 23, 3)) +
  scale_colour_manual(values = c(
    "Direct RF"         = "#E07B54",
    "FPCA-OLS"          = "#3A9BAD",
    "Official forecast" = "#A86BB5"
  )) +
  scale_linetype_manual(values = c(
    "Direct RF"         = "solid",
    "FPCA-OLS"          = "dashed",
    "Official forecast" = "dashed"
  )) +
  labs(title = "MAE by hour: Direct RF vs FPCA-OLS (2023)",
       x = "Hour", y = "MAE (MWh)", colour = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("output/figs/direct_RF_mae_by_hour_notemp.png", p,
       width = 9, height = 5, dpi = 150)
cat("Saved: output/figs/direct_RF_mae_by_hour_notemp.png\n")
