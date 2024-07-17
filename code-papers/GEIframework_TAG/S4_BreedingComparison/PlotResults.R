#=======================================================================
#-- Plot results
#=======================================================================

# Load necessary libraries
library(tidyverse)
library(reshape2)
library(ggpubr)

# Prepare data
temp <- read_csv("Output_PhenoModerate.csv")

temp <- 
  melt(temp, id.vars = c("year","scenario","GEI","stage"),
       measure.vars = colnames(temp)[6:dim(temp)[2]]) %>%
  group_by(year, scenario, GEI, stage, variable) %>%
  summarise(value_mean = mean(value),
            value_se   = sd(value, na.rm = TRUE) / sqrt(nReps),
            .groups    = 'drop')

temp$stage <- recode_factor(temp$stage,
                            `Parents` = "Parents",
                            `HDRW`    = "HDRW",
                            `PYT`     = "PYT",
                            `AYT`     = "AYT",
                            `EYT`     = "EYT")

plot_trends <- function(variables, title) {
  temp %>%
    filter(variable %in% variables) %>%
    droplevels() %>%
    ggplot(aes(x = year, y = value_mean)) +
    facet_grid(factor(variable) ~ factor(stage), scales = "fixed") +
    geom_line(linewidth = 0.8)
  # ylab(title) +
}

p1 <- plot_trends(c("mean_gv_tpe", "mean_gv_met"))
p2 <- plot_trends(c("var_gv_tpe", "var_gv_met"))
p3 <- plot_trends(c("acc_tpe", "acc_met", "met_tpe"))

# Plot
(p <- ggarrange(p1 + theme(axis.title = element_blank(),
                           axis.text.x.bottom = element_blank(),
                           axis.ticks.x.bottom = element_blank(),
                           plot.margin = unit(c(0, 0, 0, 0), "cm")),
                p2 + theme(strip.text.x = element_blank(),
                           axis.title = element_blank(),
                           axis.text.x.bottom = element_blank(),
                           axis.ticks.x.bottom = element_blank(),
                           plot.margin = unit(c(0, 0, 0, 0), "cm")),
                p3 + theme(axis.title = element_blank(),
                           strip.text.x = element_blank(),
                           plot.margin = unit(c(0, 0, 0, 0), "cm")), 
                nrow = 3, align = "v", heights = c(1,1,1.5), label.y = T))

# Save plot
ggsave(plot = p, filename = "PlotOutput.pdf", width = 5.5, height = 6, scale = 1.7)
