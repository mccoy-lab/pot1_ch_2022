library(tidyverse)
library(data.table)
library(pbmcapply)
library(RColorBrewer)

### Analysis code for simulation data.

#=========================#
# Read in simulation data #
#=========================#

read_data <- function(bg, i, telo_len) {
  # find path to file (each telo background is in a separate directory)
  path <- paste0(bg, "_", telo_len, "/",
                 bg, "_", telo_len, "_", i, ".out")
  
  # reformat slim output into data table
  if (file.exists(path)) {
    dt <- fread(path) %>%
      pivot_wider(names_from = "gen",
                  values_from = c("p1", "p2")) %>%
      setDT() %>%
      # calculate VAF
      .[, vaf := p2_size / (p1_size + p2_size)] %>%
      .[, telo := telo_len] %>%
      .[, run := i] %>%
      .[, id := paste0(bg, "_", telo)] %>%
      setnames(c("V1"), c("gen")) %>%
      .[, c("gen", "p1_size", "p2_size", "p1_telo", "vaf", "p2_telo", "telo", "run", "id")] %>%
      # optional - fix divide by zero errors for VAF
      .[is.na(vaf), vaf := 0]
    return(dt)
  }
}
# read in all data
sims <- rbind(pbmclapply(1:10000,
                         function(x) read_data("pot1", x, 13400)) %>%
                rbindlist(),
              pbmclapply(as.list(c(8600, 11000, 13400)),
                         function(x) pbmclapply(1:10000,
                                                function(y) read_data("wt", y, x)) %>%
                           rbindlist()) %>%
                rbindlist())
sims$id <- factor(sims$id, levels = c("wt_8600", "wt_11000", "wt_13400", "pot1_13400"))
fwrite(sims, "sims.txt", sep = "\t")

# set x ticks for plotting
xticks <- seq(0, 90*365, 365*10)
xticknames <- seq(0, 90, 10)


#========================#
# JAK2 VAF vs. age plots #
#========================#

# calculate mean VAF and CI/SD at each time point
get_mean <- function(ident) {
  mean <- sims[id == ident] %>%
    group_by(gen, id) %>%
    summarize(mean = mean(vaf),
              sd = sd(vaf),
              q1 = quantile(vaf, 0.05),
              q2 = quantile(vaf, 0.95)) %>%
    setDT()
  return(mean)
}
mean <- pbmclapply(as.list(c("wt_8600", "wt_11000", "wt_13400", "pot1_13400")),
                   function(x) get_mean(x),
                   mc.cores = 12L) %>%
  rbindlist()
fwrite(mean, "sim_mean.txt", sep = "\t")
mean$id <- factor(mean$id, levels = c("wt_8600", "wt_11000", "wt_13400", "pot1_13400"))

# plot mean and confidence intervals over time
p <- ggplot(data = mean) +
  geom_line(aes(x = gen, y = mean, color = id),
            size = 1, lineend = "round") +
  # show 95% CI in geom_ribbon
  geom_ribbon(aes(x = gen, ymin = q1, ymax = q2, fill = id),
              alpha = 0.08, linetype = "blank", show.legend = FALSE) +
  scale_x_continuous(breaks = xticks, labels = xticknames) +
  scale_color_brewer(palette = "Set1", aesthetics = c("color", "fill"),
                     labels = c("8.6 kb", "11 kb", "13.4 kb", "13.4 kb (POT1)")) +
  labs(x = "Years of age", y = "CH driver VAF", color = "Telomere length") +
  theme_classic() +
  theme(legend.position = c(.21, .70))
ggsave(paste0("mean_ci.pdf"), p,
       width = 2.2*2.2, height = 2.2*1.47)


#================================================#
# Percent of simulations with >2% VAF JAK2 clone #
#================================================#

# summarize percent of simulations with VAF > 2% each day
# doing it separately for each telomere length to avoid reading in too much data
percvaf <- sims[id == "pot1_13400"] %>%
  group_by(gen, id) %>%
  summarize(perc_over_thresh = sum(vaf > 0.02)/100) %>%
  as.data.table()
for (i in c(8600, 11000, 13400)) {
  percvaf <- rbind(percvaf,
                   sims[id == paste0("wt_", i)] %>%
                     group_by(gen, id) %>%
                     summarize(perc_over_thresh = sum(vaf > 0.02)/100) %>%
                     as.data.table()) 
}
# fix divide by zero errors
percvaf[is.na(perc_over_thresh), perc_over_thresh := 0]
fwrite(percvaf, "2perc_vaf.txt", sep = "\t")
percvaf$id <- factor(percvaf$id,
                     levels = c("wt_8600", "wt_11000", "wt_13400", "pot1_13400"))

# plot percentage over time
p <- ggplot(percvaf,
            aes(x = gen, y = perc_over_thresh / 100, color = id)) +
  geom_line() +
  theme_classic() +
  scale_x_continuous(breaks = ticks, labels = ticknames) +
  scale_color_brewer(palette = "Set1", aesthetics = c("colour", "fill"),
                     labels = c("8.6 kb", "11 kb", "13.4 kb", "13.4 kb (POT1)")) +
  labs(x = "Years of age", y = "Proportion of >2% CH VAF simulations", color = "Telomere length") +
  theme(legend.position = c(.21, .70))
ggsave("2perc_vaf.pdf", p, width = 2.2*2.2, height = 2.2*1.47)