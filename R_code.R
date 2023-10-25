library(lme4)
library(ggpubr)
library(lmerTest)
library(magrittr)
library(parallel)
library(reshape2)
library(tidyverse)
library(missForest)
library(RColorBrewer)



source("Source R code.txt")


## ========================= Load data ==================================== ##
Europe_Forest_raw = read.csv("Europe_Forest_raw.csv")

Europe_Forest_species = read.csv("Europe_Forest_species.csv", sep = '"')[,c('block', 'plot', 'full_species_original', 'X.6')]
colnames(Europe_Forest_species)[4] = "basal_area"
Europe_Forest_species$block = unlist(strsplit(Europe_Forest_species$block, split = " "))
Europe_Forest_species$basal_area = as.numeric(Europe_Forest_species$basal_area)

tmp = Europe_Forest_species[Europe_Forest_species$block == "Germany",]
for (i in 1:3) {
  tmp[,i] = factor(tmp[,i])
}

Europe_Forest_species[Europe_Forest_species$block == "Germany", "basal_area"] = missForest(tmp[,-1])$ximp$basal_area
rm(tmp)



for (x in c("FIN02", "FIN03", "FIN04", "FIN07", "FIN08", "FIN12", 
            "FIN13", "FIN15", "FIN20", "FIN24", "FIN25", "FIN26", "FIN28")) {
  
  Europe_Forest_species[Europe_Forest_species$plot == x & Europe_Forest_species$full_species_original == "Betula.pendula", "basal_area"] =
    sum(Europe_Forest_species[Europe_Forest_species$plot == x & Europe_Forest_species$full_species_original %in% c("Betula.pendula", "Betula.pubescens"), "basal_area"])
  
  Europe_Forest_species = Europe_Forest_species[ -which(Europe_Forest_species$plot == x & Europe_Forest_species$full_species_original == "Betula.pubescens"),]
  
}

Europe_Forest_species[Europe_Forest_species$plot == "GER05" & Europe_Forest_species$full_species_original == "Quercus.petraea", "basal_area"] =
  sum(Europe_Forest_species[Europe_Forest_species$plot == "GER05" & Europe_Forest_species$full_species_original %in% c("Quercus.petraea", "Quercus.robur"), "basal_area"])

Europe_Forest_species = Europe_Forest_species[ -which(Europe_Forest_species$plot == "GER05" & Europe_Forest_species$full_species_original == "Quercus.robur"),]

Europe_Forest_species$full_species_original[Europe_Forest_species$full_species_original == "Betula.pubescens"] = "Betula.pendula"
Europe_Forest_species$full_species_original[Europe_Forest_species$full_species_original == "Quercus.robur"] = "Quercus.petraea"


##
variables = colnames(Europe_Forest_raw)[-(1:5)]

Europe_Forest_raw = standardize.Europe_Forest_raw(Europe_Forest_raw)
variables.std <- paste0(variables, ".std")


correlation = Europe_Forest_raw[,variables.std]
correlation = cor(correlation)
distM = sqrt(1 - abs(correlation))



## ============================= Plot Figure 2 ======================================== ##
Europe_Forest_raw <- Europe_Forest_raw %>%
  mutate(mf_Chao_0 = apply(Europe_Forest_raw[,variables.std], 1, function(x) MF.uncor(x, rep(1, length(x)), 0)$qMF),
         mf_Chao_1 = apply(Europe_Forest_raw[,variables.std], 1, function(x) MF.uncor(x, rep(1, length(x)), 1)$qMF),
         mf_Chao_2 = apply(Europe_Forest_raw[,variables.std], 1, function(x) MF.uncor(x, rep(1, length(x)), 2)$qMF),
         
         mf_Chao_AUC_0 = apply(Europe_Forest_raw[,variables.std], 1, function(x) MF.cor(x, rep(1, length(x)), distM, q = 0) %>% filter(tau == 'AUC') %>% select(qMF) %>% as.numeric),
         mf_Chao_AUC_1 = apply(Europe_Forest_raw[,variables.std], 1, function(x) MF.cor(x, rep(1, length(x)), distM, q = 1) %>% filter(tau == 'AUC') %>% select(qMF) %>% as.numeric),
         mf_Chao_AUC_2 = apply(Europe_Forest_raw[,variables.std], 1, function(x) MF.cor(x, rep(1, length(x)), distM, q = 2) %>% filter(tau == 'AUC') %>% select(qMF) %>% as.numeric)
         )

mf = Europe_Forest_raw %>%
  select(plotid, target_species_richness, mf_Chao_0:mf_Chao_AUC_2) %>%
  pivot_longer(cols = c(mf_Chao_0:mf_Chao_AUC_2), names_to = "method") %>%
  mutate(method = fct_inorder(method))

## Compute species diversity
mf = mf %>% mutate(target_species_richness = sapply(unique(Europe_Forest_species$plot), function(x) 
  rep( qD( (Europe_Forest_species %>% filter(plot == x))$basal_area, q = c(0,1,2)), 2)) %>% as.vector)

mf = mf %>% mutate(Order.q = rep(c(0, 1, 2), 2*nrow(Europe_Forest_raw)),
                   Div = rep(rep(c('Uncorrected', 'Corrected'), each = 3), nrow(Europe_Forest_raw)))



## Fit linear mixed model
mf = mf %>% arrange(method) 
mf = mf %>% group_by(method) %>% 
  do(lmer(formula = value ~ 1 + target_species_richness + (1 + target_species_richness | plotid), data = . ) %>% 
       predict %>% tibble(fit = .)) %>% 
  ungroup %>% select(fit) %>% bind_cols(mf) %>% mutate("p_value" = 0, 'Sig' = "Significant slope (P < 0.05)")


mf$Order.q = as.factor(mf$Order.q)


slope = lapply(1:length(unique(mf$method)), function(j) {
  
  myout_ <- mf %>% filter(method == unique(mf$method)[j])
  
  model = lmer(formula = value ~ 1 + target_species_richness + (1 + target_species_richness | plotid), 
               data    = myout_)
  
  ##
  slope = coef(model)$plotid$target_species_richness
  slope = ifelse(round(abs(slope), 3) < 0.001, slope %>% round(., 4) %>% paste("Slope = ", ., sep = ""), 
                 ifelse(round(abs(slope), 2) < 0.01, slope %>% round(., 3) %>% paste("Slope = ", ., sep = ""),
                        format(slope %>% round(., 2), nsmall = 2) %>% paste("Slope = ", ., sep = "")))
  
  tmp = myout_[ !duplicated(myout_[, c('plotid', 'Order.q', 'Div')]),] %>% select(c('plotid', 'Order.q', 'Div'))
  
  ##
  slope_lmm = summary(model)$coefficients[2, 'Estimate']
  slope_lmm = ifelse(round(abs(slope_lmm), 3) < 0.001, slope_lmm %>% round(., 4) %>% paste("slope = ", ., sep = ""), 
                 ifelse(round(abs(slope_lmm), 2) < 0.01, slope_lmm %>% round(., 3) %>% paste("slope = ", ., sep = ""),
                        format(slope_lmm %>% round(., 2), nsmall = 2) %>% paste("slope = ", ., sep = "")))
  
  tmp_lmm = myout_[ !duplicated(myout_[, c('Order.q', 'Div')]),] %>% select(c('Order.q', 'Div'))
  
  ##
  rbind(cbind(tmp, x = max(myout_$target_species_richness), y = max(myout_$value), Slope = slope),
        cbind(plotid = 'Linear mixed', tmp_lmm, x = max(myout_$target_species_richness), y = max(myout_$value), Slope = slope_lmm))
  
  }) %>% do.call(rbind,.)

slope = slope %>% mutate(plotid = fct_inorder(plotid)) %>% arrange(plotid)


mf = rbind(mf,
           
           lapply(unique(mf$method), function(x) {
             lmm.data = mf %>% filter(method == x)
             
             model <- lmer(formula = value ~ 1 + target_species_richness + (1 + target_species_richness | plotid), 
                           data    = lmm.data)
             
             data.frame(fit = predict(model, re.form = NA), 
                        plotid = 'Linear mixed', 
                        target_species_richness = lmm.data$target_species_richness, 
                        value = 0, 
                        method = x, 
                        Order.q = unique(lmm.data$Order.q), 
                        Div = unique(lmm.data$Div),
                        p_value = summary(model)$coefficients[2, 'Pr(>|t|)'],
                        Sig = ifelse( summary(model)$coefficients[2, 'Pr(>|t|)'] < 0.05, 
                                      'Significant slope (P < 0.05)', 'Nonsignificant slope' ))
             }) %>% do.call(rbind,.)
           )


fig.theme = theme(legend.position = "bottom", 
                  legend.box = "vertical", 
                  legend.key.width = unit(1.2, "cm"), 
                  text = element_text(size = 16), 
                  plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
                  strip.text = element_text(size = 12, face = "bold"))

guide = guides(colour = guide_legend(title = "Plot ID", override.aes = list(linewidth = 1.5)),
               lty = guide_legend(override.aes = list(linewidth = 1, colour = "red")),
               size = "none")

col_manual = scale_colour_manual(breaks = c("FIN", "GER", "ITA", "POL", 
                                            "ROM", "SPA", "Linear mixed"),
                                 label = c("Finland (3)", "Germany (5)", "Italy (5)",  
                                           "Poland (5)", "Romania (4)", "Spain (4)", "Linear mixed"),
                                 values = c("FIN" = "black", "GER" = "purple2", "ITA" = "darkorange", "POL" = "steelblue1", 
                                            "ROM" = "blue", "SPA" = "gray55", "Linear mixed" = "red"))


size_manual = scale_size_manual(values = c("FIN" = 0.5, "GER" = 0.5, "ITA" = 0.5, "POL" = 0.5, 
                                           "ROM" = 0.5, "SPA" = 0.5, "Linear mixed" = 1.9))

lty_manual = scale_linetype_manual(values = c("Nonsignificant slope" = "dashed", "Significant slope (P < 0.05)" = "solid"), name = NULL, drop = FALSE)


## Figure 2
Fig2 = ggplot() +
  geom_line(data = mf %>% mutate(Div = fct_inorder(Div)),
            aes(x = target_species_richness, y = fit, col = plotid, lty = Sig, size = plotid)) +
  geom_point(data = mf %>% filter(plotid != 'Linear mixed') %>% 
               mutate(Div = fct_inorder(Div)),
             aes(x = target_species_richness, y = value, color = plotid), alpha = 0.2) +
  geom_text(data = slope %>% mutate(Div = fct_inorder(Div),
                                    x = max(x) - rep(c(3, 1, 3, 1, 3, 1, 3), each = 6), 
                                    y = max(y) + c(rep(0.7, 12), rep(0.2, 12), rep(-0.3, 12), rep(-0.8, 6))),
            aes(x = x, y = y, label = Slope, color = plotid), size = 5, show.legend = FALSE) +
  theme_bw() +
  facet_grid(Div ~ Order.q, scale = "fixed", 
             labeller = labeller(Order.q = c(`0` = "q = 0", `1` = "q = 1", `2` = "q = 2"),
                                 Div = c(`Uncorrected` = "Uncorrected for correlations", 
                                         `Corrected` = "Corrected for correlations"))) +
  labs(x = "Species diversity", y = "Multifunctionality") +
  col_manual +
  size_manual +
  lty_manual +
  fig.theme +
  guide +
  coord_cartesian(ylim = c(7.3, 15.2))

ggsave("Figure 2.pdf", plot = Fig2, width = 10, height = 10, dpi = 1000)



## =========================== Plot Figure 3, 4 ========================================= ##
cpu.cores <- detectCores()-1
cl <- makeCluster(cpu.cores)
clusterExport(cl, varlist = c("Europe_Forest_raw", "variables.std", "distM", "qD", "MF.uncor", "MF.cor", "beta.MF.uncor", "beta.MF.cor", "Europe_Forest_species"), 
              envir = environment())
clusterEvalQ(cl, c(library(dplyr), library(tidyr), library(reshape2)))


beta.result.uncor = parLapply(cl, unique(Europe_Forest_raw$plotid)[-1], function(x) {
  
  country.data = Europe_Forest_raw %>% filter(plotid == x)
  comb = combn(1:nrow(country.data), 2)
  
  species.data = Europe_Forest_species[substr(Europe_Forest_species$plot, 1, 3) == x,]
  
  lapply(1:ncol(comb), function(i) {
    
    index = comb[,i]
    
    data = country.data[index, ]
    N = length(index)
    
    Species = species.data %>% filter(plot %in% unique(species.data$plot)[index]) %>% acast(., full_species_original ~ plot, value.var = "basal_area")
    Species[is.na(Species)] = 0
    
    Gamma.div = qD(rowSums(Species), q = c(0,1,2))
    Alpha.div = qD(as.vector(Species), q = c(0,1,2)) / N
    Beta.div = Gamma.div / Alpha.div
    
    out = beta.MF.uncor(data[,variables.std] %>% t, q = c(0,1,2))
    
    out %>% filter(Type %in% c('Gamma','Alpha','Beta')) %>% 
      mutate(Gamma = rep(Gamma.div, 3), Alpha = rep(Alpha.div, 3), Beta = rep(Beta.div, 3)) %>%
      pivot_longer(cols = c(Gamma:Beta), names_to = 'Dissimilarity')
    
  }) %>% do.call(rbind,.) %>% mutate(plotid = x)
  
}) %>% do.call(rbind,.)


beta.result.cor = parLapply(cl, unique(Europe_Forest_raw$plotid)[-1], function(x) {
  
  country.data = Europe_Forest_raw %>% filter(plotid == x)
  comb = combn(1:nrow(country.data), 2)
  
  species.data = Europe_Forest_species[substr(Europe_Forest_species$plot, 1, 3) == x,]
  
  lapply(1:ncol(comb), function(i) {
    
    index = comb[,i]
    
    data = country.data[index, ]
    N = length(index)
    
    Species = species.data %>% filter(plot %in% unique(species.data$plot)[index]) %>% acast(., full_species_original ~ plot, value.var = "basal_area")
    Species[is.na(Species)] = 0
    
    Gamma.div = qD(rowSums(Species), q = c(0,1,2))
    Alpha.div = qD(as.vector(Species), q = c(0,1,2)) / N
    Beta.div = Gamma.div / Alpha.div
    
    out = beta.MF.cor(data[,variables.std] %>% t, distM = distM, q = c(0,1,2)) %>% filter(tau == 'AUC') %>% select(-tau)
    
    out %>% filter(Type %in% c('Gamma','Alpha','Beta')) %>% 
      mutate(Gamma = rep(Gamma.div, 3), Alpha = rep(Alpha.div, 3), Beta = rep(Beta.div, 3)) %>%
      pivot_longer(cols = c(Gamma:Beta), names_to = 'Dissimilarity')
    
  }) %>% do.call(rbind,.) %>% mutate(plotid = x)
  
}) %>% do.call(rbind,.)


stopCluster(cl)




## Figure 3
fig3 = fig_alpha_gamma_beta(beta.result.uncor, type = "uncorrected")

fig3$fig.alpha  ## figure 3a
fig3$fig.beta   ## figure 3b
fig3$fig.gamma  ## figure 3c

ggsave("Figure 3a.pdf", plot = fig3$fig.alpha, width = 9,  height = 6.5, dpi = 1000)
ggsave("Figure 3b.pdf", plot = fig3$fig.beta,  width = 10, height = 6.5, dpi = 1000)
ggsave("Figure 3c.pdf", plot = fig3$fig.gamma, width = 9,  height = 7,   dpi = 1000)


library(ggpubr)
ggsave("Figure 3.pdf", 
       plot = ggarrange(fig3$fig.alpha + ggtitle("(a)") + theme(plot.title = element_text(size = 20, face = "bold")),
                        fig3$fig.beta  + ggtitle("(b)") + theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.06)),
                        fig3$fig.gamma + ggtitle("(c)") + theme(plot.title = element_text(size = 20, face = "bold")),
                        nrow = 3, heights = c(6, 6, 8)), 
       width = 10,  height = 21, dpi = 1000)



## figure 4
fig4 = fig_alpha_gamma_beta(beta.result.cor,   type = "corrected")

fig4$fig.alpha  ## figure 4a
fig4$fig.beta   ## figure 4b
fig4$fig.gamma  ## figure 4c


ggsave("Figure 4a.pdf", plot = fig4$fig.alpha, width = 9,  height = 6.5, dpi = 1000)
ggsave("Figure 4b.pdf", plot = fig4$fig.beta,  width = 10, height = 6.5, dpi = 1000)
ggsave("Figure 4c.pdf", plot = fig4$fig.gamma, width = 9,  height = 7,   dpi = 1000)


ggsave("Figure 4.pdf", 
       plot = ggarrange(fig4$fig.alpha + ggtitle("(a)") + theme(plot.title = element_text(size = 20, face = "bold")),
                        fig4$fig.beta  + ggtitle("(b)") + theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.06)),
                        fig4$fig.gamma + ggtitle("(c)") + theme(plot.title = element_text(size = 20, face = "bold")),
                        nrow = 3, heights = c(6, 6, 8)), 
       width = 10,  height = 21, dpi = 1000)



