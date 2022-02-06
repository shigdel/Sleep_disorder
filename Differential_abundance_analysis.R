

#ANCOM-BC

pseq <- makePhyloseqFromTreeSummarizedExperiment(se)

# Exclude all participants who have used antibiotics in the last four weeks
pseq = subset_samples(pseq, antibiotics != "last4weeks")

# Split by age
otu_data_a = subset_samples(pseq, adult == "1")
otu_data_a


m <- meta(otu_data_a)

ls(m)


#===============================================================================
# snoring 
#===============================================================================
# Running differential abundance analysis using ANCOM-BC, at the Genus level

# Aggregate to Genus level

Genus_data = aggregate_taxa(otu_data_a, "Genus")

tax_mat = as(tax_table(Genus_data), "matrix")
# Run ANCOM-BC
feature_table = abundances(Genus_data); meta_data = meta(Genus_data)
# ANCOM-BC requires an id column for metadata
out = ancombc(phyloseq = Genus_data, formula = "snoring + bmi+ age + gender+ smoke + educationo",
              p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
              group = "snoring", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
out

# Outputs
res = out$res
res_global = out$res_global
res_global

res_direct = out$res_direct
res_direct

res_pattern = out$res_pattern
res_pattern

tab_coef = res$beta
tab_coef
col_name = c("less than once a week", "once or twice a week", "3-5 nights/days a week", "bmi", "age", "gender", "smoke1",
             "smoke2", "educationo")
colnames(tab_coef) = col_name
tab_coef %>% datatable(caption = "Coefficients from the Primary Result") %>%
  formatRound(col_name, digits = 2)




# SEs
tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name, digits = 2)

#Test statistics
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name, digits = 2)


#P-values

tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name, digits = 2)

# Adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name, digits = 2)

#Differentially abundant taxa

tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>%
  datatable(caption = "Differentially Abundant Taxa
            from the Primary Result")

# Waterfall plot
tab_coef = res$beta
col_name = c("less than once a week", "once or twice a week", "3-5 nights/days a week",  "weight", "smoke",
             "smokeprevious", "educationo2")
colnames(tab_coef) = col_name

tab_se = res$se
colnames(tab_se) = col_name

tab_w = res$W
colnames(tab_w) = col_name

tab_p = res$p_val
colnames(tab_p) = col_name

tab_q = res$q
colnames(tab_q) = col_name

tab_diff = res$diff_abn
colnames(tab_diff) = col_name

tab_zero = out$zero_ind
tab_zero = tab_zero - tab_zero[, 1]
tab_zero = abs(tab_zero[, -1])
tab_zero = data.frame(tab_zero, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")


# Coerce the SE of structural zero to be zero
tab_se[, grepl("snoring", colnames(tab_se))] = tab_se[, grepl("snoring", colnames(tab_se))] *
  (1 - tab_zero[, grepl("snoring", colnames(tab_zero))])

# less than once a week
df1 = data.frame(genus = rownames(tab_coef),
                 lfc = tab_coef$`less than once a week`,
                 se = tab_se$`less than once a week`,
                 p = tab_p$`less than once a week`,
                 q = tab_q$`less than once a week`) %>%
  filter(q < 0.05) %>%
  arrange(desc(lfc)) %>%
  mutate(type = if_else(lfc > 0, "g1", "g2"),
         star = case_when(p < .001 ~ "***",
                          p < .01 ~ "**",
                          TRUE ~ "*"),
         pos = if_else(type == "g1", 
                       lfc + se + 0.2,
                       lfc - se - 0.2)
  )

df1$genus = factor(df1$genus, levels = df1$genus)
p1 = df1 %>%
  ggplot(aes(x = genus, y = lfc, 
             fill = type, color = type)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc - se, 
                    ymax = lfc + se), 
                width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  geom_text(aes(y = pos, label = star), 
            vjust = .7, color = "black", 
            position = position_dodge(width = 0.05)) +
  labs(x = NULL, y = "Log fold change") +
  guides(color = FALSE, drop = FALSE) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold",
                                   angle = 60, hjust = 1)) + 
  labs(title = "less than once a week") +
  scale_fill_brewer(palette = "Set1", drop = FALSE,
                    name = NULL,
                    label = c("g1" = "Increased abundance",
                              "g2" = "Reduced abundance")) +
  scale_color_brewer(palette = "Set1", drop = FALSE) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 

ggsave(plot = p1, filename = "less than once a week.jpeg", height = 5, width = 6.25, dpi = 300)



# once or twice a week


df2 = data.frame(genus = rownames(tab_coef),
                 lfc = tab_coef$`once or twice a week`,
                 se = tab_se$`once or twice a week`,
                 p = tab_p$`once or twice a week`,
                 q = tab_q$`once or twice a week`) %>%
  filter(q < 0.05) %>%
  arrange(desc(lfc)) %>%
  mutate(type = if_else(lfc > 0, "g1", "g2"),
         star = case_when(p < .001 ~ "***",
                          p < .01 ~ "**",
                          TRUE ~ "*"),
         pos = if_else(type == "g1", 
                       lfc + se + 0.2,
                       lfc - se - 0.2)
  )

df2$genus = factor(df2$genus, levels = df2$genus)
p2 = df2 %>%
  ggplot(aes(x = genus, y = lfc, 
             fill = type, color = type)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc - se, 
                    ymax = lfc + se), 
                width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  geom_text(aes(y = pos, label = star), 
            vjust = .7, color = "black", 
            position = position_dodge(width = 0.05)) +
  labs(x = NULL, y = "Log fold change") +
  guides(color = FALSE, drop = FALSE) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold",
                                   angle = 60, hjust = 1)) + 
  labs(title = "once or twice a week") +
  scale_fill_brewer(palette = "Set1", drop = FALSE,
                    name = NULL,
                    label = c("g1" = "Increased abundance",
                              "g2" = "Reduced abundance")) +
  scale_color_brewer(palette = "Set1", drop = FALSE) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 

ggsave(plot = p2, filename = "once or twice a week.jpeg", height = 5, width = 6.25, dpi = 300)


#3-5 nights/days a week

df3 = data.frame(genus = rownames(tab_coef),
                 lfc = tab_coef$`3-5 nights/days a week`,
                 se = tab_se$`3-5 nights/days a week`,
                 p = tab_p$`3-5 nights/days a week`,
                 q = tab_q$`3-5 nights/days a week`) %>%
  filter(q < 0.05) %>%
  arrange(desc(lfc)) %>%
  mutate(type = if_else(lfc > 0, "g1", "g2"),
         star = case_when(p < .001 ~ "***",
                          p < .01 ~ "**",
                          TRUE ~ "*"),
         pos = if_else(type == "g1", 
                       lfc + se + 0.2,
                       lfc - se - 0.2)
  )

df3$genus = factor(df3$genus, levels = df3$genus)
p3 = df3 %>%
  ggplot(aes(x = genus, y = lfc, 
             fill = type, color = type)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc - se, 
                    ymax = lfc + se), 
                width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  geom_text(aes(y = pos, label = star), 
            vjust = .7, color = "black", 
            position = position_dodge(width = 0.05)) +
  labs(x = NULL, y = "Log fold change") +
  guides(color = FALSE, drop = FALSE) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold",
                                   angle = 60, hjust = 1)) + 
  labs(title = "3-5 nights/days a week") +
  scale_fill_brewer(palette = "Set1", drop = FALSE,
                    name = NULL,
                    label = c("g1" = "Increased abundance",
                              "g2" = "Reduced abundance")) +
  scale_color_brewer(palette = "Set1", drop = FALSE) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 
p3

ggsave(plot = p3, filename = "3-5 nights/days a week.jpeg", height = 3, width = 2, dpi = 300) # cant save this figure 

