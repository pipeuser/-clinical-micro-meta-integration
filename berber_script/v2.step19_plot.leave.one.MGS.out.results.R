###
### This is step19_plot.leave.one.MGS.out.results.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

###
### Step 19 - Plotting leave-one-MGS-out results
###
### Having computed these results, we plot the top driver taxa for the
### phenotype-associated gene functions (Figure 4).
### These include sub-plots for:
### - Distribution/density plots of correlations for KOs in a KEGG module vs 
###   all other KOs
### - Distribution of correlations when leave-one-MGS-out. 
### - Distribution of correlations when leave-one-MGS-out. bg.adj. SCC 
###   (i.e. what was shown in Figure 3c+d in Pedersen et al., 2016)
### Median SCC for KO within a module (red) and all other remaining KOs (green)
### are indicated in the first two plots.
###

### Specify pdf-file for output.
pdf (file = "../output/density_plot_SCC_HOMA.IR.pdf", width = 6 * 1, height = 3 * 3)
KOsets = koann [keggmodules2test] 
for (m in names (KOsets)) {
  
  ### Distribution of correlations for KOs in modules vs all other KOs
  incat =     na.omit (KO_cor_HOMA [names (KO_cor_HOMA) %in% KOsets [[m]] ]) ### select all correlations between HOMA-IR and KOs in the KEGG module
  incat2 =    data.frame ("KO_cor_HOMA" = incat, "cat" = "KOs in module")
  notincat =  na.omit (KO_cor_HOMA [! (names (KO_cor_HOMA) %in% KOsets [[m]])]) ### select all correlations between HOMA-IR and KOs NOT in the KEGG module   
  notincat2 = data.frame ("KO_cor_HOMA" = notincat, "cat" = "KOs not in module")
  tmp.df = rbind (incat2, notincat2)
  cdf <- ddply (tmp.df, "cat", summarise, rating.median = median (KO_cor_HOMA)) ### find median for each group
  g1 <- ggplot (tmp.df, aes (x = KO_cor_HOMA, fill = cat)) + 
    geom_density (alpha = 0.4) +
    geom_vline (data = cdf, aes (xintercept = rating.median, colour = cat), linetype = "dashed", size = 0.8) +
    ggtitle (paste0 ("KEGG module: ", m, 
                     "\n", strsplit (module_mapping [m], split = "\\[|,")[[1]][1], 
                     "\n", length (KOsets [[m]]), " KOs in module vs all remaining ", (length (KO_cor_HOMA) - length (KOsets[[m]])), " KOs", "\n")) +
    xlab ("SCC for KOs and HOMA-IR") +
    theme_classic () + theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
  
  
  
  g = plot_grid (g1, ncol = 1, nrow = 1, align = "v", labels = c ("a"), axis = "rl")
  # g = plot_grid (g1, g3, ncol = 1, nrow = 2, align = "v", labels = c ("a", "b", "c"), axis = "rl") ### use this if deleting g2
  plot (g)
  
}

dev.off ()  
