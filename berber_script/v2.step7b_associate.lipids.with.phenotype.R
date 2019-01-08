###
### This is step7b_associate.lipids.with.phenotype.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

###
### Step 7b - Link individual molecular lipids to phenotype of interest
###
### Analogous to Step 6.
### The resulting metabolite and lipid clusters are thereafter merged into a
### combined dataset, collectively termed 'metabolite clusters', for downstream
### analyses.
###
phenotypes2 <- phenotypes[rownames(lipidomic), ]
#phenotypes2 <- phenotypes2[, lable1]
tmpMat = array (NA, c (ncol (lipidomic), 2,  2))
dimnames (tmpMat) [[1]] = colnames (lipidomic)
dimnames (tmpMat) [[2]] = c (lable1,  lable2) 
dimnames (tmpMat) [[3]] = c ("estimate", "p.value")

### Associating individual molecular lipids with HOMA-IR without de-confounding
### for BMI
tmpMat [, lable1, c ("estimate", "p.value")] =
  t (apply (lipidomic , MARGIN = 2, FUN = function (x) 
    unlist (cor.test (phenotypes2 [, lable1], x,
                      method = cor_method, use = "pairwise.complete.obs") [c ("estimate", "p.value")])))       

### Associating individual molecular lipids with HOMA-IR while de-confounding
### for BMI
tmpMat[, lable2, c ("estimate", "p.value")] =
  t (apply (lipidomic , MARGIN = 2, FUN = function (x)
    unlist (pcor.test (x = x, y = phenotypes2 [, lable1], z = phenotypes2 [, c("age", "sex", "BMI")],
                       method = cor_method) [c ("estimate", "p.value")])))       

### Sort by cluster_name and then decreasing values of kIN
output2 = cbind (tmpMat[, lable1, c ("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, lable1, "p.value"], method = "BH"), 
                 tmpMat[, lable2, c ("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, lable2, "p.value"], method = "BH"))
colnames (output2) = paste (rep (c (lable1, lable2), each = 3), colnames (output2), sep = "_")
output3 = cbind (output, output2)
#output3 = output3 [with (output3, order (cluster_name, -kIN)), ]

### Write to file
write.csv (output3, file = paste0(outdir, "/individual_bileacid.csv"), quote = F, row.names = T)
rm (output, output2, output3, tmpMat, dat, dat_tmp)

### Create a joint data frame of metabolite/lipid cluster eigengene
### equivalents/effective abundances ('MEsMetLip')
namesinter <- intersect(rownames(MEsMeta), rownames(MEsLipid))
MEsMetLip = cbind (MEsMeta [namesinter,], MEsLipid [namesinter,])
### Rename module names (columnames) from 'colors' to numbers
#colnames (MEsMetLip) = cluster_mapping_file [c (paste0 ("M_", colnames (MEsMeta)), paste0 ("L_", colnames (MEsLipid))), "New_Name"]
### Exclude the two bin-clusters (i.e. "M_remaining" and "L_remaining")
#MEsMetLip = subset (MEsMetLip, select = c (-M_remaining, -L_remaining))
### Save MEsMetLip to file, order by module number
write.table (MEsMetLip [ , order (colnames (MEsMetLip))], file = paste0(outdir, "/MEs_metabolite_clusters.txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
