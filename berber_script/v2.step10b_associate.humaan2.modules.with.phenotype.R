###
### This is step7b_associate.lipids.with.phenotype.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

###
### Step 7b - Link individual pathway (from humann2)to phenotype of interest
###
### Analogous to Step 6.
### The resulting metabolite and lipid clusters are thereafter merged into a
### combined dataset, collectively termed 'metabolite clusters', for downstream
### analyses.
###
pathway <- pathway[ctrl.no.na, ]
phenotypes2 <- phenotypes[rownames(pathway), ]
tmpMat = array (NA, c (ncol (pathway), 2,  2))
dimnames (tmpMat) [[1]] = colnames (pathway)
dimnames (tmpMat) [[2]] = c (lable1, lable2) 
dimnames (tmpMat) [[3]] = c ("estimate", "p.value")

### Associating individual molecular pathway  with HOMA-IR without de-confounding
### for BMI
tmpMat [, lable1, c ("estimate", "p.value")] =
  t (apply (pathway , MARGIN = 2, FUN = function (x) 
    unlist (cor.test (phenotypes2 [, lable1], x,
                      method = cor_method, use = "pairwise.complete.obs") [c ("estimate", "p.value")])))       

### Associating individual molecular  pathway  with HOMA-IR while de-confounding
### for BMI
tmpMat[, lable2, c ("estimate", "p.value")] =
  t (apply (pathway , MARGIN = 2, FUN = function (x) 
    unlist (pcor.test (x = x, y = phenotypes2 [, lable1], z = phenotypes2 [, c("age", "sex", "BMI")],
                       method = cor_method) [c ("estimate", "p.value")])))       

# output2 = cbind (tmpMat[, "HbA1C", c ("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, "HbA1C", "p.value"], method = "BH"), 
#                  tmpMat[, "HbA1C.adj.partial", c ("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, "HbA1C.adj.partial", "p.value"], method = "BH"))
# colnames (output2) = paste (rep (c ("HbA1C", "HbA1C.adj.partial"), each = 3), colnames (output2), sep = "_")

#output3 = output3 [with (output3, order (cluster_name, -kIN)), ]

### Write to file
cor_HbA1C [["humaan2"]] <- tmpMat
