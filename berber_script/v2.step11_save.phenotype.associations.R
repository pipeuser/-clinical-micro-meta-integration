###
### This is step11_save.phenotype.associations.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

###
### Step 11 - Save phenotype associations
###
### Save the (BMI corrected) HbA1C association of metabolite clusters, MGSs
### and KEGG modules calculated in Step 8-10.
###

tmp.excel.file = paste0(outdir, "/",lable1,"_associations.xlsx")
write.xlsx (paste0 ("Associations with ,", lable1, "and" ,lable1 ,"adjusted for BMI"), 
            sheetName = "info", file = tmp.excel.file, row.names = F, col.names = F)

for (tmp_name in names (cor_HbA1C)) {
  tmpMat = cor_HbA1C [[paste (tmp_name)]]
  out = cbind (tmpMat [, lable1, c("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, lable1,"p.value"], method = "BH"), 
               tmpMat [, lable2, c ("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [, lable2, "p.value"], method = "BH"))
  colnames (out) = paste (rep (c (lable1,  lable2), each = 3), colnames (out), sep = "_")
 # if (tmp_name == "metlip") { rownames (out) = cluster_mapping_file [match( rownames (out), cluster_mapping_file$New_Name), "label"] }
  write.xlsx (x = out, file = tmp.excel.file, append = T, sheetName = paste (tmp_name, sep = "_"), row.names = T, col.names = T)            
}
