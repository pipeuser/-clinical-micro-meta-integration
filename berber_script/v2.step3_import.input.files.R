###
### This is step3_import.input.files.R, part of the Nature Protocol: 
### 'A computational framework for conducting a three-pronged association study
### to integrate high-throughput ‘-omics’ datasets for the identification of
### potential mechanistic links'.
### 26/7-2018
###

###
### Step 3 - Import input files
### 
### The following input data is assumed to exist within the top/data subdirectory:
###
### Data-files: 
### -	phenotypes.tab: 
###   File with clinical phenotypes (columns) per
###   individual (rows). Used to test for associations with, or for confounder
###   analysis. Individuals are labeled ‘idv’ followed by a number, e.g. ‘idv001’.
###
### -	metabolomic.tab / lipidomic.tab: 
###   Input data matrix for abundance of 325
###   polar metabolites or 876 molecular lipids per individual. Note that no
###   additional normalization is done in this script, so data is assumed to be
###   comparable in these regards. Such different data types are eventually merged
###   into a single set of metabolite cluster abundances. Individual
###   metabolite/lipids are named M or L (for specifying a polar metabolite or
###   molecular lipid, respectively), followed by a number and lastly the
###   annotation or ‘unknown’ in case of unannotated metabolites/lipids, e.g.
###   ‘M_20_Valine’.
###
### -	MGS_abundance.tab: 
###   File with the abundance (e.g. median gene abundance) of
###   MGSs (columns) per individual (rows). These are assumed to have been
###   rarefied to comparable depth or otherwise normalized. For historical
###   reasons, the MGSs are labeled ‘T2DCAG’ followed by a number, e.g.
###   ‘T2DCAG00001’.
###
### -	KO_abundance.tab: 
###   File with the abundance of each KO (columns) per
###   individual (rows). The data is assumed to be rarefied to comparable depth or
###   otherwise normalized. For this, the software tool rtk58 can be used.
###
### - gene_abundance_sub.tab: 
###   File with the abundance of each catalog gene
###   (subset-version) in each individual assumed to be rarefied to comparable
###   depth or otherwise normalized.
### 
### Annotation-files:
### -	cluster_mapping_file.tab: 
###   Input file with annotation for metabolite
###   clusters, as available from curation of data in the specific dataset. The
###   WGCNA clustering algorithm names the generated clusters with color codes.
###   This mapping file simply facilitate renaming to more meaningful cluster
###   descriptions. Here, the serum polar metabolite and serum molecular lipid
###   clusters are labelled M01–M35 and L01–L39, respectively, and collectively
###   termed metabolite clusters.
###
### -	MGS_taxonomy.tab: 
###   File with taxonomic annotation of the MGSs (rows) used
###   in the analysis. Each row contains the following information:
###   species_taxonomy, species_pct, genus_taxonomy, genus_pct, family_taxonomy,
###   family_pct, order_taxonomy, order_pct, phylum_taxonomy, phylum_pct, where
###   x_pct is the percentage of the MGS genes that can be annotated (by sequence
###   similarity) to the taxonomy of the MGS. If no taxonomy can be assigned to
###   the MGS, the value will be NA.
###
### -	KEGG_modules.tab: 
###   File containing definition of KEGG gene functional
###   modules (rows) specifying which KO gene groups constitute each KEGG module.
###   The first two columns contain KEGG module entry (number) and name, the third
###   column list all KOs separated by semicolon. Any other functional annotation
###   used analogously could be swapped in instead of the name. This file can be
###   obtained by downloading KEGG modules from
###   www.genome.jp/kegg-bin/get_htext?ko00002.keg (www.genome.jp/kegg/ --> KEGG
###   MODULE --> KEGG modules); download htext and then running the provided
###   script ‘parse_kegg.pl’ (after changing input filename).
###
### -	KO_to_MGS.tab: 
###   File listing for each KO (rows) the MGSs (space-separated)
###   it is a member of. This is based on the KO annotation of the gene catalog
###   (gene_to_KO.tab) and the information, what gene from the gene catalog is
###   within each MGS as given in the list MGS_to_gene.tab.
###
### -	gene_to_KO.tab: 
###   File containing the KO annotation (if any) per gene (rows)
###   in the gene catalog (constituting 7,328,469 genes). Genes are labeled
###   ‘RefCat620’ followed by a number (1…7328469), e.g. ‘RefCat620.1’. Used to
###   create KO_to_MGS.tab. Note, for the purpose of this protocol and to
###   considerably reduce size of input data, only the subset of catalogue genes
###   with KO annotation (n = 2,205,769) are provided in files with gene abundance
###   or annotation (gene_to_KO.tab, MGS_to_gene.tab and gene_abundance_sub.tab)
###   as only those genes are used in the driver-species analysis.
###
### -	MGS_to_gene.tab: 
###   List of genes binned into a given MGS (rows). Used to
###   create KO_to_MGS.tab. Rather than gene names the file contains the index
###   value of the gene. i.e. position of the gene in the gene catalogue
###   (subset-version, thus 1…2205769).
###
### All demonstration files are tab-delimited text files, but other formats
### would equally work after modifying the respective file-import commands in
### the R-script. 
###
### For all demonstration data, pseudonymised sample names were re-randomised to
### generate anonymised data.
###

options (stringsAsFactors = FALSE)

### - phenotypes.tab

phenotypes = read.csv (file = "../../../Berberine2/Result/00.Phenotype/BerBer_phe_filter2.csv", row.names = 1, header = T)
phenotypes <- phenotypes[!is.na(phenotypes[, phenotype]), ]
ctrl.no.na <- rownames(phenotypes[phenotypes$Time==3 & phenotypes$Group==1, ])
### - metabolomic.tab

metabolomic = read.csv (file = "../../../Berberine2/Result/19.metabolism2/metabolism.pro.csv", row.names = 1, header = T,  check.names = F)

### - bileacid 
lipidomic = read.csv(file = "../../../Berberine2/Result/17.metabolism/Bileacid.metabolism.v2.csv", row.names = 1, header = T)

### - species 
mgs_abundance <- read.csv(file = "../../../Berberine2/Result/03.profile/Species/Berber.species.filter.csv", row.names = 1, header = T, check.names = F)
filter <- read.csv("../../../Berberine2/Result/14.Cutoff of genus and species/filter.species.list.csv")
mgs_abundance <- mgs_abundance[, as.character(filter[,1])]

### - KEGG_modules.tab
tmp = read.table ("../data/KEGG_modules.tab", sep = "\t")
koann = strsplit (tmp[,3], split = ";")
names (koann) = tmp[,1]

module_mapping = tmp[,2] ### description of Kegg modules
names (module_mapping) = tmp[,1] ; rm (tmp)

### remove KEGG references in square brackets for more clean names for plotting
module_mapping_clean = sapply(module_mapping, function(x) strsplit(x, " \\[")[[1]][1])

### - KO_abundance.tab
ko_abundance = read.table (file = "../../../Berberine2/Result/03.profile/Ko/Berber2.update.ko.pro", sep = "\t", row.names = 1, header = T)
ko_abundance <- t(ko_abundance)

### - humann2 /pathway 
pathway = read.table(file = "../../../Berberine2/Result/03.profile/05.humaan2/humann2.path.pro", header = T, row.names = 1, sep = "\t", quote = "",
                     na.strings = "asdfk")
pathway <- t(pathway)

### - output





