orthos <- fread(
	"/home/ryandward/Downloads/linear_orthos.tsv", 
	header = FALSE, 
	col.names = c(
		"orthogroup",
		"strain",
		"locus_tag"
	))


orthos %>% 
	filter(strain %like% "19606") %>% 
	left_join(orthos %>% filter(strain %like% "AB030"), by = "orthogroup") %>% 
	left_join(curated_names, by = c("locus_tag.x" = "AB19606")) %>% 
	select(locus_tag.x, locus_tag.y) %>% 
	rename(AB19606 = locus_tag.x, AB030 = locus_tag.y) %>% 
	left_join(curated_names, by = "AB19606") %>% 
	rename(AB030_full = AB030.x, AB030_essentials = AB030.y) %>% 
	fwrite("curated_list_expanded.tsv", sep = "\t")


#######

library(data.table)
library(tidyverse)
library(ggrepel)
library(hrbrthemes)
library(ggallin)

doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")

orthos <- fread("N0.tsv")

orthos <- melt(
	orthos, 
	id.vars = c("HOG", "OG", "Gene Tree Parent Clade"), 
	variable.name = "genome", 
	value.name = "orthologs") %>%
	filter(!genome == "Caulobacter_vibrioides") %>%
	mutate(orthologs = strsplit(as.character(orthologs), ", ")) %>% 
	unnest(cols = "orthologs") %>%
	data.table %>%
	rename(locus_tag = orthologs)

curated_list_full <- orthos %>% 
	select(HOG, genome, locus_tag) %>% 
	filter(genome == "Acinetobacter_baumannii_19606") %>% 
	left_join(orthos %>% select(HOG, genome, locus_tag) %>% filter(genome == "Acinetobacter_baumannii_AB030"), by = "HOG") %>% 
	rename(AB19606 = locus_tag.x, AB030 = locus_tag.y) %>% 
	select(-genome.x, -genome.y) %>% 
	left_join(curated_names %>% select(-AB030), by = "AB19606") %>%
	na_if("unknown")

GCF_009759695 <- fread(
	"/home/ryandward/R/acinetobacter_baumannii_CRISPRi_seq/GCF_009759695.1.bed",
	header = FALSE,
	col.names = c(
		"chromosome",
		"left",
		"right",
		"locus_tag",
		"gene_name",
		"strand",
		"coding",
		"completeness"
	)
)

curated_list_full <- curated_list_full %>% 
	left_join(GCF_009759695 %>% rename(AB19606 = locus_tag)) %>% 
	na_if(".") %>% 
	mutate(
		genes = case_when(
			is.na(genes) ~ gene_name, 
			!is.na(genes) & !is.na(gene_name) & !(str_detect(gene_name, genes)) ~ paste(genes, gene_name, sep = ", "), 
			!is.na(genes) ~ genes))

string_info <- fread("470.protein.info.v11.5.txt")

curated_list_full <- curated_list_full %>% 
	left_join(
		string_info %>% rename("string_protein_id" = "#string_protein_id") %>% 
			mutate(string_protein_id = gsub("^470.", "", string_protein_id)), 
		by = c("AB030" = "string_protein_id")) %>% 
	select(HOG, genes, AB19606, AB030, preferred_name, unique_name, annotation, citation) %>% 
	unique %>% 
	mutate(unique_name = case_when(is.na(unique_name) ~ preferred_name, !is.na(unique_name) ~ unique_name)) %>% 
	rename(
		curated_name = unique_name,
		AB030_annotation = annotation) 

curated_list_full %>% fwrite("curated_list_full.tsv", sep = "\t")

