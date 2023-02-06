# Load data into melted_results
aba_all <- fread(
	"all_counts_seal.tsv.gz",
	header = FALSE,
	col.names = c("spacer", "count", "condition")
)

curated_names <- fread("curated_names.tsv")

aba_all <- aba_all %>% 
	group_by(condition) %>% 
	mutate(cpm = 10^6 * count / sum(count)) %>% 
	left_join(aba_key) %>% 
	left_join(
		curated_names %>% 
			rename("locus_tag" = "AB19606")
	) %>% 
	left_join(melted_results %>% select(unique_name, Pathway) %>% unique)

aba_multi_normalized <- aba_all %>% 
	select(spacer, count, locus_tag, type) %>% 
	inner_join(aba_design %>% filter(experiment == "tube")) %>% 
	left_join(curated_names %>% rename("locus_tag" = "AB19606")) %>% 
	select(spacer, count, locus_tag, timing, drug, dose, rep, type, unique_name) %>% 
	filter(drug %in% c("None", "Colistin", "Imipenem", "Rifampicin", "Meropenem")) %>%
	group_by(condition) %>% 
	mutate(cpm = 10^6*count/sum(count)) %>% 
	mutate(cpm.rel = cpm/mean(cpm[type == "control"])) %>% 
	ungroup %>% group_by(spacer) %>% 
	mutate(logfc = log2(cpm.rel/mean(cpm.rel[timing == "T0"]))) %>% 
	ungroup %>% group_by(condition) %>% 
	mutate(
		logfc.adj = logfc - mean(logfc[type == "control"]),
		fc.adj = 2^logfc.adj,
		treatment = case_when(drug == "None" ~ "none", TRUE ~ "drug"),
		treatment = factor(treatment, levels = c("none", "drug"))
	) %>% 
	select(drug, dose, spacer, timing, rep, type, unique_name, treatment, fc.adj)

# Add T0 data to all combinations of drug, dose, and treatment
# T0 data serves as the baseline for comparison in subsequent plots
aba_multi_normalized <- crossing(
	aba_multi_normalized %>% ungroup %>% select(drug, dose, treatment) %>% unique,
	aba_multi_normalized %>% ungroup %>% filter(timing == "T0") %>% select(-drug, -dose, -treatment) %>% unique) %>%
	rbind(aba_multi_normalized) %>% unique

# Add data for "None" treatment to all combinations of drug and dose
# This allows for comparison with the "None" treatment in every plot
aba_multi_normalized <- crossing(
	aba_multi_normalized %>% ungroup %>% select(drug, dose) %>% unique,
	aba_multi_normalized %>% ungroup %>% filter(drug == "None") %>% select(-drug, -dose) %>% unique) %>%
	rbind(aba_multi_normalized) %>% unique

# Remove the "None" drug, as it is now represented in every plot
aba_multi_normalized <- aba_multi_normalized %>% filter(drug != "None")

# Reshape the data for easier handling of missing values
# Reshaping the data into a more manageable format makes it easier to fill in any missing values
aba_multi_normalized <- aba_multi_normalized %>% 
	data.table %>%
	dcast(drug + dose + treatment + spacer + unique_name ~ timing + rep, value.var = "fc.adj", fill = 0) %>%
	melt(id.vars = c("drug", "dose", "treatment", "spacer", "unique_name"), variable.name = "rep", value.name = "fc.adj") %>%
	separate(rep, c("timing", "rep"))


aba_multi_normalized %>%	
	filter(
		spacer %in% c(
			"ACAATTTGACGATCTTGCAA", 
			"AGAATCATTGCCGCAAGTAA", 
			"CCAGCACTACCATCCATAAT", 
			"TTGCTGCGCAGAATCGACAG", 
			"AGCTTGAAAGTTTGCATGTA",
			"GCATATTATCGGCATACGGA",
			"ATAACAGCAGCTTTTAAAAT")) %>%
	ggplot(aes(x = timing, y = fc.adj, group = interaction(unique_name, treatment))) +
	geom_point(aes(colour = unique_name, alpha = treatment)) +
	stat_summary(
		aes(linetype = treatment, alpha = treatment),
		fun.data = "mean_cl_boot",
		geom = "line",
		lwd = 0.75
	) +
	facet_grid(unique_name ~ drug + dose, scales = "free_y") +
	scale_colour_manual(
		values = c(
			"control" = "dark grey",
			"nuoB" = "#6A3D9A",
			"lpxC" = "#33A02C",
			"murA" = "#FF7F00",
			"glnS" = "#1F78B4",
			"sdhB" = "#12dfda",
			"atpG" = "#ffaffa",
			"psd" =  "#5555ff"
		)
	) +
	scale_linetype_manual(
		values = c(
			"none" = "solid",
			"drug" = "twodash"
		)
	) +
	scale_alpha_manual(
		values = c(
			"drug" = 1,
			"none" = 0.25
		)
	) + 
	doc_theme

