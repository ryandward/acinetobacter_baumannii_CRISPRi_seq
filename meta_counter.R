# Load several packages from CRAN and Bioconductor
require('pacman')
p_load(
	data.table,
	scales,
	edgeR,
	statmod,
	poolr,
	pheatmap,
	svglite,
	ggplot2,
	ggrepel,
	Rtsne,
	pracma,
	colourpicker,
	RColorBrewer,
	vegan,
	tidyverse,
	magrittr
)

# Read in a file containing curated names
curated_names <- fread("curated_names.tsv")

# Read in a BED file
aba_bed <- fread(
	"CP046654.1.bed",
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

# Read in a file with counts and conditions
aba <- fread(
	"all_counts_seal.tsv.gz",
	header = FALSE,
	col.names = c("spacer", "count", "condition")
)

# Read in a file with spacer-related information
aba_key <- fread("aba_key.tsv")

# Read in a file with experimental design information
aba_design <- fread("ABA1_experimental_design.tsv", na.strings = "NA")

# Merge aba_bed and aba_key data frames
aba_genome <- aba_bed[aba_key[, .(spacer, type, locus_tag, y_pred, target, offset)], on = .(locus_tag)]

# define the experimental design space to only take into consideration "tubes"
aba_design <- aba_design[experiment == "tube"]

# Replace T with t and put parentheses around reps
publication_design <- copy(aba_design)[, c("timing", "rep") := .(gsub("T", "t", timing), paste0("(", rep, ")"))]


# keep only the counts that are in the experimental design space
aba <- aba[condition %in% aba_design$condition]

# convert single column into a table
# Convert data to a "wide" format with one column per condition
aba_grid <-
	data.table::dcast(
		aba,
		spacer ~ factor(condition, levels = unique(condition)),
		value.var = "count",
		fill = 0
	)

# Convert the data to a matrix, with spacer names as row names
aba_grid_matrix <- data.matrix(aba_grid[,-c("spacer")]) %>%
	set_rownames(aba_grid$spacer)

# Create a factor variable for the group labels
aba_group <-	factor(
	aba_design[,  paste(drug, dose, timing, sep = "_")],
	levels = unique(aba_design[,  paste(drug, dose, timing, sep = "_")]))

# Create a design matrix for the groups
aba_permut <- model.matrix(~ 0 + aba_group) %>%
	set_colnames(levels(aba_group))

# Create a DGEList object with the count data and group labels
aba_y <- DGEList(
	counts = aba_grid_matrix,
	group = aba_group,
	genes = row.names(aba_grid_matrix))

# Filter the DGEList object using the design matrix and group labels
# filterByExpr is part of the SummarizedExperiment package
aba_keep <-	aba_y %>% filterByExpr(design = aba_permut, group = aba_group)

# Subset the DGEList object
aba_y <- aba_y[aba_keep, , keep.lib.sizes = FALSE] 

# Normalize, and estimate dispersion
aba_y <- aba_y %>% calcNormFactors %>% estimateDisp(aba_permut)

##########################################################################################

# Fit a generalized linear model to the data
aba_fit <- aba_y %>% glmQLFit(aba_permut, robust = TRUE)

# Calculate counts per million (CPM) and set column names based on group labels
aba_CPM <- cpm(aba_y, prior.count = 0) %>%
	set_colnames(factor(aba_design[,  paste(drug, dose, timing, sep = "_")]))

##########################################################################################

# Create a data table called "contrast_levels" containing all pairs of levels for the
# "aba_permut" data table where the two levels are different from each other
contrast_levels <- CJ(
	level2 = colnames(aba_permut),
	level1 = colnames(aba_permut))[level2 != level1, paste(level2, level1, sep = " - ")]

# Create a contrast matrix called "aba_contrast" for the "aba_permut" data table using the
# pairs of levels contained in the "contrast_levels" data table
aba_contrast <- makeContrasts(contrasts = contrast_levels, levels = aba_permut)


##########################################################################################

# Calculate preliminary results to discard ineffective reads
results_prelim <- glmQLFTest(aba_fit, contrast = aba_contrast) %>% topTags(n = Inf) %>% 
	`$`(table) %>% data.table

# Join the results_prelim table with the aba_key table
results_prelim <- aba_key[results_prelim, on = .(spacer == genes)]

# Find guides that are ineffective in any condition
ineffective_guides <- results_prelim[type == "perfect" & FDR > 0.05, .(target)]

##########################################################################################

# Initialize data.table objects to store results
results_FDR <- aba_key[, .(genes = unique(spacer))]
results_LFC <- aba_key[, .(genes = unique(spacer))]

# Loop through each column in aba_contrast
for (i in 1:ncol(aba_contrast)) {
	
	# Perform a generalized linear model test using the contrast in aba_contrast[,i]
	results <- glmQLFTest(aba_fit, contrast = aba_contrast[, i])
	
	# Select the top tags from the results
	results <- topTags(results, n = Inf)
	
	# Convert the results to a data.table object
	results <- data.table(results$table)
	
	# Print a message indicating which contrast is being processed
	print(paste("Processing results for", contrast_levels[i], "..."))
	
	# Merge the results with the results_FDR object and rename the "FDR" column
	results_FDR <- results[, .(genes, FDR)][results_FDR, on = .(genes)]
	setnames(results_FDR, "FDR", contrast_levels[i])
	
	# Merge the results with the results_LFC object and rename the "logFC" column
	results_LFC <-
		results[, .(genes, logFC)][results_LFC, on = .(genes)]
	setnames(results_LFC, "logFC", contrast_levels[i])
}

##########################################################################################

# Join the results_FDR data frame with the aba_genome data frame on the spacer and genes columns
results_FDR <- aba_genome[results_FDR, on = .(spacer == genes)]

# Join the curated_names data frame with the results_FDR data frame on the AB19606 and locus_tag columns
# Then select the AB19606, AB030, and unique_name columns from the resulting data frame
results_FDR <- curated_names[, .(AB19606, AB030, unique_name)][results_FDR, on = .(AB19606 == locus_tag)]

# Join the results_LFC data frame with the aba_genome data frame on the spacer and genes columns
results_LFC <- aba_genome[results_LFC, on = .(spacer == genes)]

# Join the curated_names data frame with the results_LFC data frame on the AB19606 and locus_tag columns
# Then select the AB19606, AB030, and unique_name columns from the resulting data frame
results_LFC <- curated_names[, .(AB19606, AB030, unique_name)][results_LFC, on = .(AB19606 == locus_tag)]

##########################################################################################

# Melt the results_FDR data.table object and store the result in melted_results_FDR
melted_results_FDR <- results_FDR %>%
	data.table::melt(
		id.vars = c(
			"AB19606",
			"AB030",
			"unique_name",
			"spacer",
			"type",
			"y_pred",
			"target",
			"offset"
		),
		variable.name = "condition",
		value.name = "FDR",
		measure.vars = contrast_levels
	)

# Melt the results_LFC data.table object and store the result in melted_results_LFC
melted_results_LFC <- results_LFC %>%
	data.table::melt(
		id.vars = c(
			"AB19606",
			"AB030",
			"unique_name",
			"spacer",
			"type",
			"y_pred",
			"target",
			"offset"
		),
		variable.name = "condition",
		value.name = "LFC",
		measure.vars = contrast_levels
	)

# Join melted_results_LFC with melted_results_FDR using the "AB19606", "AB030",
# "unique_name", "spacer", "type", "y_pred", "target", "offset", and "condition" columns
melted_results <- melted_results_LFC %>%
	melted_results_FDR[, on = c(
		"AB19606",
		"AB030",
		"unique_name",
		"spacer",
		"type",
		"y_pred",
		"target",
		"offset",
		"condition"
	)]

# Filter melted_results to only include rows where the "FDR" and "LFC" columns are not NA
melted_results <- melted_results[!is.na(FDR) & !is.na(LFC)]

##########################################################################################

# Set the key for the melted_results table to the "condition" column
setkey(melted_results, condition)

# Calculate the LFC.adj column by subtracting the median value of the control sgRNA
# spacers from the LFC values for all sgRNA spacers in each condition
melted_results[, LFC.adj := LFC - median(LFC[type == "control"]), by = condition]

##########################################################################################

# Calculate the "medLFC" and "FDR" columns in a new data.table object using the "AB19606",
# "AB030", "unique_name", "type", and "condition" columns
median_melted_results <- melted_results[, .(medLFC = median(LFC.adj), FDR = stouffer(FDR)$p), by = .(AB19606, AB030, unique_name, type, condition)]

##########################################################################################

# Create a character vector containing a list of conditions
conditions <-
	c(
		"None_0_T1 - None_0_T0",
		"None_0_T2 - None_0_T0",
		"Colistin_0.44_T1 - None_0_T1",
		"Colistin_0.44_T2 - None_0_T2",
		"Rifampicin_0.34_T1 - None_0_T1",
		"Rifampicin_0.34_T2 - None_0_T2",
		"Imipenem_0.06_T1 - None_0_T1",
		"Imipenem_0.06_T2 - None_0_T2",
		"Imipenem_0.09_T1 - None_0_T1",
		"Imipenem_0.09_T2 - None_0_T2",
		"Meropenem_0.11_T1 - None_0_T1",
		"Meropenem_0.11_T2 - None_0_T2",
		"Meropenem_0.17_T1 - None_0_T1",
		"Meropenem_0.17_T2 - None_0_T2"
	)

# Create a data table with a single column called "condition" and assign the list of conditions to this column
interest <- data.table(condition = conditions)

##########################################################################################

# Create results_LFC and results_FDR data.tables
results_LFC <-
	dcast(
		melted_results[condition %in% interest$condition],
		AB19606 + AB030 + unique_name + type + spacer + y_pred + target + offset ~ condition,
		value.var = "LFC"
	)
results_FDR <-
	dcast(
		melted_results[condition %in% interest$condition],
		AB19606 + AB030 + unique_name + type + spacer + y_pred + target + offset ~ condition,
		value.var = "FDR"
	)

# Create median_results_LFC and median_results_FDR data.tables
median_results_LFC <-
	dcast(median_melted_results[condition %in% interest$condition],
				AB19606 + AB030 + unique_name + type ~ condition,
				value.var = "medLFC")
median_results_FDR <-
	dcast(median_melted_results[condition %in% interest$condition],
				AB19606 + AB030 + unique_name + type ~ condition,
				value.var = "FDR")

# Write data.tables to files
fwrite(results_LFC, "Results/results_LFC.tsv.gz", sep = "\t")
fwrite(results_FDR, "Results/results_FDR.tsv.gz", sep = "\t")
fwrite(median_results_LFC, "Results/median_results_LFC.tsv.gz", sep = "\t")
fwrite(median_results_FDR, "Results/median_results_FDR.tsv.gz", sep = "\t")

##########################################################################################

CellWall <- fread('CL704.tsv') #Cell wall biogenesis/degradation, and Cell Wall/PG
map03010 <- fread('map03010.tsv') #Ribosome, why is rplY and rpmB not being painted in plots?
LOS <- fread('LOS.tsv') # CL:3059; Glycolipid metabolic process, and lipopolysaccharide transport
NADH <- fread('NADH.tsv') # CL:852; NADH dehydrogenase activity
GO0004812 <- fread('GO0004812.tsv') #Aminoacyl-tRNA synthetase, GO

##########################################################################################

melted_results[, Pathway := case_when(
	AB030 %in% CellWall$AB030 ~ "Cell Wall/PG",
	AB030 %in% map03010$AB030 ~ "Ribosome",
	AB030 %in% LOS$AB030 ~ "LOS",
	AB030 %in% NADH$AB030 ~ "NADH",
	AB030 %in% GO0004812$AB030 ~ "tRNA Ligase",
	TRUE ~ "Other"
)]

melted_results[, AB19606 := ifelse(AB19606 == "", NA, AB19606)]
melted_results[, c("shift", "base") := tstrsplit(condition, " - ")]

melted_results[type == "perfect", y_pred := 1]
melted_results[, y_pred := y_pred - min(y_pred, na.rm = T)]


# melted_results[type != "control", lin_pred := (y_pred - min(y_pred) / (max(y_pred) - min(y_pred)))]

setorder(melted_results, FDR)

##########################################################################################

median_melted_results[, Pathway := case_when(
	AB030 %in% CellWall$AB030 ~ "Cell Wall/PG",
	AB030 %in% map03010$AB030 ~ "Ribosome",
	AB030 %in% LOS$AB030 ~ "LOS",
	AB030 %in% NADH$AB030 ~ "NADH",
	AB030 %in% GO0004812$AB030 ~ "tRNA Ligase",
	TRUE ~ "Other"
)]

median_melted_results[, AB19606 := ifelse(AB19606 == "", NA, AB19606)]
median_melted_results[, c("shift", "base") := tstrsplit(condition, " - ")]

setorder(median_melted_results, FDR)

##########################################################################################

median_melted_results[, gene_name_stylized := ifelse(
	unique_name %like% "GO593",
	unique_name,
	paste0("italic('", unique_name, "')"))]

##########################################################################################

fwrite(melted_results, "Results/melted_results.tsv.gz", sep = "\t")
fwrite(median_melted_results, "Results/median_melted_results.tsv.gz", sep = "\t")
fwrite(interest, "interest.tsv", sep = "\t")