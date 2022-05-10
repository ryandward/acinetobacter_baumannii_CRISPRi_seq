require('pacman');

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
	vegan)

curated_names <- fread(
	"curated_names.tsv")

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
		"completeness"))

aba <- fread(
	"all_counts_seal.tsv.gz",
	header = FALSE,
	col.names = c(
		"spacer",
		"count",
		"condition"))

aba_key <- fread(
	"aba_key.tsv")

aba_design <- fread(
	"ABA1_experimental_design.tsv",
	na.strings = c("NA"))

aba_genome <- aba_bed[
	aba_key[, .(spacer, type, locus_tag, y_pred, target, offset)],
	on = .(locus_tag)]

# define the experimental design space to only take into consideration "tubes"
aba_design <- aba_design[experiment == "tube"]

publication_design <- copy(aba_design)

publication_design[, timing := gsub("T", "t", timing)]

publication_design[, rep := paste0("(", rep, ")")]

# keep only the counts that are in the experimental design space
aba <- aba[condition %in% aba_design$condition]

# convert single column into a table
aba_grid <-
	data.table::dcast(
		aba,
		spacer ~ factor(condition, levels = unique(condition)),
		value.var = "count",
		fill = 0)

aba_grid_matrix <-
	data.matrix(aba_grid[, -c("spacer")])

row.names(aba_grid_matrix) <- aba_grid$spacer

aba_group <-
	factor(
		aba_design[,  paste(drug, dose, timing, sep = "_")])

aba_permut <-
	model.matrix( ~ 0 + aba_group)

colnames(aba_permut) <-
	levels(aba_group)

aba_y <-
	DGEList(
		counts = aba_grid_matrix,
		group = aba_group,
		genes = row.names(aba_grid_matrix))

aba_keep <-
	filterByExpr(
		y = aba_y,
		design = aba_permut,
		group = aba_group)

aba_y <- aba_y[aba_keep, , keep.lib.sizes = FALSE]

aba_y <- calcNormFactors(aba_y)

aba_y <- estimateDisp(aba_y, aba_permut)

aba_fit <- glmQLFit(aba_y, aba_permut, robust = TRUE)

aba_CPM <- cpm(aba_y, prior.count = 0)

colnames(aba_CPM) <- factor(aba_design[,  paste(drug, dose, timing, sep = "_")])

########################

contrast_levels <- CJ(
	level2 = colnames(aba_permut),
	level1 = colnames(aba_permut))[
		level2 != level1,
		paste(level2, level1, sep = " - ")]

aba_contrast <- makeContrasts(
	contrasts = contrast_levels,
	levels = aba_permut)

########################

results_prelim <- glmQLFTest(aba_fit, contrast = aba_contrast)

results_prelim <- topTags(results_prelim, n = Inf)

results_prelim <- data.table(results_prelim$table)

results_prelim <- aba_key[results_prelim, on = .(spacer == genes)]

ineffective_guides <- results_prelim[type == "perfect" & FDR > 0.05, .(target)]

########################

results_FDR <- aba_key[, .(genes = unique(spacer))]

results_LFC <- aba_key[, .(genes = unique(spacer))]

########################

for (i in 1:ncol(aba_contrast)) {

	results <- glmQLFTest(aba_fit, contrast = aba_contrast[,i])

	results <- topTags(results, n = Inf)

	results <- data.table(results$table)

	print(paste("Processing results for", contrast_levels[i], "..."))

	results_FDR <- results[, .(genes, FDR)][results_FDR, on = .(genes)]

	setnames(results_FDR, "FDR", contrast_levels[i])

	results_LFC <- results[, .(genes, logFC)][results_LFC, on = .(genes)]

	setnames(results_LFC, "logFC", contrast_levels[i])

}

########################

results_FDR <- aba_genome[results_FDR, on = .(spacer == genes)]

results_FDR <- curated_names[, .(AB19606, AB030, unique_name)][results_FDR, on = .(AB19606 == locus_tag)]


results_LFC <- aba_genome[results_LFC, on = .(spacer == genes)]

results_LFC <- curated_names[, .(AB19606, AB030, unique_name)][results_LFC, on = .(AB19606 == locus_tag)]

########################

melted_results_FDR <-
	data.table::melt(
		results_FDR,
		id.vars = c(
			"AB19606",
			"AB030",
			"unique_name",
			"spacer",
			"type",
			"y_pred",
			"target",
			"offset"),
		variable.name = "condition",
		value.name = "FDR",
		measure.vars = contrast_levels)

melted_results_LFC <-
	data.table::melt(
		results_LFC,
		id.vars = c(
			"AB19606",
			"AB030",
			"unique_name",
			"spacer",
			"type",
			"y_pred",
			"target",
			"offset"),
		variable.name = "condition",
		value.name = "LFC",
		measure.vars = contrast_levels)

########################

melted_results <-
	melted_results_LFC[
		melted_results_FDR,
		on = .(
			AB19606,
			AB030,
			unique_name,
			spacer,
			type,
			y_pred,
			target,
			offset,
			condition)]

melted_results <- melted_results[!is.na(FDR) & !is.na(LFC)]

##########################################################################################

control_melted_results_by_condition <-
	melted_results[
		type == "control",
		.(med_LFC = median(LFC)),
		keyby = .(condition)]

setkey(melted_results, condition)

melted_results[, LFC.adj := control_melted_results_by_condition[
	melted_results, LFC - med_LFC, by = .EACHI]$V1]

##########################################################################################

median_melted_results <-
	melted_results[
		, .(medLFC = median(LFC.adj), FDR = stouffer(FDR)$p),
		by = .(AB19606, AB030, unique_name, type, condition)]

############################################################################################################

conditions <-
	c("None_0_T1 - None_0_T0",
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
		"Meropenem_0.17_T2 - None_0_T2")

interest <- data.table(condition = conditions)
############################################################################################################

results_LFC <- dcast(
	melted_results[condition %in% interest$condition],
	AB19606 + AB030 + unique_name + type + spacer + y_pred + target + offset ~ condition,
	value.var = "LFC")

results_FDR <- dcast(
	melted_results[condition %in% interest$condition],
	AB19606 + AB030 + unique_name + type + spacer + y_pred + target + offset ~ condition,
	value.var = "FDR")

median_results_LFC <- dcast(
	median_melted_results[condition %in% interest$condition],
	AB19606 + AB030 + unique_name + type ~ condition,
	value.var = "medLFC")

median_results_FDR <- dcast(
	median_melted_results[condition %in% interest$condition],
	AB19606 + AB030 + unique_name + type ~ condition,
	value.var = "FDR")

fwrite(results_LFC, "Results/results_LFC.tsv.gz", sep = "\t")

fwrite(results_FDR, "Results/results_FDR.tsv.gz", sep = "\t")

fwrite(median_results_LFC, "Results/median_results_LFC.tsv.gz", sep = "\t")

fwrite(median_results_FDR, "Results/median_results_FDR.tsv.gz", sep = "\t")

############################################################################################################
CL2682 <- fread('CL2687.tsv') #Cell wall biogenesis/degradation, and Cell Wall/PG

map03010 <- fread('map03010.tsv') #Ribosome, why is rplY and rpmB not being painted in plots?

LPS <- fread('LPS.tsv')
# CL:3059; Glycolipid metabolic process, and lipopolysaccharide transport

NADH <- fread('NADH.tsv')
# CL:852; NADH dehydrogenase activity

GO0004812 <- fread('GO0004812.tsv') #Aminoacyl-tRNA synthetase, GO

melted_results[AB030 %in% CL2682$AB030, Pathway := "Cell Wall/PG"]

melted_results[AB030 %in% map03010$AB030, Pathway := "Ribosome"]

melted_results[AB030 %in% LPS$AB030, Pathway := "LPS"]

melted_results[AB030 %in% NADH$AB030, Pathway := "NADH"]

melted_results[AB030 %in% GO0004812$AB030, Pathway := "tRNA Ligase"]

melted_results[is.na(Pathway), Pathway := "Other"]

melted_results[AB19606 == "", AB19606 := NA]

melted_results[, c("shift", "base") := tstrsplit(condition, " - ")]

melted_results[type == "perfect", y_pred := 1]

melted_results[, y_pred := y_pred - min(y_pred, na.rm = T)]


# melted_results[type != "control", lin_pred := (y_pred - min(y_pred) / (max(y_pred) - min(y_pred)))]

setorder(melted_results, FDR)

##########################################################################################


median_melted_results[AB030 %in% CL2682$AB030, Pathway := "Cell Wall/PG"]

median_melted_results[AB030 %in% map03010$AB030, Pathway := "Ribosome"]

median_melted_results[AB030 %in% LPS$AB030, Pathway := "LPS"]

median_melted_results[AB030 %in% NADH$AB030, Pathway := "NADH"]

median_melted_results[AB030 %in% GO0004812$AB030, Pathway := "tRNA Ligase"]

median_melted_results[is.na(Pathway), Pathway := "Other"]

median_melted_results[AB19606 == "", AB19606 := NA]

median_melted_results[, c("shift", "base") := tstrsplit(condition, " - ")]

setorder(median_melted_results, FDR)

##########################################################################################

median_melted_results[!unique_name %like% "GO593", gene_name_stylized := paste0("italic('", unique_name, "')")]

median_melted_results[unique_name %like% "GO593", gene_name_stylized := paste0("bold('", unique_name, "')")]


##########################################################################################

fwrite(melted_results, "Results/melted_results.tsv.gz", sep = "\t")

fwrite(median_melted_results, "Results/median_melted_results.tsv.gz", sep = "\t")

fwrite(interest, "interest.tsv", sep = "\t")
