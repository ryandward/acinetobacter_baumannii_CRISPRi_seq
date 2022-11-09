# source('counter_in_progress.R')


library(data.table)
library(pheatmap)
library(poolr)
library(pracma)

# modify this number to the number of genes you want to follow up on
target_k = 25

# modify this to the pathway you want to follow up on
this_pathway = "Cell Wall/PG"


# change the conditions to match (approximately) the condition you are interested in
selected_LFC <- 
	melted_results[Pathway == this_pathway
		& (condition %like% "Mero" | condition %like% "Imi")
		& type == "perfect"]


LFC_guides <- dcast(
	selected_LFC, 
	AB030 + unique_name + spacer + Pathway ~ condition, 
	value.var = "LFC", 
	fun.aggregate = median)



LFC_grid <- data.matrix(
	LFC_guides[, c(-1:-4)])

rownames(LFC_grid) <- LFC_guides[, paste(spacer)]

# format = Pathway + italic(Gene Name)
newnames <- lapply(
	LFC_guides[, .I],
	function(x) bquote(.(LFC_guides[x]$Pathway) ~ italic(.(LFC_guides[x]$unique_name)) ~ .(LFC_guides[x]$spacer)))

plot_matrix <- LFC_grid


break_halves <- length(unique(as.vector(plot_matrix)))

breaks <- c(
	seq(min(plot_matrix, na.rm = T), 
			0, 
			length.out = break_halves)[-break_halves], 
	seq(0, 
			max(plot_matrix, na.rm = T), 
			length.out = break_halves))

plot_colors <- c(
	colorRampPalette(c("#ba000d", "white"))(break_halves)[-break_halves], 
	colorRampPalette(c("white", "#007ac1"))(break_halves)[-1])

to_plot_title <- paste("Log Fold Change Fitness by Pathway")

to_plot <- pheatmap(plot_matrix,
										col = plot_colors,
										breaks = breaks,
										border_color = NA,
										cellwidth = 45,
										cellheight = 14,
										# cutree_rows = 2,
										# cutree_cols = 2,
										main = to_plot_title,
										angle_col = 45,
										fontsize_col = 16,
										fontsize_row = 12,
										clustering_method = "ward.D2",
										clustering_distance_rows = "maximum",
										clustering_distance_cols = "maximum",
										labels_row = as.expression(newnames))


target_h <-
	sort(
		to_plot$tree_row$height, 
		decreasing = T)[target_k]


plot(
	to_plot$tree_row, 
	cex = 0.35)

abline(
	h = target_h, 
	col="red", 
	lty = 2, 
	lwd = 2)


spacer_performance <-
data.matrix(
	sort(
		cutree(
			to_plot$tree_row, 
			h = target_h)))

spacer_performance <-
	data.table(spacer_performance, keep.rownames = "spacer")

setnames(spacer_performance, "V1", "k_group")

stouffer_FDR <- selected_LFC[, .(FDR = stouffer(FDR)$p), by = .(AB030, unique_name, spacer)]

spacer_performance <-
	stouffer_FDR[spacer_performance, on = .(spacer)]

best_guides <-
spacer_performance[spacer_performance[, .(FDR = min(FDR)), by = .(k_group)], on = .(FDR, k_group)]

best_guides <- selected_LFC[, .(medLFC = median(LFC), meanLFC = mean(LFC), stdLFC = std(LFC)), by = .(AB030, unique_name, spacer)][best_guides, on = .(AB030, unique_name, spacer)]

best_guides <- best_guides[best_guides[, .(FDR = min(FDR)), by = .(unique_name)], on = .(FDR, unique_name)]

print(best_guides)

################################################################################
# Second heat map

LFC_guides <-
	LFC_guides[spacer %in% best_guides$spacer]


LFC_grid <- data.matrix(
	LFC_guides[, c(-1:-4)])

rownames(LFC_grid) <- LFC_guides[, paste(spacer)]

# format = Pathway + italic(Gene Name)
newnames <- lapply(
	LFC_guides[, .I],
	function(x) bquote(.(LFC_guides[x]$Pathway) ~ italic(.(LFC_guides[x]$unique_name)) ~ .(LFC_guides[x]$spacer)))

plot_matrix <- LFC_grid


break_halves <- length(unique(as.vector(plot_matrix)))

breaks <- c(
	seq(min(plot_matrix), 
			0, 
			length.out = break_halves)[-break_halves], 
	seq(0, 
			max(plot_matrix), 
			length.out = break_halves))

plot_colors <- c(
	colorRampPalette(c("#ba000d", "white"))(break_halves)[-break_halves], 
	colorRampPalette(c("white", "#007ac1"))(break_halves)[-1])

to_plot_title <- paste(paste(nrow(best_guides), "Most Consistent Diverse Guides for Pathway"))

to_plot <- pheatmap(plot_matrix,
										col = plot_colors,
										breaks = breaks,
										border_color = NA,
										cellwidth = 45,
										cellheight = 14,
										# cutree_rows = 2,
										# cutree_cols = 2,
										main = to_plot_title,
										angle_col = 45,
										# fontsize_col = 16,
										# fontsize_row = 12,
										clustering_method = "ward.D2",
										clustering_distance_rows = "maximum",
										clustering_distance_cols = "maximum",
										cluster_cols = FALSE,
										labels_row = as.expression(newnames))

