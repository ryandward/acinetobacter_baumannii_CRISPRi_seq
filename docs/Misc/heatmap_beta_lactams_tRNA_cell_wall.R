# source('meta_counter.R')

selected_LFC <- 
	melted_results[
		((shift %like% "T1" & base %like% "T1") | (shift %like% "T2" & base %like% "T2"))
		& base %like% "None"
		& !shift %like% "None" 
		& !Pathway %like% "Ribosome" 
		& !Pathway %like% "Other" 
		& (shift %like% "Mero" | shift %like% "Imi")
		& !Pathway  %like% "Quinone"
		& !Pathway %like% "LPS"
		& type == "perfect"]


median_LFC <- dcast(
	selected_LFC, 
	AB030 + unique_name + Pathway ~ condition, 
	value.var = "LFC", 
	fun.aggregate = median)


LFC_grid <- data.matrix(
	median_LFC[, c(-1:-3)])

rownames(LFC_grid) <- median_LFC[, paste(unique_name)]

# format = Pathway + italic(Gene Name)
newnames <- lapply(
	median_LFC[, .I],
	function(x) bquote(.(median_LFC[x]$Pathway) ~ italic(.(median_LFC[x]$unique_name))))

plot_matrix <- LFC_grid

# colnames(plot_matrix) <- c(
# 	# "Colistin", 
# 	"Imipenem (low)", 
# 	"Imipenem (high)", 
# 	"Meropenem (low)", 
# 	"Meropenem (high)"
# 	# "Rifampicin"
# )

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

to_plot_title <- paste("Log Fold Change Fitness by Pathway")

to_plot <- pheatmap(plot_matrix,
										col = plot_colors,
										breaks = breaks,
										border_color = NA,
										cellwidth = 45,
										cellheight = 14,
										cutree_rows = 2,
										cutree_cols = 2,
										main = to_plot_title,
										angle_col = 45,
										fontsize_col = 16,
										fontsize_row = 12,
										clustering_method = "ward.D2",
										# clustering_distance_rows = "maximum",
										# clustering_distance_cols = "maximum",
										cluster_cols = FALSE,
										gaps_col = c(2,4,6,8),
										labels_row = as.expression(newnames))


print(to_plot)

