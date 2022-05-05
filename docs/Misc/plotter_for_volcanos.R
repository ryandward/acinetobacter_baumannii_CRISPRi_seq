library(pacman)
p_load(ggplot2, data.table, pheatmap, ggrepel, hrbrthemes, viridis)

melted_results <- fread("Results/melted_results.tsv.gz", sep = "\t")

interest <- fread("interest.tsv", sep = "\t")

#T1

to_plot <-
	median_melted_results[
		condition == "None_0_T1 - None_0_T0" & type == "perfect"]

plot_object <-
	ggplot(data = to_plot,
				 aes(x = medLFC,
				 		y = -log10(FDR))) +
	geom_point(aes(),
						 color = "grey") +
	geom_hline(yintercept = 1.30103,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	geom_text_repel(
		data = to_plot[unique_name == "GO593_00515"],
		aes(label = gene_name_stylized),
		size = 5,
		box.padding = unit(1, "lines"),
		point.padding = unit(0.25, "lines"),
		parse = TRUE) +
	theme_ipsum() +
	# ggtitle("Effect of Library Induction at T[1]")
	ggtitle(bquote(bold(Effect~of~Library~Induction~at~t[1])))


print(plot_object)

# T2

to_plot <-
	median_melted_results[
		condition == "None_0_T2 - None_0_T0" & type == "perfect"]

plot_object <-
	ggplot(data = to_plot,
				 aes(x = medLFC,
				 		y = -log10(FDR))) +
	geom_point(aes(),
						 color = "grey") +
	geom_hline(yintercept = 1.30103,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	geom_text_repel(
		data = to_plot[unique_name == "GO593_00515"],
		aes(label = gene_name_stylized),
		size = 5,
		box.padding = unit(1, "lines"),
		point.padding = unit(0.25, "lines"),
		parse = TRUE) +
	theme_ipsum() +
	# ggtitle("Effect of Library Induction at T[2]")
	ggtitle(bquote(bold(Effect~of~Library~Induction~at~t[2])))


print(plot_object)


##############################	##############################	#################

to_plot <-
	median_melted_results[
		condition == "Imipenem_0.09_T1 - None_0_T1" & type == "perfect"]

to_plot[Pathway %like% "Cell Wall", `Cell Wall Genes` := TRUE]

plot_object <-
	ggplot(data = to_plot,
				 aes(x = medLFC,
				 		y = -log10(FDR))) +
	geom_point(aes(color = `Cell Wall Genes`)) +
	geom_hline(yintercept = 1.30103,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	geom_text_repel(
		data = to_plot[`Cell Wall Genes` == TRUE],
		aes(label = gene_name_stylized),
		size = 5,
		box.padding = unit(0.5, "lines"),
		point.padding = unit(0.25, "lines"),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf) +
	theme_ipsum() +
	theme(legend.position = "bottom") +
	scale_color_manual(values = c("#33a02c"), na.value = "grey") +
	# ggtitle("Effect of Library Induction + Imipenem at T[1]")
	ggtitle(bquote(bold(Effect~of~Library~Induction~"+"~Imipenem~at~t[1])))


print(plot_object)

to_plot[Pathway %like% "tRNA", `tRNA Genes` := TRUE]

plot_object <-
	ggplot(data = to_plot,
				 aes(x = medLFC,
				 		y = -log10(FDR))) +
	geom_point(aes(color = `tRNA Genes`)) +
	geom_hline(yintercept = 1.30103,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	geom_text_repel(
		data = to_plot[`tRNA Genes` == TRUE],
		aes(label = gene_name_stylized),
		size = 5,
		box.padding = unit(0.5, "lines"),
		point.padding = unit(0.25, "lines"),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf) +
	theme_ipsum() +
	theme(legend.position = "bottom") +
	scale_color_manual(values = c("#6a3d9a"), na.value = "grey") +
	# ggtitle("Effect of Library Induction + Imipenem at T[1]")
	ggtitle(bquote(bold(Effect~of~Library~Induction~"+"~Imipenem~at~t[1])))


print(plot_object)

################################################################################
to_plot <-
	median_melted_results[
		condition == "Rifampicin_0.34_T2 - None_0_T2" & type == "perfect"]

to_plot[Pathway %like% "LPS", `LPS Genes` := TRUE]
to_plot[Pathway %like% "NADH", `NADH Genes` := TRUE]

plot_object <-
	ggplot(data = to_plot,
				 aes(x = medLFC,
				 		y = -log10(FDR))) +
	geom_point(aes(color = `LPS Genes`)) +
	geom_hline(yintercept = 1.30103,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	geom_text_repel(
		data = to_plot[`LPS Genes` == TRUE],
		aes(label = gene_name_stylized),
		size = 5,
		box.padding = unit(0.5, "lines"),
		point.padding = unit(0.25, "lines"),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf) +
	theme_ipsum() +
	theme(legend.position = "bottom") +
	scale_color_manual(values = c("#1f78b4"), na.value = "grey") +
	# ggtitle("Effect of Library Induction + Rifampicin at T[2]")
	ggtitle(bquote(bold(Effect~of~Library~Induction~"+"~Rifampicin~at~t[2])))

print(plot_object)

plot_object <-
	ggplot(data = to_plot,
				 aes(x = medLFC,
				 		y = -log10(FDR))) +
	geom_point(aes(color = `NADH Genes`)) +
	geom_hline(yintercept = 1.30103,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	geom_text_repel(
		data = to_plot[`NADH Genes` == TRUE],
		aes(label = gene_name_stylized),
		size = 5,
		box.padding = unit(0.5, "lines"),
		point.padding = unit(0.25, "lines"),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf) +
	theme_ipsum() +
	theme(legend.position = "bottom") +
	scale_color_manual(values = c("#ff7f00"), na.value = "grey") +
	# ggtitle("Effect of Library Induction + Rifampicin at T[2]")
	ggtitle(bquote(bold(Effect~of~Library~Induction~"+"~Rifampicin~at~t[2])))

print(plot_object)

################################################################################

to_plot <-
	median_melted_results[
		condition == "Colistin_0.44_T2 - None_0_T2" & type == "perfect"]

to_plot[Pathway %like% "LPS", `LPS Genes` := TRUE]
to_plot[Pathway %like% "NADH", `NADH Genes` := TRUE]

plot_object <-
	ggplot(data = to_plot,
				 aes(x = medLFC,
				 		y = -log10(FDR))) +
	geom_point(aes(color = `LPS Genes`)) +
	geom_hline(yintercept = 1.30103,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	geom_text_repel(
		data = to_plot[`LPS Genes` == TRUE],
		aes(label = gene_name_stylized),
		size = 5,
		box.padding = unit(0.5, "lines"),
		point.padding = unit(0.25, "lines"),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf) +
	theme_ipsum() +
	theme(legend.position = "bottom") +
	scale_color_manual(values = c("#1f78b4"), na.value = "grey") +
	# ggtitle("Effect of Library Induction + Colistin at T[2]")
	ggtitle(bquote(bold(Effect~of~Library~Induction~"+"~Colistin~at~t[2])))

print(plot_object)

plot_object <-
	ggplot(data = to_plot,
				 aes(x = medLFC,
				 		y = -log10(FDR))) +
	geom_point(aes(color = `NADH Genes`)) +
	geom_hline(yintercept = 1.30103,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	geom_text_repel(
		data = to_plot[`NADH Genes` == TRUE],
		aes(label = gene_name_stylized),
		size = 5,
		box.padding = unit(0.5, "lines"),
		point.padding = unit(0.25, "lines"),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf) +
	theme_ipsum() +
	theme(legend.position = "bottom") +
	scale_color_manual(values = c("#ff7f00"), na.value = "grey") +
	# ggtitle("Effect of Library Induction + Colistin at T[2]")
	ggtitle(bquote(bold(Effect~of~Library~Induction~"+"~Colistin~at~t[2])))


print(plot_object)

################################################################################

to_plot <-
	median_melted_results[
		condition == "Meropenem_0.17_T1 - None_0_T1" & type == "perfect"]

to_plot[Pathway %like% "Cell Wall", `Cell Wall Genes` := TRUE]

plot_object <-
	ggplot(data = to_plot,
				 aes(x = medLFC,
				 		y = -log10(FDR))) +
	geom_point(aes(color = `Cell Wall Genes`)) +
	geom_hline(yintercept = 1.30103,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	geom_text_repel(
		data = to_plot[`Cell Wall Genes` == TRUE],
		aes(label = gene_name_stylized),
		size = 5,
		box.padding = unit(0.5, "lines"),
		point.padding = unit(0.25, "lines"),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf) +
	theme_ipsum() +
	theme(legend.position = "bottom") +
	scale_color_manual(values = c("#33a02c"), na.value = "grey") +
	# ggtitle("Effect of Library Induction + Meropenem at T[1]")
	ggtitle(bquote(bold(Effect~of~Library~Induction~"+"~Meropenem~at~t[1])))


print(plot_object)

to_plot[Pathway %like% "tRNA", `tRNA Genes` := TRUE]

plot_object <-
	ggplot(data = to_plot,
				 aes(x = medLFC,
				 		y = -log10(FDR))) +
	geom_point(aes(color = `tRNA Genes`)) +
	geom_hline(yintercept = 1.30103,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	geom_text_repel(
		data = to_plot[`tRNA Genes` == TRUE],
		aes(label = gene_name_stylized),
		size = 5,
		box.padding = unit(0.5, "lines"),
		point.padding = unit(0.25, "lines"),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf) +
	theme_ipsum() +
	theme(legend.position = "bottom") +
	scale_color_manual(values = c("#6a3d9a"), na.value = "grey") +
	# ggtitle("Effect of Library Induction + Meropenem at T[1]")
	ggtitle(bquote(bold(Effect~of~Library~Induction~"+"~Meropenem~at~t[1])))


print(plot_object)

################################################################################

to_plot <-
	median_melted_results[
		condition == "Colistin_0.44_T1 - None_0_T1" & type == "perfect"]

to_plot[Pathway %like% "LPS", `LPS Genes` := TRUE]
to_plot[Pathway %like% "NADH", `NADH Genes` := TRUE]

plot_object <-
	ggplot(data = to_plot,
				 aes(x = medLFC,
				 		y = -log10(FDR))) +
	geom_point(aes(color = `LPS Genes`)) +
	geom_hline(yintercept = 1.30103,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	geom_text_repel(
		data = to_plot[`LPS Genes` == TRUE],
		aes(label = gene_name_stylized),
		size = 5,
		box.padding = unit(0.5, "lines"),
		point.padding = unit(0.25, "lines"),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf) +
	theme_ipsum() +
	theme(legend.position = "bottom") +
	scale_color_manual(values = c("#1f78b4"), na.value = "grey") +
	# ggtitle("Effect of Library Induction + Colistin at T[1]")
	ggtitle(bquote(bold(Effect~of~Library~Induction~"+"~Colistin~at~t[1])))

print(plot_object)

plot_object <-
	ggplot(data = to_plot,
				 aes(x = medLFC,
				 		y = -log10(FDR))) +
	geom_point(aes(color = `NADH Genes`)) +
	geom_hline(yintercept = 1.30103,
						 linetype = "dashed",
						 color = "red") +
	geom_vline(xintercept = 0,
						 linetype = "dotted",
						 color = "black") +
	geom_text_repel(
		data = to_plot[`NADH Genes` == TRUE],
		aes(label = gene_name_stylized),
		size = 5,
		box.padding = unit(0.5, "lines"),
		point.padding = unit(0.25, "lines"),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf) +
	theme_ipsum() +
	theme(legend.position = "bottom") +
	scale_color_manual(values = c("#ff7f00"), na.value = "grey") +
	# ggtitle("Effect of Library Induction + Colistin at T[1]")
	ggtitle(bquote(bold(Effect~of~Library~Induction~"+"~Colistin~at~t[1])))


print(plot_object)


