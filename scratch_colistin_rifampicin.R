
median_melted_results %>% 
	filter(type == "perfect") %>%
	mutate(
		Response = case_when(
			medLFC < -8 & FDR < 0.05 ~ "Very Vulnerable (LFC < -8)",
			medLFC < -1 & FDR < 0.05 ~ "Vulnerable (LFC < -1)",
			medLFC > 1 & FDR < 0.05 ~ "Resistant (LFC > 1)",
			TRUE ~ "No Response")) %>%
	select(AB19606, condition, Response) %>% 
	arrange(desc("Response")) %>%
	mutate(condition = case_when(
		condition == "Colistin_0.44_T1 - None_0_T1" ~ "Significant at t1",
		condition == "Colistin_0.44_T2 - None_0_T2" ~ "Significant at t2"))	%>%
	filter(!is.na(condition)) %>%
	pivot_longer(!c(condition, AB19606)) %>% 
	pivot_wider(id_cols = c(AB19606), names_from = condition, values_from = value) %>%
	make_long(`Significant at t1`, `Significant at t2`) %>% 
	ggplot(aes(
		x = x, 
		next_x = next_x,
		node = node,
		next_node = next_node,
		fill = factor(node),
		label = node)) +
	geom_sankey(
		flow.colour = "black",
		flow.alpha = 0.25, 
		node.color = "black") +
	scale_fill_manual(values = c(
		"grey", "#039be5", "#e53935", "#ffb300")) +
	geom_sankey_label(aes(colour = "node"),
										size = 3.5, color = 1) +
	theme_sankey(base_size = 16) +
	guides(
		fill = guide_legend(title = "Relative Response to Knockdown")) +
	theme(
		axis.title.x = element_blank(),
		legend.position = "bottom") +
	ggtitle("Essential Fitness Phenotypes in Colistin (beyond induction)") +
	theme(legend.position = "none")


################################################################################

median_melted_results %>% 
	filter(condition %in% c("Colistin_0.44_T1 - None_0_T1", "Colistin_0.44_T2 - None_0_T2")) %>%
	filter(type == "perfect") %>%
	mutate(
		Response = case_when(
			medLFC < -8 & FDR < 0.05 ~ "Very Vulnerable (LFC < -8)",
			medLFC < -1 & FDR < 0.05 ~ "Vulnerable (LFC < -1)",
			medLFC > 1 & FDR < 0.05 ~ "Resistant (LFC > 1)",
			TRUE ~ "No Response")) %>%
	inner_join(z) %>% 
	select(AB030, unique_name, condition, Response, medLFC, FDR, vuln.est, vuln.kd_50, vuln.p) %>% 
	filter(Response %like% "Very")

################################################################################
################################################################################

median_melted_results %>% 
	filter(type == "perfect") %>%
	mutate(
		Response = case_when(
			medLFC < -8 & FDR < 0.05 ~ "Very Vulnerable (LFC < -8)",
			medLFC < -1 & FDR < 0.05 ~ "Vulnerable (LFC < -1)",
			medLFC > 1 & FDR < 0.05 ~ "Resistant (LFC > 1)",
			TRUE ~ "No Response")) %>%
	select(AB19606, condition, Response) %>% 
	arrange(desc("Response")) %>%
	mutate(condition = case_when(
		condition == "Rifampicin_0.34_T1 - None_0_T1" ~ "Significant at t1",
		condition == "Rifampicin_0.34_T2 - None_0_T2" ~ "Significant at t2"))	%>%
	filter(!is.na(condition)) %>%
	pivot_longer(!c(condition, AB19606)) %>% 
	pivot_wider(id_cols = c(AB19606), names_from = condition, values_from = value) %>%
	make_long(`Significant at t1`, `Significant at t2`) %>% 
	ggplot(aes(
		x = x, 
		next_x = next_x,
		node = node,
		next_node = next_node,
		fill = factor(node),
		label = node)) +
	geom_sankey(
		flow.colour = "black",
		flow.alpha = 0.25, 
		node.color = "black") +
	scale_fill_manual(values = c(
		"grey", "#039be5", "#ffb300")) +
	geom_sankey_label(aes(colour = "node"),
										size = 3.5, color = 1) +
	theme_sankey(base_size = 16) +
	guides(
		fill = guide_legend(title = "Relative Response to Knockdown")) +
	theme(
		axis.title.x = element_blank(),
		legend.position = "bottom") +
	ggtitle("Essential Fitness Phenotypes in Rifampicin (beyond induction)") +
	theme(legend.position = "none")

################################################################################
################################################################################

median_melted_results %>% 
	filter(type == "perfect") %>%
	mutate(
		Response = case_when(
			medLFC < -8 & FDR < 0.05 ~ "Very Vulnerable (LFC < -8)",
			medLFC < -1 & FDR < 0.05 ~ "Vulnerable (LFC < -1)",
			medLFC > 1 & FDR < 0.05 ~ "Resistant (LFC > 1)",
			TRUE ~ "No Response")) %>%
	select(AB19606, condition, Response) %>% 
	arrange(desc("Response")) %>%
	mutate(condition = case_when(
		condition == "Colistin_0.44_T2 - None_0_T0" ~ "Colistin at t2",
		condition == "Rifampicin_0.34_T2 - None_0_T0" ~ "Rifampicin at t2"))	%>%
	filter(!is.na(condition)) %>%
	pivot_longer(!c(condition, AB19606)) %>% 
	pivot_wider(id_cols = c(AB19606), names_from = condition, values_from = value) %>%
	make_long(`Colistin at t2`, `Rifampicin at t2`) %>% 
	ggplot(aes(
		x = x, 
		next_x = next_x,
		node = node,
		next_node = next_node,
		fill = factor(node),
		label = node)) +
	geom_sankey(
		flow.colour = "black",
		flow.alpha = 0.25, 
		node.color = "black") +
	scale_fill_manual(values = c(
		"grey", "#039be5", "#e53935", "#ffb300")) +
	geom_sankey_label(aes(colour = "node"),
										size = 3.5, color = 1) +
	theme_sankey(base_size = 16) +
	guides(
		fill = guide_legend(title = "Relative Response to Knockdown")) +
	theme(
		axis.title.x = element_blank(),
		legend.position = "bottom") +
	ggtitle("Essential Fitness Phenotypes in Rifampicin/Colistin (beyond induction)") +
	theme(legend.position = "none")
