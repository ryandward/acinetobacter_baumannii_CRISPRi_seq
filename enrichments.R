enrichments <- fread("Enrichments/summary.tsv")
# enrichments <- fread("Enrichments/induction_enrichment_2.tsv")


enrichments <- enrichments[
	, .(
		AB030 = unlist(strsplit(`matching proteins in your input (IDs)`, split = ",")))
		,
	by = .(
		condition, 
		`term ID`, 
		`term description`, 
		`enrichment score`, 
		`false discovery rate`, 
		direction, 
		`genes mapped`)]

enrichments[, AB030 := gsub("470.", "", AB030)]

enrichments <- 
	enrichments %>% 
	inner_join(
		aba.genome %>% 
			select(AB030, unique_name) %>% 
			unique)

enrichments %>%
	group_by(unique_name) %>% 
	filter(`enrichment score` == max(`enrichment score`)) %>% 
	arrange(desc(`enrichment score`)) %>% View

enrichments %>% 
	inner_join(fread("Enrichments/induction_enrichment_2.tsv") %>% group_by(`term ID`) %>% tally %>% filter(n==2) %>% ungroup %>% select(`term ID`)) %>%
	select(condition, `enrichment score`, `false discovery rate`, `term ID`, `term description`, direction)	%>%
	unique %>%
	filter(condition %like% "None_0_T0") %>%
	mutate(
		Pathway = case_when(
			`term description` %like% "Ribosome" ~ "Ribosome",
			`term description` %like% "Oxidative" | `term description` %like% "NADH" ~ "Respiration")) %>%
	mutate(
		`enrichment score` = case_when(
			direction == "bottom" ~ -`enrichment score`,
			direction == "top" ~ `enrichment score`
		)) %>%
	filter(!is.na(`enrichment score`)) %>%
	ggplot(aes(x = `enrichment score`, y = `false discovery rate`)) + 
	geom_point(aes(colour = Pathway)) + 
	facet_wrap(~condition) +
	doc_theme +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) 

+
	# geom_text_repel(
	# 	data = . %>% filter(`term ID` %in% c("CL:69")),
	# 	aes(label = `term description`),
	# 	min.segment.length = unit(0, 'lines'), 
	# 	# label.size = 0.15,
	# 	max.overlaps = Inf,
	# 	point.padding = 5,
	# 	# nudge_x = 1,
	# 	# nudge_y = 1,
	# 	segment.curvature = -1e-20,
	# 	arrow = arrow(length = unit(0.015, "npc")),
	# 	colour = "black")
