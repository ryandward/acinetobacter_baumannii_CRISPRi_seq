gene.guide.tally <- melted_results %>% 
	select(spacer, unique_name) %>% 
	unique %>% 
	group_by(unique_name) %>% 
	tally %>% 
	ungroup %>% 
	filter(unique_name != "")  %>% 
	select(n) 

gene.guide.tally %>% 
	ggplot() + 
	geom_bar(aes(n)) + 
	scale_x_continuous(
		minor_breaks = seq(1,14), 
		breaks = seq(1,14), 
		limits = c(1,15)) + 
	geom_text(
		data = . %>% group_by(n) %>% tally, 
		aes(y = nn + 15, x = n, label = nn)) + 
	xlab("Guides per gene") +
	ylab("Number of genes") +
	doc_theme
