

mapply(
	function(x, y) which(x != y)[1], 
	strsplit(aba.genome[unique_name == "lpxC"]$target, ""), 
	strsplit(aba.genome[unique_name == "lpxC"]$spacer, "")) %>% 
	data.table(diff = .) %>% 
	cbind(aba.genome[unique_name == "lpxC"]) ->
	lpxC.stuff


lpxC.stuff <- lpxC.stuff %>% 
	mutate(
		beginning = paste0("'", substr(spacer, 1, diff - 1), "'"), 
		middle = paste0("bold(underline(", tolower(substr(spacer, diff, diff)), "))"), 
		end = case_when(
			diff == 20 ~ "",
			TRUE ~ paste0("*",substr(spacer, diff + 1, 20))), 
		recon = paste0(beginning ,"*" , middle, end)) %>%  
	mutate(spacer = case_when(
		type == "mismatch" ~ recon, 
		TRUE ~ paste0("bold(", spacer, ")"))) %>%
	mutate(y_pred = case_when(type == "perfect" ~ 1, TRUE ~ y_pred)) %>%
	mutate(diff = ifelse(is.na(diff), NA_real_, diff)) %>%
	mutate(diff = case_when(is.na(diff) ~ 0, TRUE ~ diff)) %>%
	group_by(target) %>%
	arrange(y_pred)

setorder(lpxC.stuff, offset, y_pred)

lpxC.stuff <- data.table(lpxC.stuff)

lpxC.stuff[, index := .I]

lpxC.stuff %>% 
	ggplot(aes(
		y = reorder(spacer, index), 
		x = y_pred)) + 
	geom_col(
		aes(fill = target),
		colour = "black") +
	doc_theme +
	scale_fill_viridis(discrete = T) +
	scale_y_discrete(labels = function(L) parse(text = L)) +
	# scale_linetype_manual(values = c("dashed", "solid")) +
	theme_bw() +
	theme(axis.text.y = element_text(size = 16, family = "Courier"),
		legend.position = "none") +
	xlab("Predicted Knockdown") + 
	ylab("Guide Sequence") +
	


# mutate(
# 	y_pred = case_when(is.na(y_pred) ~ 1, TRUE ~ y_pred), 
# 	spacer = fct_reorder(spacer, y_pred)) %>% 
