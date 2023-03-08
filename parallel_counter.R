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
	magrittr,
	parallel,
	doParallel,
	foreach
)

my.cluster <- parallel::makeCluster(
	parallel::detectCores() - 1, 
	type = "FORK"
)

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# Loop through each column in aba_contrast
results <- foreach (i = 1:1:ncol(aba_contrast)) %dopar% {
	
	glmQLFTest(aba_fit, contrast = aba_contrast[, i]) %>% 
		topTags(n = Inf) %>% use_series(table) %>% 
		mutate(
			condition = colnames(aba_contrast)[i],
			logFC = round(logFC, 2)) %>% 
		rename(spacer = genes) %>%
		data.table
}

parallel::stopCluster(cl = my.cluster)

if(length(x) > 0) {
	melted_results <- rbindlist(results)
}

melted_results <- melted_results %>% inner_join(aba_genome %>% select(spacer, type)) %>% data.table

melted_results[, logFC.adj := logFC - median(logFC[type == "control"]), by = condition]
