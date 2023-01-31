# Load data into melted_results
melted_results <- fread("Results/melted_results.tsv.gz")

# Define the conditions to be compared to T0
tSNE_conditions <- c(
	"None_0_T1 - None_0_T0",
	"None_0_T2 - None_0_T0",
	"Colistin_0.44_T1 - None_0_T0",
	"Colistin_0.44_T2 - None_0_T0",
	"Rifampicin_0.34_T1 - None_0_T0",
	"Rifampicin_0.34_T2 - None_0_T0",
	"Imipenem_0.06_T1 - None_0_T0",
	"Imipenem_0.06_T2 - None_0_T0",
	"Imipenem_0.09_T1 - None_0_T0",
	"Imipenem_0.09_T2 - None_0_T0",
	"Meropenem_0.11_T1 - None_0_T0",
	"Meropenem_0.11_T2 - None_0_T0",
	"Meropenem_0.17_T1 - None_0_T0",
	"Meropenem_0.17_T2 - None_0_T0")

# Extract timing from conditions
tSNE_timing <- tSNE_conditions %>% stringr::str_extract("T[12]")

# Create metadata data.table
tSNE_meta <- data.table(condition = tSNE_conditions, timing = tSNE_timing)

# Join metadata with melted results
tSNE_results <- tSNE_meta %>% inner_join(melted_results)

# Get rid of timing information from the "shift" column
tSNE_results <- tSNE_results %>% mutate(shift = gsub("_T[12]", "", shift))

# Spread data for tSNE input
tSNE_results_spread <- tSNE_results %>% 
	data.table::dcast(timing + spacer + type + Pathway + unique_name ~ shift, value.var = "LFC")

# Load Rtsne library
library(Rtsne)

# Perform t-SNE
tsne_out <- tSNE_results_spread %>%
	select(-timing, -spacer, -type, -Pathway, -unique_name) %>%
	data.matrix %>%
	Rtsne(perplexity = 200)

# Combine t-SNE results with metadata
tsne_dt <- tSNE_results_spread %>%
	select(timing, spacer, type, Pathway, unique_name) %>%
	cbind(tsne_out$Y %>% as_tibble)

# Define a function for color scaling
scm = function(palette = cols) {
	scale_color_manual(values = palette, na.value = alpha("light gray", 0.25))
}

# Join data frames and create a new column "Pathway"
joined_dt <- inner_join(tsne_dt, aba_key) %>%
	mutate(
		Pathway = case_when(
			unique_name %like% "atp" ~ "ATP",
			unique_name %like% "cyo" ~ "CYO",
			unique_name %like% "nuo" ~ "NADH"))


# Plot the data
ggplot(joined_dt, aes(x = V1, y = V2)) + 
	geom_point(aes(colour = Pathway), alpha = 0.5) +
	scale_color_manual(values = c("red", "blue", "green")) +
	doc_theme +
	facet_grid(timing ~ type)

