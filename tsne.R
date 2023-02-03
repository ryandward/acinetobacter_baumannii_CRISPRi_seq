# Load data into melted_results
aba_all <- fread(
  "all_counts_seal.tsv.gz",
  header = FALSE,
  col.names = c("spacer", "count", "condition"))

curated_names <- fread("curated_names.tsv")

aba_all <-
	aba_all %>% 
	group_by(condition) %>% 
	mutate(cpm = 10^6 * count / sum(count)) %>% 
	left_join(aba_key) %>% 
	left_join(
		curated_names %>% 
		rename("locus_tag" = "AB19606")) %>% 
	left_join(melted_results %>% select(unique_name, Pathway) %>% unique)

# Define the conditions to be compared to T0
aba_all_spread <- aba_all %>% 
	ungroup() %>% 
	filter(condition %in% (aba_design %>% filter(experiment == "tube") %>% pull(condition))) %>%
	data.table() %>%
	dcast(Pathway + spacer + type + unique_name ~ condition, value.var = "count", fill = 0)

# Load Rtsne library
library(Rtsne)

# Perform t-SNE
tsne_out <- aba_all_spread %>%
  select(-Pathway, -spacer, -type, -unique_name) %>%
  data.matrix %>%
  Rtsne()

# Combine t-SNE results with metadata
tsne_dt <- aba_all_spread %>%
  select(Pathway, spacer, type, unique_name) %>%
  cbind(tsne_out$Y %>% as_tibble)

# Join data frames and create a new column "Pathway"

# Plot the data
tsne_dt %>% 
	filter(Pathway != "Other" | is.na(Pathway))  %>% 
	ggplot(aes(x = V1, y = V2)) + 
  geom_point(aes(colour = Pathway), alpha = 0.5) +
  doc_theme +
  facet_grid( ~ type)

