median_melted_results %>% filter(medLFC < -1 & FDR < 0.05) %>% filter(condition %in% c("None_0_T1 - None_0_T0")) %>% filter(type == "perfect") -> T1_stuff

median_melted_results %>% filter(medLFC < -1 & FDR < 0.05) %>% filter(condition %in% c("None_0_T1 - None_0_T0") & unique_name %in% T1_stuff$unique_name) %>% filter(type == "perfect") %>% pull(medLFC) %>% density -> T1_shape

median_melted_results %>% filter(medLFC < -1 & FDR < 0.05) %>% filter(condition %in% c("None_0_T2 - None_0_T0") & unique_name %in% T1_stuff$unique_name) %>% filter(type == "perfect") %>% pull(medLFC) %>% density -> T2_shape

plot(T2_shape, col = "red", ylim = c(0, 1), xlim = c(-15, 0))

lines(T1_shape, col = "blue")


###########

# genes with phenotypes
median_melted_results %>% filter(condition %in% (interest %>% mutate(condition = gsub("- None_0_T1", "- None_0_T0", condition)) %>% mutate(condition = gsub("- None_0_T2", "- None_0_T0", condition)) %>% pull("condition"))) %>% filter(FDR < 0.05) %>% filter(abs(medLFC) > 1) %>% filter(type != "perfect") %>% select(unique_name) %>% unique

