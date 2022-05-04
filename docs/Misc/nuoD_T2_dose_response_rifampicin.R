#https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0146021.s001

# source("meta_counter.R")

library(pacman)
p_load(data.table, drc, protti, dplyr, magrittr, ggplot2)


this_gene <-  "nuoD"

this_shift_1 <- "None_0_T2"
this_shift_2 <- "Rifampicin_0.34_T2"


selected_results <- 
	melted_results[
		type == "mismatch" &
			y_pred > 0 &
			shift %in% c(this_shift_1, this_shift_2) &
			base == "None_0_T0" &
			unique_name %like% this_gene]


min_response = min(selected_results$LFC)
max_response = max(selected_results$LFC)


selected_results_1 <-
	aba_key[, .(target, offset)][selected_results[shift == this_shift_1], on = .(target)]


model_1 <- drm(
	LFC ~ y_pred, 
	data = selected_results_1, 
	fct = LL.4(
		names = c(
			"Slope", 
			"Lower Limit", 
			"Upper Limit", 
			"ED50")))


plot(model_1, 
		 broken = F,
		 ylab = "Log 2-fold Change", 
		 xlab = "Predicted Guide Efficacy", 
		 lwd = 2, 
		 cex = 1.2, 
		 cex.axis = 1.2, 
		 cex.lab = 1.2, 
		 pch = 20, 
		 col = "green", 
		 add = F,
		 ylim = c (min_response, max_response),
		 xlim = c (0, 1.2),
		 main =  ("nuoD in rifampicin after ten doublings."))

################################################################################

selected_results_2 <-
	aba_key[, .(target, offset)][selected_results[shift == this_shift_2], on = .(target)]


model_2 <- drm(
	LFC ~ y_pred, 
	data = selected_results_2, 
	fct = LL.4(
		names = c(
			"Slope", 
			"Lower Limit", 
			"Upper Limit", 
			"ED50")))


plot(model_2, 
		 broken = F,
		 ylab = "Log 2-fold Change", 
		 xlab = "Predicted Guide Efficacy", 
		 lwd = 2, 
		 cex = 1.2, 
		 cex.axis = 1.2, 
		 cex.lab = 1.2, 
		 pch = 20, 
		 col = "red", 
		 add = T,
		 ylim = c (min_response, max_response),
		 xlim = c (0, 1.2))

abline(
	h = 0, 
	lty = 'dotted')
################################################################################

legend("bottomleft", 
			 legend = c("Induction", "Induction + Rifampicin"),
			 col = c("green", "red"), 
			 lwd = 2,
			 cex = 0.8)
