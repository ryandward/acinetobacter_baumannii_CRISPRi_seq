#https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0146021.s001

# source("../meta_counter.R")

library(pacman)
p_load(data.table, drc, protti, dplyr, magrittr, ggplot2)


this_gene <-  "murC"

this_shift_1 <- "None_0_T1"
this_shift_2 <- "Meropenem_0.11_T1"


selected_results <- 
	melted_results[
		# type == "mismatch" &
		y_pred > 0 &
			shift %in% c(this_shift_1, this_shift_2) &
			base == "None_0_T0" &
			unique_name %like% this_gene]

selected_results <- selected_results[target %in% selected_results[, .N, by = .(target, condition)][N>2, .N, by = .(target)][N==2]$target]

targets <- selected_results[, .N, by = .(target, condition)][N>2][, unique(target)]

min_response = min(selected_results$LFC)
max_response = max(selected_results$LFC)

plot(model_1, 
		 broken = F,
		 ylab = "Log 2-fold Change", 
		 xlab = "Predicted Guide Efficacy", 
		 col = "white", 
		 add = F,
		 ylim = c (min_response, max_response),
		 xlim = c (0, 1.2),
		 main =  ("nuo genes in colistin after five generations"))


for(i in targets){
	
	target_results <- selected_results[target == i]
	
target_results_1 <-
	aba_key[, .(target, offset)][target_results[shift == this_shift_1], on = .(target)]


model_1 <- drm(
	LFC ~ y_pred, 
	data = target_results_1, 
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
		 lwd = 10, 
		 cex = 0.5, 
		 cex.axis = 1.2, 
		 cex.lab = 1.2, 
		 pch = 20, 
		 col = alpha("blue", 0.15), 
		 add = T)

################################################################################

target_results_2 <-
	aba_key[, .(target, offset)][target_results[shift == this_shift_2], on = .(target)]


model_2 <- drm(
	LFC ~ y_pred, 
	data = target_results_2, 
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
		 lwd = 10, 
		 cex = 0.5, 
		 cex.axis = 1.2, 
		 cex.lab = 1.2, 
		 pch = 20, 
		 col = alpha("red", 0.15), 
		 add = T)

abline(
	h = 0, 
	lty = 'dotted')
################################################################################

legend("bottomleft", 
			 legend = c("Induction", "Induction + Colistin"),
			 col = c("blue", "red"), 
			 lwd = 2,
			 cex = 0.8)
}

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
		 lwd = 6, 
		 cex = 0, 
		 cex.axis = 1.2, 
		 cex.lab = 1.2, 
		 pch = 20, 
		 col = "blue", 
		 add = T)

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
		 lwd = 6, 
		 cex = 0, 
		 cex.axis = 1.2, 
		 cex.lab = 1.2, 
		 pch = 20, 
		 col = "red", 
		 add = T)
