# Brain-Cousens Hormetic DRC Logistic Model (Not LLogistic)
# Ryan Ward
# Tuesday, April 11, 2023

BC.5.logistic <-
	function (fixed = c(NA, NA, NA, NA, NA),
						names = c("b", "c", "d", "e", "f"),
						...)
	{
		numParm <- 5
		if (!is.character(names) | !(length(names) == numParm)) {
			stop("Not correct 'names' argument")
		}
		if (!(length(fixed) == numParm)) {
			stop("Not correct length of 'fixed' argument")
		}
		return(braincousens_logistic(
			names = names,
			fixed = fixed,
			fctName = as.character(match.call()[[1]]),
			...
		))
	}


braincousens_logistic <- function (
		fixed = c(NA, NA, NA, NA, NA),
		names = c("b", "c", "d", "e", "f"),
		method = c("1", "2", "3", "4"),
		ssfct = NULL,
		fctName,
		fctText) {
	numParm <- 5
	if (!is.character(names) | !(length(names) == numParm)) {
		stop("Not correct 'names' argument")
	}
	if (!(length(fixed) == numParm)) {
		stop("Not correct 'fixed' argument")
	}
	notFixed <- is.na(fixed)
	parmVec <- rep(0, numParm)
	parmVec[!notFixed] <- fixed[!notFixed]
	parmVec1 <- parmVec
	parmVec2 <- parmVec
	fct <- function(dose, parm) {
		parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
		parmMat[, notFixed] <- parm
		parmMat[, 2] + (parmMat[, 3] + parmMat[, 5] * dose - parmMat[, 2]) / (1 + exp(parmMat[, 1] * (dose -parmMat[, 4])))
	}
	if (FALSE) {
		ssfct <- function(dataFra) {
			dose2 <- dataFra[, 1]
			resp3 <- dataFra[, 2]
			startVal <- rep(0, numParm)
			startVal[3] <- max(resp3) + 0.001
			startVal[2] <- min(resp3) - 0.001
			startVal[5] <- 0
			if (length(unique(dose2)) == 1) {
				return((c(NA, NA, startVal[3], NA, NA))[notFixed])
			}
			indexT2 <- (dose2 > 0)
			if (!any(indexT2)) {
				return((rep(NA, numParm))[notFixed])
			}
			dose3 <- dose2[indexT2]
			resp3 <- resp3[indexT2]
			logitTrans <- log((startVal[3] - resp3) / (resp3 -
																								 	startVal[2] + 0.001))
			logitFit <- lm(logitTrans ~ log(dose3))
			startVal[4] <- exp((-coef(logitFit)[1] / coef(logitFit)[2]))
			startVal[1] <- coef(logitFit)[2]
			return(startVal[notFixed])
		}
	}
	if (!is.null(ssfct)) {
		ssfct <- ssfct
	}
	else {
		ssfct <- function(dframe) {
			initval <- llogistic()$ssfct(dframe)
			initval[5] <- 0
			return(initval[notFixed])
		}
	}
	names <- names[notFixed]
	deriv1 <- function(dose, parm) {
		parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
		parmMat[, notFixed] <- parm
		t1 <- parmMat[, 3] - parmMat[, 2] + parmMat[, 5] * dose
		t2 <- exp(parmMat[, 1] * (dose - parmMat[, 4]))
		t3 <- 1 + t2
		t4 <- (1 + t2) ^ (-2)
		cbind(-t1 * t2 * t4, 1 - 1 / t3, 1 / t3, t1 * t2 * parmMat[, 1] * t4, dose /
						t3)[, notFixed]
	}
	deriv2 <- NULL
	edfct <- function(
		parm,
		respl,
		reference,
		type,
		lower = 0.001,
		upper = 1000,
		...) {
		interval <- c(lower, upper)
		parmVec[notFixed] <- parm
		p <- EDhelper(parmVec, respl, reference, type)
		tempVal <- (100 - p) / 100
		helpEqn <- function(dose) {
			expVal <- exp(parmVec[1] * (log(dose) - log(parmVec[4])))
			parmVec[5] * (1 + expVal * (1 - parmVec[1])) - 
				(parmVec[3] - parmVec[2]) * expVal * parmVec[1] / dose
		}
		maxAt <- uniroot(helpEqn, interval)$root
		eqn <- function(dose) {
			tempVal * (1 + exp(parmVec[1] * (log(dose) - log(parmVec[4])))) -
				(1 + parmVec[5] * dose / (parmVec[3] - parmVec[2]))
		}
		EDp <- uniroot(eqn, lower = maxAt, upper = upper)$root
		EDdose <- EDp
		tempVal1 <- exp(parmVec[1] * (log(EDdose) - log(parmVec[4])))
		tempVal2 <- parmVec[3] - parmVec[2]
		derParm <-
			c(
				tempVal * tempVal1 * (log(EDdose) - log(parmVec[4])),
				-parmVec[5] * EDdose / ((tempVal2) ^ 2),
				parmVec[5] *
					EDdose / ((tempVal2) ^ 2),
				-tempVal * tempVal1 *
					parmVec[1] / parmVec[4],
				-EDdose / tempVal2
			)
		derDose <- tempVal * tempVal1 * parmVec[1] / EDdose -
			parmVec[5] / tempVal2
		EDder <- derParm / derDose
		return(list(EDp, EDder[notFixed]))
	}
	maxfct <- function(parm,
										 lower = 0.001,
										 upper = 1000) {
		parmVec[notFixed] <- parm
		if (parmVec[1] < 1) {
			stop("Brain-Cousens model with b<1 not meaningful")
		}
		if (parmVec[5] < 0) {
			stop("Brain-Cousens model with f<0 not meaningful")
		}
		optfct <- function(t) {
			expTerm1 <- parmVec[5] * t
			expTerm2 <- exp(parmVec[1] * (log(t) - log(parmVec[4])))
			return(parmVec[5] * (1 + expTerm2) - (parmVec[3] -
																							parmVec[2] + expTerm1) * expTerm2 * parmVec[1] /
						 	t)
		}
		ED1 <- edfct(parm, 1, lower, upper)[[1]]
		doseVec <- exp(seq(log(1e-06), log(ED1), length = 100))
		maxDose <- uniroot(optfct, c((doseVec[optfct(doseVec) >
																						0])[1], ED1))$root
		return(c(maxDose, fct(maxDose, matrix(
			parm, 1, length(names)
		))))
	}
	returnList <- list(
		fct = fct,
		ssfct = ssfct,
		names = names,
		deriv1 = deriv1,
		deriv2 = deriv2,
		edfct = edfct,
		maxfct = maxfct,
		name = ifelse(missing(fctName), as.character(match.call()[[1]]),
									fctName),
		text = ifelse(missing(fctText), "Brain-Cousens (hormesis)",
									fctText),
		noParm = sum(is.na(fixed))
	)
	class(returnList) <- "braincousens_logistic"
	invisible(returnList)
}
