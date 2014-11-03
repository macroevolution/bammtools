assignColorBreaks <- function(rates, NCOLORS = 64, spex = "s", log = FALSE, method = c("linear","quantile")) {
	method = match.arg(method, c("linear", "quantile"));
	if (mode(rates) == "numeric") {
		if (log == FALSE) {
			if (method == "quantile")
				bks <- quantile(rates, seq(0,1, length.out=(NCOLORS+1)))
			else
				bks <- seq(min(rates), max(rates), length.out = (NCOLORS+1));
		}
		else {
			if (method == "quantile")
				bks <- quantile(log(rates), seq(0,1, length.out=(NCOLORS+1)))
			else
				bks <- seq(min(log(rates)), max(log(rates)), length.out = (NCOLORS+1));
		}	
	}
	else if (mode(rates) == "list") {
		if (tolower(spex) == "s") {
			if (log == FALSE) {
				if (method == "quantile")
					bks <- quantile(rates[[1]], seq(0,1, length.out=(NCOLORS+1)))
				else
					bks <- seq(min(rates[[1]]), max(rates[[1]]), length.out = (NCOLORS+1))
			}
			else {
				if (method == "quantile")	
					bks <- quantile(log(rates[[1]]), seq(0,1, length.out=(NCOLORS+1)))
				else
					bks <- seq(min(log(rates[[1]])), max(log(rates[[2]])), length.out=(NCOLORS+1));
			}
		}
		else if (tolower(spex) == "e") {
			if (log == FALSE) {
				if (method == "quantile")
					bks <- quantile(rates[[2]], seq(0,1, length.out=(NCOLORS+1)))
				else
					bks <- seq(min(rates[[2]]), max(rates[[2]]), length.out = (NCOLORS+1));
			}
			else {
				if (method == "quantile")
					bks <- quantile(log(rates[[2]]), seq(0,1, length.out=(NCOLORS+1)))
				else
					bks <- seq(min(log(rates[[2]])), max(log(rates[[2]])), length.out=(NCOLORS+1));
			}
		}
		else {
			if (log == FALSE) {
				if (method == "quantile")	
					bks <- quantile(rates[[1]] - rates[[2]], seq(0,1, length.out=(NCOLORS+1)))
				else
					bks <- seq(min(rates[[1]] - rates[[2]]), max(rates[[1]] - rates[[2]]), length.out=(NCOLORS+1));
			}
			else { 
				if (method == "quantile")
					bks <- quantile(log(rates[[1]] - rates[[2]]), seq(0,1, length.out=(NCOLORS+1)))
				else
					bks <- seq(min(log(rates[[1]] - rates[[2]])), max(min(log(rates[[1]] - rates[[2]]))), length.out=(NCOLORS+1) );
			}
		}
	}
	return(bks);
}
