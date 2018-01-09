###############################################
## Plot a 2-Way, 3-Way or 4-Way Venn Diagram ##
###############################################
## Author: Thomas Girke
## Last update: May 4, 2007
## Utility: Plots a non-proportional 2-, 3- or 4-way venn diagram based on provided key vectors
## How to run the script:
##	source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/vennDia.R") 	
## 	x <- sample(letters, 18); y <- sample(letters, 16); z <- sample(letters, 20); w <- sample(letters, 22) # Creates a test sample.
##	qlist <- venndiagram(x=x, y=y, z=z, w=w, unique=T, title="4-Way Venn Diagram", labels=c("Sample1", "Sample2", "Sample3", "Sample4"), plot=T, lines=c(2,3,4,6), lcol=c(2,3,4,6), tcol=1, lwd=3, cex=1, type="4") # Generates 4-way venn diagram 	
##	qlist <- venndiagram(x=x, y=y, z=z, unique=T, title="3-Way Venn Diagram", labels=c("Sample1", "Sample2", "Sample3"), plot=T, lines=c(2,3,4), lcol=c(2,3,4), lwd=3, cex=1.3, type="3") # Generates 3-way venn diagram 	
##	qlist <- venndiagram(x=x, y=y, unique=T, title="2-Way Venn Diagram", labels=c("Sample1", "Sample2"), plot=T, lines=c(2,3), lcol=c(2,3), lwd=3, cex=1.3, type="2") # Generates 2-way venn diagram 	
##	olReport(qlist, missing=".", type=1) # Returns the overlap results in different formats: 'type=1' returns a data frame with the overlap keys and 'type=2' the overlap counts.  
##	venndiagram(x=x, y=y, z=z, w=w, unique=T, labels=c("Sample1", "Sample2", "Sample3", "Sample4"), lines=c(2,3,4,6), lcol=c(2,3,4,6), tcol=1, lwd=3, cex=1, type="4map") # Plots 4-way mapping venn diagram 
##	venndiagram(x=x, y=y, z=z, unique=T, labels=c("Sample1", "Sample2", "Sample3"), lines=c(2,3,4), lcol=c(2,3,4), lwd=3, cex=1.3, type="3map") # Plots 3-way mapping venn diagram 
##	venndiagram(x=x, y=y, unique=T, labels=c("Sample1", "Sample2"), lines=c(2,3), lcol=c(2,3), lwd=3, cex=1.3, type="2map") # Plots 2-way mapping venn diagram 

## Define venndiagram function
venndiagram <- function(x=x, y=y, z=z, w=w, unique=T, title="Venn Diagram", labels=c("x", "y", "z", "w"), lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="3", subtitle=T, ...) {
	## Remove duplicates and NA fields in x, y, z and w
	if(unique==T) {
		x <- unique(x); x <- as.vector(na.omit(x))
		y <- unique(y); y <- as.vector(na.omit(y))
		if(!missing("z")) {
			z <- unique(z); z <- as.vector(na.omit(z))
		}
		if(!missing("w")) {
			w <- unique(w); w <- as.vector(na.omit(w))
		}
	}
	
	## Check valid type selection
	if(!type %in% c("2", "2map", "3", "3map", "4", "4map")) {
		return("Error: the 'type' argument can only be set to one of these values: 2, 2map, 3, 3map, 4, 4map.")	
	}

	## Plot a 2-way venn diagram
	if(type=="2") {
		# Define ovelap queries 
		q1 <- x[x %in% y]
		q2 <- x[!x %in% y]
		q3 <- y[!y %in% x]
		
		## Store query vectors in list
		qlist <- list(q1=q1, q2=q2, q3=q3)
		
		## Perfom query counts
		count <- unlist(lapply(qlist, length))
		countDF <- data.frame(query=names(count) , count=as.vector(count))
		olDF <- data.frame(x=c(5.0, 3.1, 7.0), y=c(6.1, 6.1, 6.1), count=countDF$count)
		mysub <- ifelse(subtitle, paste(paste("N unique: xy =", length(unique(c(x,y)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), sep=""), "")
		
		if(plot==T) {
			## Plot the 2-way venn diagram
			symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...); 
			text(olDF$x, olDF$y, olDF$count, col=tcol, ...); text(c(2.0, 8.0), c(8.8, 8.8), labels[1:2], col=lcol, ...)
		}
		
		## Return query list 
		return(qlist)
	}
	
	## Plot 2-way mapping venn diagram
	if(type=="2map") {
		olDFdebug <- data.frame(x=c(5.0, 3.1, 7.0), y=c(6.1, 6.1, 6.1), count=paste("q", 1:3, sep=""), ...)
		symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...); 
		text(olDFdebug$x, olDFdebug$y, olDFdebug$count, col=tcol, ...); text(c(2.0, 8.0), c(8.8, 8.8), paste(labels[1:2], "=", c("x","y")), col=lcol, ...)
	}
	
	## Plot a 3-way venn diagram
	if(type=="3") {
		## Define ovelap queries 
		q1 <- x[x %in% y & x %in% z]
		q2 <- x[x %in% z]; q2 <- q2[!q2 %in% y]
		q3 <- y[y %in% z]; q3 <- q3[!q3 %in% x]
		q4 <- x[x %in% y]; q4 <- q4[!q4 %in% z]
		q5 <- x[!x %in% y]; q5 <- q5[!q5 %in% z]
		q6 <- y[!y %in% z]; q6 <- q6[!q6 %in% x]
		q7 <- z[!z %in% x]; q7 <- q7[!q7 %in% y]
		
		## Store query vectors in list
		qlist <- list(q1=q1, q2=q2, q3=q3, q4=q4, q5=q5, q6=q6, q7=q7)
		
		## Perfom query counts
		count <- unlist(lapply(qlist, length))
		countDF <- data.frame(query=names(count) , count=as.vector(count))
		olDF <- data.frame(x=c(5.0, 3.8, 6.3, 5.0, 3.0, 7.0, 5.0), y=c(5.6, 4.6, 4.6, 6.9, 6.5, 6.5, 3.0), count=countDF$count)
		mysub <- ifelse(subtitle, paste(paste("N unique: xyz =", length(unique(c(x,y,z)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), paste("; z =", length(unique(z))), sep=""), "")
		
		if(plot==T) {
			## Plot the 3-way venn diagram
			symbols(x=c(4, 6, 5), y = c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...); 
			text(olDF$x, olDF$y, olDF$count, col=tcol, ...); text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), labels[1:3], col=lcol, ...)
		}
			
		## Return query list 
		return(qlist)
	}
	
	## Plot 3-way mapping venn diagram
	if(type=="3map") {
		olDFdebug <- data.frame(x=c(5.0, 3.8, 6.3, 5.0, 3.0, 7.0, 5.0), y=c(5.6, 4.6, 4.6, 6.9, 6.5, 6.5, 3.0), count=paste("q", 1:7, sep=""), ...)
		symbols(x=c(4, 6, 5), y = c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...); 
		text(olDFdebug$x, olDFdebug$y, olDFdebug$count, col=tcol, ...); text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), paste(labels[1:3], "=", c("x","y","z")), col=lcol, ...)
	}
	
	## Plot a 4-way venn diagram
	if(type=="4") {
		## Define ovelap queries 
		xy <- x[x %in% y]; xz <-x[x %in% z]; xw <- x[x %in% w]; yz <- y[y %in% z]; yw <- y[y %in% w]; zw <- z[z %in% w]
		q1 <- xy[xy %in% zw]
		q2 <- xw[xw %in% z]; q2 <- q2[!q2 %in% y]
		q3 <- yz[yz %in% w]; q3 <- q3[!q3 %in% x]
		q4 <- yz[yz %in% x]; q4 <- q4[!q4 %in% w]
		q5 <- xw[xw %in% y]; q5 <- q5[!q5 %in% z]
		q6 <- xy[!xy %in% z]; q6 <- q6[!q6 %in% w]
		q7 <- zw[!zw %in% x]; q7 <- q7[!q7 %in% y]
		q8 <- xz[!xz %in% y]; q8 <- q8[!q8 %in% w] 
		q9 <- yw[!yw %in% x]; q9 <- q9[!q9 %in% z]
		q10 <- x[!x %in% c(y,z,w)]
		q11 <- y[!y %in% c(x,z,w)]
		q12 <- z[!z %in% c(x,y,w)]
		q13 <- w[!w %in% c(x,y,z)]
		q14 <- xw[!xw %in% y]; q14 <- q14[!q14 %in% z]
		q15 <- yz[!yz %in% x]; q15 <- q15[!q15 %in% w]
		
		## Store query vectors in list
		qlist <- list(q1=q1, q2=q2, q3=q3, q4=q4, q5=q5, q6=q6, q7=q7, q8=q8, q9=q9, q10=q10, q11=q11, q12=q12, q13=q13, q14=q14, q15=q15)
			
		## Perfom query counts
		count <- unlist(lapply(qlist, length))
		countDF <- data.frame(query=names(count) , count=as.vector(count))
		olDF <- data.frame(x=c(4.8, 3.9, 5.7, 3.9, 5.7, 4.8, 4.8, 3.0, 6.5, 3.0, 6.5, 3.0, 6.5, 4.8, 4.8), y=c(5.2, 4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 3.2, 3.2, 1.0, 0.4), count=countDF$count)
		mysub <- ifelse(subtitle, paste(paste("N unique: xyzw =", length(unique(c(x,y,z,w)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), paste("; z =", length(unique(z))), paste("; w =", length(unique(w))), sep=""), "")
		
		if(plot==T) {
			## Plot the 4-way venn diagram
			symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...); 
			text(olDF$x[1:13], olDF$y[1:13], olDF$count[1:13], col=tcol, ...) # rows 14-15 of olDF are printed in last step
			text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), labels, col=lcol, ...)
		
			text(c(3.8, 3.8), c(1.0, 0.4), c(paste("Only in ", labels[1], " & ", labels[4], ": ", olDF$count[14], sep=""), paste("Only in ", labels[2], " & ", labels[3], ": ", olDF$count[15], sep="")), col=diacol, ...)
		}
			
		## Return query list 
		return(qlist)
	}
	
	## Plot 4-way mapping venn diagram
	if(type=="4map") {
		olDFdebug <- data.frame(x=c(4.8, 3.9, 5.7, 3.9, 5.7, 4.8, 4.8, 3.0, 6.5, 3.0, 6.5, 3.0, 6.5, 4.8, 4.8), y=c(5.2, 4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 3.2, 3.2, 1.0, 0.4), count=paste("q", 1:15, sep=""), ...)
		symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...); 
		#text(olDFdebug$x[1:13], olDFdebug$y[1:13], olDFdebug$count[1:13], col=tcol, ...); text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), paste(labels, "=", c("x","y","z","w")), col=lcol, ...)
		#text(c(3.8, 3.8), c(0.97, 0.36), c(paste("Only in ", labels[1], " & ", labels[4], ": ", olDFdebug$count[14], sep=""), paste("Only in ", labels[2], " & ", labels[3], ": ", olDFdebug$count[15], sep="")), col=tcol, ...)
	}
}
	
## Generate overlap reports
olReport <- function(qlist=qlist, missing=".", type=1) {
	## Check valid type selection
	if(!type %in% c(1, 2)) {
		return("Error: the 'type' argument can only be set to one of these values: 1 or 2.")	
	}
	
	## Output data frame with overlap keys
	if(type==1) {
		ids <- sort(unique(as.vector(unlist(qlist))))
		qDF <- matrix(ids, nrow=length(ids), ncol=length(qlist), dimnames=list(1:length(ids), names(qlist)))
		lqDF <- as.data.frame(lapply(names(qlist), function(x) qDF[,x] %in% qlist[[x]]))
		colnames(lqDF) <- colnames(qDF)
		lqDF <- as.matrix(lqDF)
		qDF[!lqDF] <- missing
		qDF <- data.frame(IDs=ids, qDF)
		return(qDF)
	}
	
	## Output data frame with overlap counts
	if(type==2) {
		qStat <- data.frame(count=sapply(qlist, length))
		return(qStat)
	}
}

# source("./Scripts/vennDiagram.R")
# jnk<-venndiagram(x=rownames(temp9), y=rownames(b9), type="2", labels=c("old", "new"))
# jnk$q3