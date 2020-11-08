#This code contains functions for drawing NROY implausibility matrices. The tricky part in using it is having the list of NROY densities and Implausibility matrices in the right format. #Separate files will include code for generating one of these given emulators and data etc. There will also be an example list of lists and separate code for drawing general maps in layout form.
rain <- function(n){
	rainbow(n,st=1-0.01/12,en=2/6,alpha=1)[n:1]
}
#library(rgl)
library(shape)
rain <- function(n){
	intpalette(c("darkgreen","green","yellow","orange","red"),numcol=n)
	}
#dev.new()
#pie(rep(1,6),col=rain(6))
#Breaks for implausibility colours
breakVec <- sort(c(seq(from=20,to=3,len=4),seq(from=3.49999,to=2,len=31),seq(from=1.9999999,to=1,len=46),seq(from=0.9999999,to=0,len=20)))
breakVec <- sort(c(20,seq(from=2.9999,to=2.3,len=8),seq(from=2.29999999,to=1,len=22),seq(from=0.9999999,to=0,len=8)))
breakVec <- sort(c(200,3,2.5,2,1.5,1,0.5,0))
cold1 <- function(n,fra1=0.1,fra2=0.25,fra3=0.1){
#cold1 <- function(n,fra1=0.15,fra2=0.25,fra3=0.2,ga1=1){
    n4 <- round(n*fra3,0)
    n1 <- 1
    n2 <- round(n*fra2,0)
    n3 <- n-n1-n2-n4
    co1 <- hsv(0.68,0.25,0.55)
    co2 <- rainbow(n2,st=3/6,en=3.001/6,alpha=seq(0,0.9,len=n2))
    co3 <- rainbow(n3,st=3/6,en=5.05/6)
    co4 <- rainbow(n4,st=5.1/6,en=5.3/6)
    c(co1,co2,co3,co4)
}
cold1 <- function(n,fra1=0.1,fra2=0.25,fra3=0.3,n1=1){
#cold1 <- function(n,fra1=0.15,fra2=0.25,fra3=0.2,ga1=1){
    n4 <- round(n*fra3,0)
    n2 <- round(n*fra2,0)
    n3 <- n-n1-n2-n4
    co1 <- rep(hsv(0.68,0.25,0.55),n1)
    co2 <- rainbow(n2,st=3/6,en=3.001/6,alpha=seq(0,0.9,len=n2))
    co3 <- rainbow(n3,st=3/6,en=1,alpha=seq(0.7,0.8,len=n3))
    co4 <- rainbow(n4,st=0,en=1/6,alpha=seq(0.8,1,len=n4))
    c(co1,co2,co3,co4)
}

#Implausibility NROY plots for x on [0,1]
	#ImpList is a list of lists of minImplau/optDepth matrices corresponding to the implausibility diagnostics and comes in the following order. 
	#ImpList[[1]] is a list of minImp/optDep for variable 1 and contains each of the interactions with that variable in the order of varnames
	#ImpList[[2]] is a similar list for variable 2 with varnames having lost a member.
imp.layout01 <- function(ImpList,VarNames,VariableDensity=TRUE,newPDF=FALSE,the.title=NULL,newPNG=FALSE,newJPEG=FALSE,newEPS=FALSE,Points=NULL){
	if(newPDF){
		if(is.null(the.title))
		    stop("Require a Title for the PDF")
		else
		    pdf(the.title,compress=FALSE)
	}
	else if(newPNG){
		if(is.null(the.title))
		    stop("Require a Title for the PDF")
		else
		    png(the.title)
	}
	else if(newJPEG){
		if(is.null(the.title))
			stop("Require a JPEG title")
		else
			jpeg(the.title,quality=100)
	}
	else if(newEPS){
		if(is.null(the.title))
			stop("Require an eps title")
		else
			postscript(the.title)
	}
	else{
	#	dev.new()
	}
	tBig <- 0
    n <- length(VarNames)
    n1 <- n-1
    for(j in 1:n1){
		for(i in 1:(max(n1-j+1))){
			tBig <- max(tBig,max(ImpList[[j]][[i]][2,]))
		}
	}
	Mainlevs <- pretty(c(0,tBig),100)
	Maincols <- cold1(length(Mainlevs)-1,n1=1)
    mar <- par("mar")
    tmar <- mar
    tmar[4] <- tmar[2]
    tmar[2] <- 1
    tmar <- tmar/4
    tmar[4] <- 3
    tmar[1] <- 0.3
    tmar[3] <- tmar[1]
    key <- rep(1,n)
    key2 <- rep(2,n)
    numMatrix <- diag(3:(n+2))
    currentNum <- n+3
    for(i in 1:(n-1)){
    	for(j in (i+1):n){
    		numMatrix[i,j] <- currentNum
    		currentNum <- currentNum + 1
    		numMatrix[j,i] <- currentNum
    		currentNum <- currentNum + 1
    	}
    }
    layout(cbind(key2,numMatrix,key),widths=c(lcm(1.6),rep(0.3,n),lcm(1.6)))
    par(mar=tmar,las=1)
    plot.new()
    plot.window(c(0,1),range(Mainlevs),xaxs="i",yaxs="i")
    rect(0,Mainlevs[-length(Mainlevs)],1,Mainlevs[-1],col=Maincols)
    axis(4)
    par(mar=c(0.3,2.2,0.3,0.8))
    plot.new()
    plot.window(c(0,1),c(0,3.45),xaxs="i",yaxs="i")
    bV <- breakVec
    bV[8] <- 3.5
    rect(0,bV[-length(bV)],1,bV[-1],col=rain(7))
    axis(2)
	par(mar=rep(0.1,4),usr=c(0,1,0,1))
	for(k in 1:n){
		plot.new()
		plot.window(c(0,1),c(0,1))
		rect(0,0,1,1)
		text(0.5,0.5,labels=VarNames[k],cex=5/n)
		if(k==1){
			arrows(0,0,0,1.01,col=1,code=2,lwd=0.8,length=0.1)
			arrows(0,0,1.01,0,col=1,code=2,lwd=0.8,length=0.1)
		}
	}
	for(i in 1:(n-1)){
		for(j in (i+1):n){
	        	zmatod <- ImpList[[i]][[j-(i-1)-1]][2,]
	        	tmpn <- sqrt(length(zmatod))
	        	dim(zmatod) <- c(tmpn,tmpn)
	        	x1 <- seq(from=1,by=1,to=tmpn)/tmpn - (1/(2*tmpn))
	        	plot.new()
	        	plot.window(c(0,1),c(0,1))
	        	#.Internal(filledcontour(x1,x1,t(zmatod),Mainlevs,Maincols))
	        	.filled.contour(x1,x1,t(zmatod),Mainlevs,Maincols)
	        	rect(0,0,1,1)
	        	if(!is.null(Points)){
					points(Points[,j],Points[,i],pch=c(16:(16-length(Points[,1]))),col=1,cex=1.2)
					#points(Points[,j],Points[,i],pch=c(16,17),col=1,cex=1.2)
	        	}
	        	if(VariableDensity){
	        		plot.new()
	        		plot.window(c(0,1),c(0,1))
	        		levs <- pretty(range(zmatod),50)
	        		cols <- cold1(length(levs)-1)
	        		#.Internal(filledcontour(x1,x1,zmatod,levs,cols))
	        		.filled.contour(x1,x1,t(zmatod),levs,cols)
	        		rect(0,0,1,1)
	        	}
	        	else{
	        	    zmatimp <- ImpList[[i]][[j-(i-1)-1]][1,]
	        	    tmpn <- sqrt(length(zmatimp))
	        	    dim(zmatimp) <- c(tmpn,tmpn)
	        	    x1 <- seq(from=1,by=1,to=tmpn)/tmpn - (1/(2*tmpn))
	        	    plot.new()
	        	    #.Internal(filledcontour(x1,x1,zmatimp,breakVec,rain(100)))
	        	    .filled.contour(x1,x1,t(zmatimp),breakVec,rain(7))
	        	    rect(0,0,1,1)
	        	    if(!is.null(Points)){
					points(Points[,j],Points[,i],pch=c(16:(16-length(Points[,1]))),col=1,cex=1.2)
					#points(Points[,j],Points[,i],pch=c(16,17),col=1,cex=1.2)
	        		}
	        	}
	        }       
	    }
	}
	
	#Implausibility NROY plots for x on [-1,1]
	
	imp.layoutm11 <- function(ImpList,VarNames,VariableDensity=TRUE,newPDF=FALSE,the.title=NULL,newPNG=FALSE,newJPEG=FALSE,newEPS=FALSE,Points=NULL, tbreakVec=breakVec){
	if(newPDF){
		if(is.null(the.title))
		    stop("Require a Title for the PDF")
		else
		    pdf(the.title,compress=FALSE)
	}
	else if(newPNG){
		if(is.null(the.title))
		    stop("Require a Title for the PDF")
		else
		    png(the.title)
	}
	else if(newJPEG){
		if(is.null(the.title))
			stop("Require a JPEG title")
		else
			jpeg(the.title,quality=100)
	}
	else if(newEPS){
		if(is.null(the.title))
			stop("Require an eps title")
		else
			postscript(the.title)
	}
	else{
		#dev.new()
	}
  par(oma=c(1,0,0,0))
	tBig <- 0
    n <- length(VarNames)
    n1 <- n-1
    for(j in 1:n1){
		for(i in 1:(max(n1-j+1))){
			#tBig <- max(tBig,max(na.omit(ImpList[[j]][[i]][2,])))
			tBig <- max(tBig,max(ImpList[[j]][[i]][2,]))
		}
	}
	Mainlevs <- pretty(c(0,tBig),100)
	Maincols <- cold1(length(Mainlevs)-1,n1=1)
    mar <- par("mar")
    tmar <- mar
    tmar[4] <- tmar[2]
    tmar[2] <- 1
    tmar <- tmar/4
    tmar[4] <- 3
    tmar[1] <- 0.3
    tmar[3] <- tmar[1]
    key <- rep(1,n)
        key2 <- rep(2,n)
    numMatrix <- diag(3:(n+2))
    currentNum <- n+3
    for(i in 1:(n-1)){
    	for(j in (i+1):n){
    		numMatrix[i,j] <- currentNum
    		currentNum <- currentNum + 1
    		numMatrix[j,i] <- currentNum
    		currentNum <- currentNum + 1
    	}
    } 
    layout(cbind(key2,numMatrix,key),widths=c(lcm(1.6),rep(0.3,n),lcm(1.6)))
    par(mar=tmar,las=1)
    plot.new()
    plot.window(c(0,1),range(Mainlevs),xaxs="i",yaxs="i")
    rect(0,Mainlevs[-length(Mainlevs)],1,Mainlevs[-1],col=Maincols)
    axis(4)
    par(mar=c(0.3,2.2,0.3,0.8))
    plot.new()
    upper <- tbreakVec[7]
    gap  <- tbreakVec[3]-tbreakVec[2]
    plot.window(c(0,1),c(0,upper+gap-gap/10),xaxs="i",yaxs="r",yaxp=c(0,upper+gap,8))
    bV <- tbreakVec
    bV[8] <- upper+gap
    rect(0,bV[-length(bV)],1,bV[-1],col=rain(7))
    axis(2)
	par(mar=rep(0.1,4),usr=c(0,1,0,1))
	for(k in 1:n){
		plot.new()
		plot.window(c(-1,1),c(-1,1))
		rect(-1,-1,1,1)
		text(0,0,labels=VarNames[k],cex=5/n)
		if(k==1){
			arrows(-1,-1,-1,1.01,col=1,code=2,lwd=0.8,length=0.1)
			arrows(-1,-1,1.01,-1,col=1,code=2,lwd=0.8,length=0.1)
		}
	}
	for(i in 1:(n-1)){
		for(j in (i+1):n){
	        	zmatod <- ImpList[[i]][[j-(i-1)-1]][2,]
	        	tmpn <- sqrt(length(zmatod))
	        	dim(zmatod) <- c(tmpn,tmpn)
	        	x1 <- 2*(seq(from=1,by=1,to=tmpn)/tmpn - (1/(2*tmpn))) - 1
	        	plot.new()
	        	plot.window(c(-1,1),c(-1,1))
	        	#.Internal(filledcontour(x1,x1,t(zmatod),Mainlevs,Maincols))
	        	.filled.contour(x1,x1,t(zmatod),Mainlevs,Maincols)
	        	rect(-1,-1,1,1)
	        	if(!is.null(Points)){
					points(Points[,j],Points[,i],pch=c(16:(16+length(Points[,1]))-1),col=1,cex=1.2)
	        	}
	        	if(VariableDensity){
	        		plot.new()
	        		plot.window(c(-1,1),c(-1,1))
	        		levs <- pretty(range(zmatod),50)
	        		cols <- cold1(length(levs)-1)
	        		#.Internal(filledcontour(x1,x1,zmatod,levs,cols))
	        		.filled.contour(x1,x1,t(zmatod),levs,cols)
	        		rect(-1,-1,1,1)
	        	    if(!is.null(Points)){
					points(Points[,j], Points[,i], pch=c(16,17), col=1, cex=1.2)
					}
	        	}
	        	else{
	        	    zmatimp <- ImpList[[i]][[j-(i-1)-1]][1,]
	        	    tmpn <- sqrt(length(zmatimp))
	        	    dim(zmatimp) <- c(tmpn,tmpn)
	        	    x1 <- 2*(seq(from=1,by=1,to=tmpn)/tmpn - (1/(2*tmpn))) -1
	        	    plot.new()
	        	    plot.window(c(-1,1),c(-1,1))
	        	    .filled.contour(x1,x1,t(zmatimp),tbreakVec,rain(7))
	        	    rect(-1,-1,1,1)
	        	    if(!is.null(Points)){
					points(Points[,j],Points[,i],pch=c(16:(16-length(Points[,1]))),col=1,cex=1.2)
					#points(Points[,j], Points[,i], pch=c(16,17), col=1, cex=1.2)
					}
	        	}
	        }       
	    }
	}


