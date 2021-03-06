# functions for plotting the empirical CDF 

Ecdf <- function(Samples, Weights=NULL) # calculate empirical cdf
{
  if (is.null(Weights)) {Weights=rep(1/length(Samples),length(Samples))}
  Matrix = data.frame(Samples, Weights)
  Matrix = Matrix[order(Matrix$Samples),]
  Matrix$Heights=cumsum(Matrix$Weights)
  return(Matrix)
}

plot_two_ecdf <- function(sample_prev, data_prev, weights){
  par(mfrow=c(1,1))
  weights = weights/sum(weights)
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="x", ylab = "Empirical cumulative distribution")
  sample_ecdf_matrix <- Ecdf(sample_prev, Weights=weights)
  data_ecdf_matrix <- Ecdf(data_prev)
  lines(stepfun(sample_ecdf_matrix$Samples, c(0,sample_ecdf_matrix$Height)), col=4) # red 
  lines(stepfun(data_ecdf_matrix$Samples, c(0,data_ecdf_matrix$Height)), col=2) # blue
  abline(h=0, lty="dashed", col="gray")
  abline(h=1, lty="dashed", col="gray")
}

add.points <- function(M,extra) {
  L<-length(extra)
  O<-data.frame(Samples=extra,Weights=rep(0,L),Heights=rep(NA,L))
  rownames(O)<-paste("e",1:length(extra))
  M<-rbind(M,O)
  M<-M[order(M$Samples),]
  for (i in 1:length(M$Samples)) {
    if (is.na(M$Heights[i])) {
      if (i==1) {
        M$Heights[i]<-0
      } else {
        M$Heights[i]<-M$Heights[i-1]
      }
    }
  }
  return(M)   
} 

TestsCdfs <- function(s1,w1,s2,w2) {
	# add normalization of weights
	w1 <- w1/sum(w1)
	w2 <- w2/sum(w2)
  Ecdf1<-add.points(Ecdf(s1,w1),s2)
  Ecdf2<-add.points(Ecdf(s2,w2),s1)
  dheight<-Ecdf1$Height-Ecdf2$Height
  dwidth<-diff(Ecdf1$Samples)
  #wh.pos<-which(dheight>0)
  #wh.neg<-which(dheight<0)
  wh.pos<-which(dheight[-length(dheight)]>0)
  wh.neg<-which(dheight[-length(dheight)]<0)
  #if (Ecdf1$Height[length(Ecdf1$Height)]!=1) {stop("Error: cdf1 does not reach 1!")} 
  # It is one when you print the value - I guess is rounding error? 
  #if (Ecdf2$Height[length(Ecdf2$Height)]!=1) {stop("Error: cdf2 does not reach 1!")}
  A.pos<-sum(dheight[wh.pos]*dwidth[wh.pos])
  A.neg<-sum(dheight[wh.neg]*dwidth[wh.neg])
  Distance<-sum((dheight[wh.pos]^2)*dwidth[wh.pos]) + sum((dheight[wh.neg]^2)*dwidth[wh.neg])
  return(c(A.pos,abs(A.neg),max(abs(dheight), na.rm=T), Distance))
}  
