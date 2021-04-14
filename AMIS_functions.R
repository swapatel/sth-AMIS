# functions for AMIS algorithm
get_second_weight <- function(sample_prev, data_prev, delta, first_weight){
  w<-sapply(1:length(sample_prev), function(j) length(which((data_prev>sample_prev[j]-delta/2) & (data_prev<=sample_prev[j]+delta/2)))/sum(first_weight[which((sample_prev>sample_prev[j]-delta/2) & (sample_prev<=sample_prev[j]+delta/2))]) ) #f/g from AMIS
  return(w)
}

get_second_weight_v2 <- function(sample_prev, data_prev, delta, first_weight){ #for checking v1 function
  w <- c()
  for(j in 1:length(sample_prev)){
  f_count <- length(which((data_prev>sample_prev[j]-delta/2) & (data_prev<=sample_prev[j]+delta/2)))
  g_count <- sum(first_weight[which((sample_prev>sample_prev[j]-delta/2) & (sample_prev<=sample_prev[j]+delta/2))])
  w <- c(w, f_count/g_count)
  }
  return(w)
}

calculate_ESS <- function(weight_matrix){
  ess<-c()
  n.ius <- ncol(weight_matrix)
  n.sam <- nrow(weight_matrix)
  for(i in 1:n.ius){
    ww <- weight_matrix[,i] #weights for IU
    if( sum(ww)==0){
      www<-0
    } else {
      www<-(sum((ww)^2))^(-1)
    }
    ess<-c(ess, www)
  }
  return(ess)
}

qa_plots <- function(paramw, grp_map_data, plot_folder, scenid){
  n.ius <- dim(grp_map_data)[1]
  n.sam <- dim(grp_map_data)[2] - 1
  IUlist <- map_data_for_grp[,1]
  for(i in 1:n.ius){
    pdf(file=paste0(plot_folder, "/ecdf.", scenid, IUlist[i], ".pdf"))
    par(mfrow=c(1,1))
    data_prev <- grp_map_data[i, 2:(n.sam+1)]
    sample_prev <- paramw[,4]
    weights <- paramw[,5+i]
    ecdf_dist <- TestsCdfs(sample_prev, weights, data_prev, rep(1, n.sampl))
    plot_two_ecdf(sample_prev, data_prev, weights)
    text(.2,.9, paste("sqrd diff=", round(ecdf_dist[4],4)))
    text(.2,.8, paste("max=", round(ecdf_dist[3],4)))
    dev.off() 
  }
}
  