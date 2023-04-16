#' There are cluster size limits and kmean that maintains the distribution of the original data
#'
#' @param dataf input data,a dataframe
#' @param K the number of the clusters
#' @param iter the number of iter
#' @param min_num Each cluster contains the fewest number of points
#' @param max_num The maximum number of points each cluster contains
#' @param b A multiple of the mean number of candidate points for each cluster center, which must be greater than 1, with a default value of 3
#' @param c Ensure that each point is assigned to at least c center points
#' @param norm Whether the metadata is positively distributed   bool
#' @param d The composition ratio of the loss function, the Euclidean distance of 1-d, and the probability error of d
#'
#' @return center and mix_matrix
#' @export igraph,mvtnorm,tidyverse,lpSolve,dplyr
#'
#' @examples
#' load("TLkmean/data/data.RData")
#' out=constract_k_mean(dataf=dataf,K=500,iter=5,min_num=8,max_num=15,b=2,c=3,TRUE,0.6)
#'
library(igraph)
library("mvtnorm")
library(tidyverse)
library(lpSolve)
library(dplyr)
constract_k_mean<-function(dataf,K,iter,min_num,max_num,b=2,c=3,norm=TRUE,d=0.5){
  density_fun<-function(dataf,x){
    out_=data.frame()
    for (i in 1:ncol(dataf)){
      temp_i=density(dataf[,i])
      out_<-rbind(out_,spline(temp_i$x,temp_i$y,xout=x[,i])$y)
    }
    return(t(out_))
  }
  cost_fun<-function(scatter_,center,b=2,c=4,norm=TRUE,d=0.5){
    sumx_sca=rowSums(scatter_^2)
    sumx_cent=rowSums(center^2)
    print("scatter_center")
    print(dim(scatter_))
    print(dim(center))
    scatter_center<-scatter_ %*% t(center)
    out1=sqrt(t(t(-2*scatter_center+sumx_sca+0.00001)+sumx_cent))
    if (norm){
      pr_sca=dmvnorm(scatter_,apply(scatter_,2,mean),cov(scatter_))
      pr_cent=dmvnorm(center,apply(center,2,mean),cov(center))
      out2=abs(log(pr_sca%*%t(1/pr_cent)))
    }
    else{
      pr_sca<-density_fun(scatter_,scatter_)
      pr_cent<-density_fun(scatter_,center)
      out2=abs(log(pr_sca%*%t(1/pr_cent)))/ncol(center)
    }
    # c(matrix(out1+out2,nrow = nrow(scatter_)*nrow(center),byrow = TRUE))
    out=(1-d)*out1+d*out2
    labe=min(ceiling(nrow(scatter_)/nrow(center))*c,nrow(scatter_))
    cost<-c()
    from<-c()
    to<-c()
    for (i in 1:nrow(center)){
      min_loc=order(out[,i],decreasing=FALSE)[1:labe]
      cost<-c(cost,out[,i][min_loc])
      to<-c(to,rep(i+nrow(scatter_),labe))
      from<-c(from,min_loc)
    }
    for (j in 1:nrow(scatter_)){
      min_loc=order(out[j,],decreasing=FALSE)[1:b]
      cost<-c(cost,out[j,][min_loc])
      to<-c(to,min_loc+nrow(scatter_))
      from<-c(from,rep(j,b))
    }
    out_=data.frame(cost,from,to)
    out_<-out_ %>% distinct(from,to, .keep_all = T)
    return(out_)
  }
  to_bulid_graph<-function(dataf,center,b=2,c=4,norm,d){
    edgelist<-cost_fun(dataf,center,b=b,c=c,norm,d)
    #  print("cost_fun--length(unique(edgelist$from))")
    #  print(length(unique(edgelist$from)))
    g <- graph_from_edgelist(as.matrix(edgelist[,c('from','to')]))
    E(g)$cost <- edgelist$cost
    edgelist$ID <- seq(1, nrow(edgelist))
    edgelist

  }
  createConstraintsMatrix <- function(edges, min_samples,max_samples) {

    names_edges <- edges$ID
    numberof_edges <- length(names_edges)

    names_nodes <- c(edges$from, edges$to) %>% unique

    numberof_nodes <- length(names_nodes)

    constraints <- list(
      lhs = NA,
      dir = NA,
      rhs = NA)


    nodeflow <- matrix(0,
                       nrow = numberof_nodes,
                       ncol = numberof_edges,
                       dimnames = list(names_nodes, names_edges))

    #t1=Sys.time()
    sourcenode_id <- edges$from %>%unique
    targetnode_id <- edges$to %>%unique

    for (i in sourcenode_id){
      nodeflow[rownames(nodeflow)==i,][edges$ID[edges$from==i]]<-1
    }
    for (i in targetnode_id){
      nodeflow[rownames(nodeflow)==i,][edges$ID[edges$to==i]]<-1
    }


    nodeflow_source <- nodeflow[rownames(nodeflow) %in% sourcenode_id,]
    print("dim(nodeflow_source)")
    print(dim(nodeflow_source))
    nodeflow_target <- nodeflow[rownames(nodeflow) %in%targetnode_id,]
    print("dim(nodeflow_target")
    print(dim(nodeflow_target))

    constraints$lhs <- nodeflow_source
    constraints$dir <- rep('==', times = nrow(nodeflow_source))
    constraints$rhs <- rep(1, times = nrow(nodeflow_source))


    constraints$lhs <- rbind(constraints$lhs,nodeflow_target)
    constraints$dir <- c(constraints$dir, rep('>=', times = nrow(nodeflow_target)))
    constraints$rhs <- c(constraints$rhs, rep(min_samples,times = nrow(nodeflow_target)))

    constraints$lhs <- rbind(constraints$lhs,nodeflow_target)
    constraints$dir <- c(constraints$dir, rep('<=', times = nrow(nodeflow_target)))
    constraints$rhs <- c(constraints$rhs, rep(max_samples,times = nrow(nodeflow_target)))

    constraints$lhs <- rbind(constraints$lhs,rep(1,numberof_edges))
    constraints$dir <- c(constraints$dir, "==")
    constraints$rhs <- c(constraints$rhs, length(sourcenode_id))





    return(constraints)
  }
  new_center<-function(edgelist,dataf){

    edgelist=edgelist[edgelist$flow==1,]
    print("dim(edgelist)")
    print(dim(edgelist))
    print("length(unique(edgelist$to))")
    print(length(unique(edgelist$to)))
    newcenter<-data.frame()
    for (i in unique(edgelist$to)){
      newcenter<-rbind(newcenter,apply(as.data.frame(dataf[edgelist[edgelist$to==i,]$from,]),2,mean))
    }
    newcenter
  }
  mix_matrix<-function(edgelist,center,dataf){


    k=dim(center)[1]
    p=dim(dataf)[1]
    out_matrix=matrix(0,nrow=k,ncol=p)
    usefull_edgelist=edgelist[edgelist$flow==1,]

    for (i in 1:nrow(usefull_edgelist)){
      num=length(usefull_edgelist$to[usefull_edgelist$to==usefull_edgelist$to[i]])
      print(usefull_edgelist$to[i]-p)
      out_matrix[usefull_edgelist$to[i]-p,usefull_edgelist$from[i]]=1/num
    }

    out<-list(center=center,
              matrix=out_matrix)
    return(out)

  }

  out<-kmeans(dataf,K,10)
  center<-as.matrix(out$centers)
  for (i in 1:iter){
    edgelist<-to_bulid_graph(dataf,center,b=b,c=c,norm,d)
    #  print("length(unique(edgelist$from))")
    print(length(unique(edgelist$from)))
    constraintsMatrix<-createConstraintsMatrix(edgelist,min_num,max_num)
    solution <- lp(
      direction = 'min',
      objective.in = edgelist$cost,
      const.mat = constraintsMatrix$lhs,
      const.dir = constraintsMatrix$dir,
      const.rhs = constraintsMatrix$rhs)
    edgelist$flow<-round(solution$solution,0)
    C=c
    while(sum(edgelist$flow==1)<K){
      C=C+1
      edgelist<-to_bulid_graph(dataf,center,b,C,norm,d)
      #  print("length(unique(edgelist$from))")
      print(length(unique(edgelist$from)))
      constraintsMatrix<-createConstraintsMatrix(edgelist,min_num,max_num)
      solution <- lp(
        direction = 'min',
        objective.in = edgelist$cost,
        const.mat = constraintsMatrix$lhs,
        const.dir = constraintsMatrix$dir,
        const.rhs = constraintsMatrix$rhs)
      edgelist$flow<-round(solution$solution,0)
    }
    center<-new_center(edgelist,dataf)
    #print(table(edgelist[edgelist$flow==1,]$to))
    #print(length(unique(edgelist[edgelist$flow==1,]$from)))

  }
  out_center_matrix=mix_matrix(edgelist,center,dataf)
  return(out_center_matrix)
}

#dataf<-read.csv("D:/Appdnowload/016261_files1/016261_files2.txt",sep = "\t")
#dataf<-as.matrix(dataf[1:5000,2:4])
#save(dataf,file="D:/rwork/TLkmean/data/data.RData")



#library(igraph)
#library("mvtnorm")
#library(tidyverse)
#library(lpSolve)
#library(dplyr)


#dataf:表型数据
#k 聚类个数
#iter 中心点迭代次数
#min_num 每个聚类最少元素个数
#max_num 每个聚类最多元素个数
#out=constract_k_mean(dataf=dataf,K=100,iter=20,min_num=8,max_num=15,b=2,c=4,TRUE,0.6)

#中心点
#out$center


#mix_matrix  在右边乘上SNP矩阵 就是 混样的SNP型
#out$matrix








