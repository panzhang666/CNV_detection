# title: "myraid_CNV"
# author: "PanZhang"
# date: "5/26/2020"
# output: word_document

  

#The statistical model for the likelihood of observing a given number of mapped reads di,j 
#at a given genomic position i for sample j with copy number ci,j is constructed following 
#the paper "Development and validation of a 36-gene sequencing assay for hereditary cancer risk assessment".

cal_p_var <- function(dij, cij, ui, uj, ri){
  # Function to calculate the probility of dij when copy number = cij based on the negative bionominal model - p(dij |cij ) = NegBinom(dij |μ = cij*μi*μj ,r = ri).
  # dij - depth
  # cij - copy number
  # μi- the average depth for that targeted location across samples 
  # μj - the average depth for that particular sample across targeted positions
  # ri - variance parameter  
  #output: probility of dij when copy number = cij
  
  mu0 <-  cij*ui*uj
  
  densi <- dnbinom(dij, size = 1/ri, mu=mu0)
  
  return(densi)
}


findruns <- function(x, k, cnv, ind, probs, ethni){
  # Function to find the a deletion (cnv =1) or duplication(cnv =3) that is any contiguous stretch of at least k well behaved probes.
  
  n <- length(x)
  runs <- NULL
  last_i <- 0
  
  for(i in 1:(n-k+1)){
    
    if(all( x[i:(i+k-1)] == cnv)){
      
      if ((i -last_i) > 1){
        
        runs <- rbind( runs , c(ind, ethni[ind], i, probs[i], 4))
        last_i <- i
        
      }
      
      else{
        
        d<- dim(runs)
        runs[d[1],d[2]] <- as.numeric(runs[d[1],d[2]]) + 1
        last_i <- i
        
      }
    }
  }
  
  return(runs)
  
}


breakpoint <- function(cnv,probs ){
  # Function to summarize breakpoint positions on a per-ethnicity basis.
  
  end_prob <- apply(cnv, 1, function(x) probs[as.numeric(x[3]) + as.numeric(x[5])-1])
  cnv <- cbind(cnv,end_prob)
  
  ethnicity <- as.character(cnv$ethnicity)
  breakpoint5 <- as.character(cnv$start_prob)
  breatpoint3 <- as.character(cnv$end_prob)
  
  new_cnv <- as.data.frame(cbind(ethnicity, breakpoint5, breatpoint3 ))
  
  out <- ddply(new_cnv,.(ethnicity,breakpoint5, breatpoint3 ),nrow)
  colnames(out) <- c("ethnicity", "5' breakpoint", "3' breakpoint", "count")
  
  return (list(out,cnv))
}



cal_frequency <- function (result,ethni){
  #Function to characterize the deletion and duplication frequencies on a per-ethnicity basis.
  
  N_A <- as.numeric(table(ethni)["A"])
  N_B <- as.numeric(table(ethni)["B"])
  N_C <- as.numeric(table(ethni)["C"])
  
  extreact <- cbind(as.character(result$index), as.character(result$ethnicity))
  extreact_dup <- as.data.frame(extreact[!duplicated(extreact),])
  
  cn_A <- as.numeric(table(extreact_dup$V2)["A"])
  cn_B <- as.numeric(table(extreact_dup$V2)["B"])
  cn_C <- as.numeric(table(extreact_dup$V2)["C"])
  
  return(list(round(cn_A/N_A*100, digits = 2), round(cn_B/N_B*100,digits = 2), round(cn_C/N_C*100, digits = 2)))
} 


summ_breakpoint <- function(breakpoint_info){
  # Function to calculation the frequency of breakpoint pairs on a per-ethnicity basis 
  new_bk <- NULL
  for (e in c("A", "B", "C")){
    bk <- breakpoint_info[which(breakpoint_info$ethnicity == e),]
    freq <- round(bk$count/sum(bk$count),digits = 3)
    bk <- cbind(bk,freq)
    bk <- bk[order(bk$freq, decreasing = T),] 
    new_bk <- rbind(new_bk, bk)
  }
  colnames(new_bk) <- c("ethnicity", "5' breakpoint", "3' breakpoint", "count", "Frequency")
  return (new_bk)
}

# read cnsl_data
setwd("/Users/panzhang/Desktop/PHD/myraid")
raw_data <- read.csv("cnsl_data.csv")

value_data <- raw_data[,3:102]
# check prob quality and remove the prob with low quality
prob_abund <- apply(value_data, 2, function(x) length(which((x < 10))))

# CNSL_probe_5 is removed since its depth is less than 10 among 503 samples ( >10 )
which (prob_abund > 10)
value_data <- value_data[, prob_abund <10]

table(value_data == NA)

#calculate parameters for negative bionominal distribution
dimen = dim(value_data)
dimen
mu_i = apply(value_data,2,mean)
f_muj <- apply(value_data,1,function(x) x/mu_i)
mu_j <- apply(f_muj, 2, mean)
r_i = (mu_i^2)/(apply(value_data,2,var) - mu_i)


# Predict copy number for each dij 
# Step 1 - calculate the density of dij in NegBinom when cn = 1,2,3, respectively
# Step 2 - Compare the 3 densities and assign the cn with largest density to dij
cnv_data2 <- matrix(NA,nrow=dimen[1], ncol=dimen[2])
for (i in seq(1, length(mu_i))){
  
  for (j in seq(1, length(mu_j))){
    # Since the cn of healthy samples and nonCNSL region is exprected to be 2, the cij = 1/2 means cn =1, cij means cn = 2, cij = 3/2 means cn =3. 
    #Here, cij is the input of function cal_p_var  
    
    dens1 <- cal_p_var(round(value_data[j,i]), 1/2, mu_i[i], mu_j[j], r_i[i])
    dens2 <- cal_p_var(round(value_data[j,i]), 1, mu_i[i], mu_j[j], r_i[i])
    dens3 <- cal_p_var(round(value_data[j,i]), 3/2, mu_i[i], mu_j[j], r_i[i]) 
    
    p <- c(dens1, dens2,dens3)
    x <- which(p==max(p),arr.ind=TRUE)
    if (length(x) == 1){
      cnv_data2[j,i] <- x
    }
  }
}
colnames(cnv_data2) <- colnames(value_data)

probs <- colnames(value_data)
ethni <- as.vector(raw_data$ethnicity)
cnv_data2[is.na(cnv_data2)] <- 2


info_gain <- NULL
info_loss <- NULL
for (h in seq(1,10000)){
  pos_gain <- findruns(cnv_data2[h,], 4, 3, h, probs, ethni)
  pos_loss <- findruns(cnv_data2[h,], 4, 1, h, probs, ethni)
  if (length(pos_gain)){
    info_gain <- rbind(info_gain,pos_gain )
  }
  if (length(pos_loss)){
    info_loss <- rbind(info_loss,pos_loss )
  }
}
gain <- as.data.frame(info_gain)
loss <- as.data.frame(info_loss)

colnames(gain) <- c("index", "ethnicity", "position", "start_prob", "N_prob")
colnames(loss) <- c("index", "ethnicity", "position", "start_prob", "N_prob")

#filter the deletion or duplication at non-CNSL region
gain <- gain[apply(gain, 1, function(x) x[3] <=50 ),]
loss <- loss[apply(loss, 1, function(x) x[3] <=50 ),]


N_A <- as.numeric(table(ethni)["A"])
N_B <- as.numeric(table(ethni)["B"])
N_C <- as.numeric(table(ethni)["C"])

duplication_freq <- cal_frequency(gain, ethni)

deletion_freq <- cal_frequency(loss, ethni)

#deletion and duplication frequencies on a per-ethnicity basis 
freq_table <- cbind(duplication_freq, deletion_freq)
rownames(freq_table) <- c("A", "B", "C")
colnames(freq_table) <- c("duplication (%)","deletion (%)" )
freq_table



#breakpoint positions on a per-ethnicity basis
gain_table <- breakpoint(gain, probs)[2]


gain_breakpoint <- as.data.frame(breakpoint(gain, probs)[1])
dup_ethn_breakpoint_freq <- summ_breakpoint(gain_breakpoint)

loss_table <- as.data.frame(breakpoint(loss, probs)[2])
loss_breakpoint <- as.data.frame(breakpoint(loss, probs)[1])

del_ethn_breakpoint_freq <- summ_breakpoint(loss_breakpoint)

write.csv(del_ethn_breakpoint_freq, "del_ethn_breakpoint_freq.csv",quote = FALSE, row.names = FALSE)
write.csv(dup_ethn_breakpoint_freq, "dup_ethn_breakpoint_freq.csv",quote = FALSE, row.names = FALSE)

