cov4gappy <- function(F1, F2=NULL){
    if(is.null(F2)){
        F1 <- as.matrix(F1)
        F1_val<-replace(F1, which(!is.na(F1)), 1)
        F1_val<-replace(F1_val, which(is.na(F1_val)), 0) 
        n_pairs=(t(F1_val)%*%F1_val)
 
        F1<-replace(F1, which(is.na(F1)), 0)
        cov_mat <- (t(F1)%*%F1)/n_pairs
        cov_mat <- replace(cov_mat, which(is.na(cov_mat)), 0)
    }
 
    if(!is.null(F2)){
        if(dim(F1)[1] == dim(F2)[1]){
            F1 <- as.matrix(F1)
            F2 <- as.matrix(F2)
 
            F1_val<-replace(F1, which(!is.na(F1)), 1)
            F1_val<-replace(F1_val, which(is.na(F1_val)), 0) 
            F2_val<-replace(F2, which(!is.na(F2)), 1)
            F2_val<-replace(F2_val, which(is.na(F2_val)), 0) 
            n_pairs=(t(F1_val)%*%F2_val)
 
            F1<-replace(F1, which(is.na(F1)), 0)
            F2<-replace(F2, which(is.na(F2)), 0)
            cov_mat <- (t(F1)%*%F2)/n_pairs
            cov_mat <- replace(cov_mat, which(is.na(cov_mat)), 0)
 
        } else {
            print("ERROR; matrices columns not of the same lengths")
        }
    }
 
    cov_mat
}
