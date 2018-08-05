#####Code to generate unique IDs for max sample you want. Later you can sample from this to create varying data sets##########
max_sample<-1000000
letter_length=4
lmatrix<-matrix(LETTERS[floor(runif(max_sample*letter_length,1,26))],ncol=letter_length)
lpart<-t(data.table(lapply(1:max_sample,function(i) paste0(data.table(lmatrix)[i,],collapse = ""))))
test_fn<-function(max_sample,lmatrix){
  X<-lapply(1:max_sample,function(i) paste0(data.table(lmatrix)[i,],collapse = ""))
  return(X)
}

npart<-floor(runif(max_sample,1000000,9999999))
ID_all<-paste0(lpart,npart)
Basic_data<-setnames(data.table(ID_all)[,Group:=rep(c("Comparison","Treatment"),c(max_sample/2,max_sample/2))],"ID_all","ID") #Initial pool from which to sample later can change the ratio of control to treatment.