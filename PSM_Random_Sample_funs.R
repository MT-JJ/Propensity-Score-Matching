####Type of variable parameters not dependent on sample size.  Even though the distribution variables have sample size (nsubjects) in them this will be automatically determined from the sample size specified in function below. They will only evaluate when you run the sample function and will use whatever sample size you specified. Created for education policy analysts but could easily be adapted for any sampling purpose.

#Stores future commands for use in generating covariates for models based on various statistical distributions.
rand_normal_GPA<-substitute(abs(rnorm(nsubjects,mean=GPA_mean)))
prob_binequal<-substitute(as.numeric(rbinom(nsubjects,1,0.5)))
prob_binbiasup<-substitute(as.numeric(rbinom(nsubjects,1,0.8)))
prob_binbiasdown<-substitute(as.numeric(rbinom(nsubjects,1,0.2)))
prob_int<-substitute(as.numeric(rpois(nsubjects,19)))
unif_age_hs<-substitute(as.numeric(ceiling(runif(nsubjects,min=13,max=18))))
unif_age_ms<-substitute(as.numeric(ceiling(runif(nsubjects,min=9,max = 13))))
unif_age_elem<-substitute(as.numeric(ceiling(runif(nsubjects,min=4,max = 9))))
SPED_dist<-substitute(sample(x=names(SPED_weights),size = nsubjects,replace = T,prob = SPED_weights))
Gender_dist<-substitute(sample(x=names(Gender_levels),size=nsubjects,replace = T,prob=Gender_levels))
Race_dist<-substitute(sample(x=Race_cats,size=nsubjects,replace = T,prob = Race_probs))
FRL_dist<-substitute(sample(x=names(FRL_cat),size = nsubjects,replace=T,prob=FRL_cat))
MCA_math_dist<-substitute(round(rnorm(nsubjects,mean=MCA_math,sd=MCA_math_sd)))
MCA_read_dist<-substitute(round(rnorm(nsubjects,mean=MCA_read,sd=MCA_read_sd)))

######Individual parameters that are called by the above expressions when they are evaluated. Currently default to demographics that reflect the state of Minnesota but can be modified as desired.

Race_cats<-c("American Indian","Asian","Black","Hispanic","White","Other NonWhite")
Race_probs<-setNames(c(0.016,0.067,0.107,0.09,0.675,0.045),Race_cats)
SPED_weights<-setNames(c(0.15,0.85),c("Y","N"))
Gender_levels<-setNames(c(0.46,0.51,0.03),c("Female","Male","Neither"))
FRL_cat<-setNames(c(0.38,0.62),c("Y","N"))
GPA_max=4
GPA_round=2
GPA_mean=2.5
NAEP8_mean=294
MCA_math<-1150.37
MCA_math_sd<-11.83
MCA_math_range<-c(1102,1195)
MCA_read<-1053.23
MCA_read_sd<-10.83
MCA_read_range<-c(1013,1094)

#Note: Var_type_str must be present in the global environment before function will work.
var_type_str<-c(Gender=Gender_dist,Age=unif_age_hs,FRL=FRL_dist,Race_Ethnicity=Race_dist,Cum_GPA=rand_normal_GPA,ACT_COMP=prob_int,SPED=SPED_dist)

##This is the function and can be specified by ratio which defaults to 3:1 comparison to treatment and you must specify total subjects or could alternatively use the numtreat and numcontrol parameters if you prefer to enter set numbers more easily. These default to NULL.

Sample_fn<-function(X,control_ratio=3,total_sub,var_type_str=eval(as.name("var_type_str"),envir = parent.frame()),numtreat=NULL,numcontrol=NULL){
  require(data.table)
  X<-data.table(X)
  if(!is.null(control_ratio)){
    if(!is.null(total_sub)){
      treat_tot<-floor(total_sub/(control_ratio+1))
      Treat<-copy(X)[Group=="Treatment"][ID%in%sample(ID,treat_tot)]
      Control<-copy(X)[Group=="Comparison"][ID%in%sample(ID,total_sub-treat_tot)]
      X<-rbind(Treat,Control)
    } else stop("You need to enter the total sample size from which to apply ratio or specify 'numtreat' and 'numcontrol' before sampling can take place")
  } else{
    if(!is.null(numtreat)&!is.null(numcontrol)){
      Treat<-copy(X)[Group=="Treatment"][ID%in%sample(ID,numtreat)]
      Control<-copy(X)[Group:="Comparison"][ID%in%sample(ID,numcontrol)]
      X<-rbind(Treat,Control)
    } else{stop("You need to enter values for both 'numtreat' and 'numcontrol' or specify a ratio and 'total_sub' (total number of subjects).")}
  }
  nsubjects <<- nrow(X)
  
  MatchingVars<-data.frame(lapply(var_type_str,eval))
  categorical_vars<-names(which(sapply(MatchingVars,is.factor)))
  cat_ref<-paste0(categorical_vars,unlist(lapply(categorical_vars,function(i) levels(copy(MatchingVars)[[i]])[1])))
  
  categorical_list<-lapply(categorical_vars,function(i) data.table(data.frame(model.matrix(as.formula(paste("~",i)),MatchingVars))))
  
  categorical_expand<-Reduce(cbind,lapply(seq_along(categorical_list),function(i) categorical_list[[i]][,cat_ref[i]:=1*(rowSums(copy(.SD))==1)][,-"X.Intercept.",with=F]))
  final_results<-data.table(cbind(copy(X),MatchingVars,categorical_expand))[,Cum_GPA:=round(ifelse(Cum_GPA>4,GPA_max,Cum_GPA),digits = GPA_round)][,ACT_COMP:=ifelse(ACT_COMP>36,36,ACT_COMP)][,-categorical_vars,with=F]
  return(setnames(final_results,make.names(colnames(final_results))))
}


#####This is used to have finer control over initial balance of covariates between treatment and comparison group. Currently only binomial, normal, and uniform distributions are supported. 'change_vars_dist_list' must be entered as a named list of those distributions with names corresponding to column names in the sample set. Ex: list(Gender="binomial"). Supports multiple changes at once as long as appropriate parameters are aligned correctly.

Update_sample_fn<-function(sample_set,grouplevel=c("Both","Treat","Comparison"),change_vars_dist_list,IDvars="ID",Groupvars="Group",new_prob=NULL, new_mean=NULL,new_sd=NULL,new_range=NULL){
  if(!is.list(change_vars_dist_list)|is.null(change_vars_dist_list)){
    stop("You must enter a named list with the name of the variable(s) you wish to change and its value a distribution of type 'binomial','normal',or uniform")
  }
  if(sum(grepl("^N",toupper(unlist(change_vars_dist_list))))>0){
    print("If using the normal distribution for GPA the upper values may exceed the maximum GPA, so adjust accordingly.")
  }
  for(i in seq_along(change_vars_dist_list)){
    switch(substr(toupper(change_vars_dist_list[[i]]),1,1),
           "B"={
             Num_cats<-sample_set[,length(grep(names(change_vars_dist_list)[i],colnames(copy(.SD))))]
             
             All_cat<-rev_indicator_converter(copy(sample_set)[,c(IDvars,Groupvars,grep(names(change_vars_dist_list)[i],colnames(sample_set),value = T)),with=F],names(change_vars_dist_list)[i],IDvars=c(IDvars,Groupvars))[copy(sample_set)[,c(IDvars,Groupvars),with=F],on=c(IDvars,Groupvars)]
             
             if(All_cat[,length(unique(eval(as.name(names(change_vars_dist_list)[i]))))]!=Num_cats){
               stop(paste0("Error in ",sQuote(names(change_vars_dist_list)[i]),". Number of indicator categories is not equal to number of levels in factor you wish to change. Check spelling or set up indicators with identical prefix before proceeding."))
             } 
             if(length(new_prob[[i]])!=Num_cats){
               if(length(new_prob[[i]])==Num_cats-1){
                 new_prob[[i]]<-c(new_prob[[i]],1-sum(new_prob[[i]]))
                } else {stop(paste0("Number of probabilities entered for ",sQuote(names(change_vars_dist_list)[i])," not equal to number of unique categories."))}
             }
             
             if(sum(new_prob[[i]])!=1){stop(paste0("Sum of the probabilities entered for ",sQuote(names(change_vars_dist_list)[i])," is not equal to 1."))}
             
            switch(substr(toupper(grouplevel[1]),1,2),
                      "TR"={
                        
                        new_sample_set<-indicator_converter(All_cat[eval(as.name(Groupvars))=="Treatment",names(change_vars_dist_list)[i]:=sample(x=sort(unique(eval(as.name(names(change_vars_dist_list)[i])))),size=.N,replace=T,prob=new_prob[[i]])],c(IDvars,Groupvars))[sample_set,on=c(IDvars,Groupvars)]
                        sample_set<-new_sample_set[,-grep("^i.",colnames(new_sample_set),value=T),with=F]
                      
                      },
                      "CO"={
                        new_sample_set<-indicator_converter(All_cat[eval(as.name(Groupvars))=="Comparison",names(change_vars_dist_list)[i]:=sample(x=sort(unique(eval(as.name(names(change_vars_dist_list)[i])))),size=.N,replace=T,prob=new_prob[[i]])],c(IDvars,Groupvars))[sample_set,on=c(IDvars,Groupvars)]
                        sample_set<-new_sample_set[,-grep("^i.",colnames(new_sample_set),value=T),with=F]
                        
                      },
                      "BO"={
                        new_sample_set<-indicator_converter(All_cat[,names(change_vars_dist_list)[i]:=sample(x=sort(unique(eval(as.name(names(change_vars_dist_list)[i])))),size=.N,replace=T,prob=new_prob[[i]])],c(IDvars,Groupvars))[sample_set,on=c(IDvars,Groupvars)]
                        sample_set<-new_sample_set[,-grep("^i.",colnames(new_sample_set),value=T),with=F]
                        
                      },
                      stop("You need to enter a valid group to apply the binomial distribution."))
           },
           
           "N"={
             
             switch(substr(toupper(grouplevel[1]),1,2),
                    "TR"={
                      sample_set[Group=="Treatment",names(change_vars_dist_list)[i]:=abs(rnorm(.N,mean = new_mean[i],sd=new_sd[i])),by="Group"]
                    },
                    "CO"={
                      sample_set[Group=="Comparison",names(change_vars_dist_list)[i]:=abs(rnorm(.N,mean = new_mean[i],sd=new_sd[i])),by="Group"]
                    },
                    "BO"={
                      sample_set[,names(change_vars_dist_list)[i]:=abs(rnorm(.N,mean=new_mean[i],sd=new_sd[i]))]
                    },
                    stop("You need to enter a valid group to apply the normal distribution."))
           },
           "U"={
             
             switch(substr(toupper(grouplevel[1]),1,2),
                    "TR"={
                      sample_set[Group=="Treatment",names(change_vars_dist_list)[i]:=as.numeric(ceiling(runif(.N,min=min(new_range[[i]]),max=max(new_range[[i]])))),by="Group"]
                    },
                    "CO"={
                      sample_set[Group=="Comparison",names(change_vars_dist_list)[i]:=as.numeric(ceiling(runif(.N,min=min(new_range[[i]]),max=max(new_range[[i]])))),by="Group"]
                    },
                    "BO"={
                      sample_set[,names(change_vars_dist_list)[i]:=as.numeric(ceiling(runif(.N,min=min(new_range[[i]]),max=max(new_range[[i]]))))]
                    },
                    stop("You need to enter a valid group to apply the uniform distribution."))
           },
           stop("You need to specify 'binomial','normal' or 'uniform' distribution for each variable you wish to change.")
    )
  }

  assign(as.character(pairlist_converter(as.list(match.call())$sample_set)),sample_set,envir = .GlobalEnv)
  return(print(paste0("Updated ",sQuote(pairlist_converter(as.list(match.call())$sample_set))," successfully.")))
  
}

