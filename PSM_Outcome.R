PSM_Outcome<-function(Finalset,Treatvar,Outcome,Combine_out=F,ID_vars=NULL,ignore_vars=NULL,Orig_data_psm=NULL,inv_weights=F,rematch=F,Match_num="Match_num",PSM_var="psmscore",verbose=T,...){
  require(data.table)
  require(Matching)
  require(survey)
  require(rbounds)
  dataframe_name<-pairlist_converter(as.list(match.call())$FinalSet)
  ll<-lapply(list(Treatvar=eval(Treatvar),PSM_var=eval(PSM_var)),as.name)
  X<-data.table(Finalset)[,-ignore_vars,with=F]
  if(Combine_out){
    X_Outcome<-tryCatch(copy(X)[data.table(Outcome),on=c(ID_vars,Treatvar)],error=function(e) stop("You indicated you wanted to join outcome variables to matched data but need to provide an ID variable connected to the outcome in order to join them. Maybe remove it from your ignore_vars argument if you included it there and re-run.",call. = F))
    X<-X_Outcome[,-ID_vars,with=F]
    Outcome<-Outcome[,-c(ID_vars,Treatvar),with=F]
  }
  X[[Treatvar]]<-cat2facNum_converter_fn(X[[Treatvar]],verbose=verbose)##Converted to factor
  orig_treat_levels<-setNames(levels(X[[Treatvar]]),c("Comparison","Treatment"))##Link between orginal levels and comparison and treatment.
  X[[Treatvar]]<-as.numeric(X[[Treatvar]])-1
  if(sum(grepl("data",class(Outcome)))>0){
    Outcome<-suppressWarnings(data.table(Outcome)[,-c(Treatvar,ID_vars),with=F])
    outcome_vars<-colnames(Outcome)
  } else if(is.character(Outcome)){
    outcome_vars<-Outcome
  } else{
    stop("Outcome must be specified as either a data frame object consisting of outcome variables, treatment variable and ID numbers to connect to matched data or a character vector of names identifying the outcome variables in your final data set.",call. = F)
  }
  #####Weighting Analysis for ATE and ATT based on Original Data Set#############
  if(inv_weights==T){
    Covariates<-colnames(X)[which(colnames(X)%notin%c(Match_num,PSM_var,outcome_vars))]
    print(Covariates)
    if(is.null(Orig_data_psm)){###If FinalSet is the Orig_data
      X_weight<-copy(X)[,w.ate:=ifelse(ll$Treatvar==1,1/eval(ll$PSM_var),1/(1-eval(ll$PSM_var)))][,w.att:=ifelse(ll$Treatvar==1,1,1/eval(ll$PSM_var))]
      design.ate<-svydesign(ids=~1,weights = ~w.ate,data=X_weight)
      design.att<-svydesign(ids=~1,weights = ~w.att,data=X_weight)
    } else { ###If user wants both Inverse weighting and other match methods
      Orig_data<-data.table(Orig_data_psm)
      Covariates<-colnames(Orig_data)[which(colnames(Orig_data)%notin%c(outcome_vars,ignore_vars,PSM_var))]
      Orig_data[,w.ate:=ifelse(as.numeric(eval(ll$Treatvar))==1,1/ll$PSM_var,1/(1-ll$PSM_var))][,w.att:=ifelse(as.numeric(eval(ll$Treatvar))==1,1,1/ll$PSM_var)]
      design.ate<-svydesign(ids=~1,weights = ~w.ate,data=Orig_data)
      design.att<-svydesign(ids=~1,weights = ~w.att,data=Orig_data)
    }
    ATE_glm<-setattr(lapply(outcome_vars,function(i) svyglm(paste(i,"~",paste0(Covariates,collapse = "+")),design=design.ate)),"names",outcome_vars)
    ATT_glm<-setattr(lapply(outcome_vars,function(i) svyglm(paste(i,"~",paste0(Covariates,collapse = "+")),design=design.att)),"names",outcome_vars)
  }
  return(list(ATE=ATE_glm,ATT=ATT_glm))
}