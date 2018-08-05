PSM_diagnostics<-function(Matchingcovs,Treatvar,ID_vars,convert_ind=T,ref_vars=NULL,outcome_vars=NULL,missing_check=F,family="binomial",plot.points=T,pretty_print=F,add_new_file=T,psm=T,OR=F,global=T,verbose=T,...)
{
  '%notin%'<-Negate('%in%')
  options(useFancyQuotes = F)
  if(global){
    suppressPackageStartupMessages(require(openxlsx))
    suppressPackageStartupMessages(require(data.table))
    suppressPackageStartupMessages(require(lattice))
    suppressPackageStartupMessages(require(latticeExtra))
    suppressPackageStartupMessages(require(gridExtra))
    suppressPackageStartupMessages(require(maptools))
    suppressPackageStartupMessages(require(DescTools))
  } else {
    require(openxlsx)
    require(data.table)
    require(lattice)
    require(latticeExtra)
    require(gridExtra)
    require(maptools)
    require(DescTools)
  }
  dataframe_name<-as.character(pairlist_converter(as.list(match.call())$Matchingcovs))
  #####Checking for missing values imputation####
  ##Not allowing missing values for Treatment variable and ID variable
  if(sum(sapply(data.table(Matchingcovs)[,c(ID_vars,Treatvar),with=F],function(i) sum(is.na(i))))!=0){
    if(verbose){
      cat("\nThe following variable(s) have missing values:")
      char_quote_fn(names(which(sapply(data.table(Matchingcovs)[,c(ID_vars,Treatvar),with=F],function(i) sum(is.na(i)))==0)))
    }
    stop("Cannot have missing values in your ID or Treatment variable",call. = F)
  }
  ##Imputation function##
  if(missing_check!=F){
    impute_method<-switch(as.character(missing_check),
                          "TRUE"={"R"},
                          {NULL})
    
    Matchingcovs<-Impute_missing_fn(Matchingcovs,impute_method = impute_method,ignore_vars = c(ID_vars,Treatvar),verbose=verbose)
  }
  #####Control flow statements for processing categorical variables####
  if(is.null(ref_vars)){
    if(convert_ind){##No ref vars provided but user wants to convert##
      tmp_X<-data.table(Matchingcovs)[,-c(ID_vars,Treatvar),with=F][,lapply(copy(.SD),function(i) if(is.character(i)|is.factor(i)){
        as.factor(i)
      } else{
        as.numeric(i)
      })]
      categorical_vars<-names(which(sapply(tmp_X,is.factor)))
      if(length(categorical_vars)>0){
        ref_vars<-paste0(categorical_vars,unlist(lapply(categorical_vars,function(i) levels(copy(tmp_X)[[i]])[1])))
        Matchingcovs<-indicator_converter(data.table(Matchingcovs),IDvars = c(Treatvar,ID_vars),verbose=verbose)
        if(verbose){
          print(paste0("Creating ",paste0(ref_vars,collapse=", ")," as the reference variables"))
        }
        
        setnames(Matchingcovs,make.names(colnames(Matchingcovs)))###Matchingcovs contains all variables and X no reference variables
        X<-copy(Matchingcovs)[,-ref_vars,with=F]
      } else{
        X<-copy(Matchingcovs)
      }
    } else {
      print("Assuming no reference variables are included in 'Matchingcovs.'  Diagnostics will not contain stats on reference variables.")
      Matchingcovs<-setnames(data.table(Matchingcovs),make.names(colnames(Matchingcovs)))
      X<-copy(Matchingcovs)
    }
  } else {
    if(convert_ind){####Ref vars indicated and user wants to convert####
      Matchingcovs<-indicator_converter(data.table(Matchingcovs),IDvars = c(Treatvar,ID_vars))
      ref_vars<-grep(paste0(ref_vars,collapse="|"),colnames(Matchingcovs),value = T)
      ref_call<-pairlist_converter(match.call()$ref_vars)
      if(length(ref_vars)!=length(ref_call)){
        ref_check<-setattr(lapply(ref_call,function(i) sum(grepl(i,colnames(Matchingcovs)))),"names",ref_call)
        if(verbose){
          cat("\nReference variable(s) not found:")
          char_quote_fn(names(which(ref_check==0)))
          cat("\n")
        }
        stop("The above reference variable name(s) not found in any of the categorical variables. Check spelling.",call. = F)
      } else {
        Matchingcovs<-setnames(data.table(Matchingcovs),make.names(colnames(Matchingcovs)))
        X<-copy(Matchingcovs)[,-ref_vars,with=F]
      }
    } else {##Ref vars indicated and user already converted.
      Matchingcovs<-setnames(data.table(Matchingcovs),make.names(colnames(Matchingcovs)))
      ref_vars<-make.names(ref_vars)
      X<-copy(Matchingcovs)[,-ref_vars,with=F]
    }
  }
  ######Separating Outcome Vars if any#####
  if(!is.null(outcome_vars)){
    Outcome_vars<-copy(Matchingcovs)[,c(ID_vars,Treatvar,outcome_vars),with=F]
    X<-X[,-outcome_vars,with=F]
  }
  #####Processing Treatment Variable#########
  ll<-lapply(as.list(match.call(expand.dots = F))[grepl("Treatvar|ID_vars|psm",names(match.call()))],as.name)
  if(verbose){
    if(is.numeric(X[[Treatvar]])){
      print(paste(sQuote(match.call()$Treatvar),"is not of class factor.  Assuming",sQuote(max(X[[Treatvar]])),"is the treatment group."))
    } else if(is.character(X[[Treatvar]])){
      print(paste(sQuote(match.call()$Treatvar),"is not of class factor.  Assuming",sQuote(sort(unique(X[[Treatvar]]))[length(unique(X[[Treatvar]]))]),"is the treatment category based on alphanumeric order. If this is incorrect either convert class to a factor specifying proper treatment variable or recode the treatment variable such that the treatment code is the highest alphanumeric value."))
    } else {
      print(paste(sQuote(match.call()$Treatvar),"is of class factor and assuming",sQuote(levels(X[[Treatvar]])[nlevels(X[[Treatvar]])]),"is the treatment group."))
    }
  }
  
  X[[Treatvar]]<-cat2facNum_converter_fn(X[[Treatvar]],verbose=F)##Converted to factor
  
  orig_treat_levels<-setNames(levels(X[[Treatvar]]),c("Comparison","Treatment"))##Link between orginal levels and comparison and treatment.
  levels(X[[Treatvar]])<-c("Comparison","Treatment")
  
  ######Creating glm object based on call to PSM##########
  if(psm==T){
    glm_call<-substitute(paste(Treatvar,"~.",paste0("-(",paste0(c(ID_vars),collapse = "+"),")")))
    glmtest<-suppressWarnings(glm(eval(glm_call),data=X,family = family))
    pseudo_test<-listrep2datatable_fn(PseudoR2(do.call("glm",list(formula=paste(Treatvar,"~."),data=copy(X)[,c(Treatvar,names(glmtest$coefficients)[-1]),with=F],family=family)),which = "all"),left_varname = "Statistic",right_varname = "Value")[,Value:=format(round(Value,digits = 4),nsmall = 4)]
    if(is.null(ref_vars)){
      tmp_data<-data.table(data.frame(glmtest$model,psmscore=predict(glmtest,type="response")))[,Inverse_weight:=ifelse(as.numeric(eval(ll$Treatvar))-1==1,1/psmscore,1/(1-psmscore))]
    } else {
      tmp_data<-data.table(data.frame(glmtest$model,psmscore=predict(glmtest,type="response")))[,Inverse_weight:=ifelse(as.numeric(eval(ll$Treatvar))-1==1,1/psmscore,1/(1-psmscore))][Matchingcovs[,c(ID_vars,ref_vars),with=F],on=ID_vars]
    }
    
    if(nrow(glmtest$model)!=nrow(glmtest$data)){
      stop(paste(length(glmtest$na.action),"observations missing values on your matching covariates. Set 'missing_check=T' and rerun."),call. = F)
    } else{final_datapsm<-tmp_data}
  } else if(psm==F){
    final_datapsm<-tryCatch(copy(X)[Matchingcovs[,c(ID_vars,ref_vars),with=F],on=ID_vars],error=function(e) copy(X))
  } else{
    final_datapsm<-tryCatch(copy(X)[,Inverse_weight:=ifelse(as.numeric(eval(ll$Treatvar))==1,1/eval(ll$psm),1/(1-eval(ll$psm)))][Matchingcovs[,c(ID_vars,ref_vars),with=F],on=ID_vars],error=function(e) copy(X)[,Inverse_weight:=ifelse(as.numeric(eval(ll$Treatvar))==1,1/eval(ll$psm),1/(1-eval(ll$psm)))])
  }
  ####Calculating Descriptive Statistics###########
  cat_vars<-Binary_indication_fn(final_datapsm)
  IndicatorNs<-data.table(t(copy(final_datapsm)[,c(Treatvar,cat_vars),with=F][,lapply(.SD,sum),by=Treatvar,.SDcols=c(cat_vars)]),keep.rownames = "Covariates")
  IndicatorNs<-setcolorder(setnames(IndicatorNs,c("V1","V2"),paste("N",unlist(IndicatorNs[1,2:3,with=F]),sep="_"))[-1,],c(1,3,2))
  ContNs<-setnames(dcast(melt(copy(final_datapsm)[,lapply(.SD,function(i) sum(!is.na(i))),by=Treatvar,.SDcols=which(colnames(final_datapsm)%notin%c(cat_vars,ID_vars,Treatvar,"Inverse_weight"))],id.vars=Treatvar),paste0("variable~",Treatvar)),colnames(IndicatorNs))
  
  Group_means_sds<-rbind(IndicatorNs,ContNs)[eval(setnames(dcast(melt(copy(final_datapsm)[,-c(ID_vars,"Inverse_weight"),with=F][,lapply(copy(.SD),function(i) list(Mean=mean(i,na.rm=T),SD=sd(i,na.rm=T))),by=Treatvar][,Stat:=rep(c("Mean","SD"),2)],id.vars = c(Treatvar,"Stat")),paste0("variable~Stat+",Treatvar)),"variable","Covariates")),on="Covariates"][,lapply(.SD,function(i) as.vector(i,mode = "character"))][,paste(rep(c("N","Mean","SD"),c(2,2,2)),rep(c("Comparison","Treatment"),3),sep="_"):=lapply(.SD,function(i) as.numeric(i)),.SDcols=paste(rep(c("N","Mean","SD"),c(2,2,2)),rep(c("Comparison","Treatment"),3),sep="_")]
  
  SD_list<-rbind(Group_means_sds[Covariates%notin%cat_vars,grep("^SD_|Covariates",colnames(Group_means_sds)),with=F],Group_means_sds[Covariates%in%cat_vars,c("Covariates","Mean_Comparison","Mean_Treatment"),with=F][,c("SD_Comparison","SD_Treatment"):=lapply(.SD,function(i) sqrt(i*(1-i))),.SDcols=c("Mean_Comparison","Mean_Treatment")][,-c("Mean_Comparison","Mean_Treatment"),with=F])##Need to sqrt this to make consistent with continuous formula. SD vs variance for binomials.
  
  ####Original Balance Results Calculation#####
  final_results<-copy(Group_means_sds)[,Raw_Bias:=abs(Mean_Treatment-Mean_Comparison)][,Pooled_SD:=sqrt((SD_Treatment^2+SD_Comparison^2)/2)][,SD_Bias:=Raw_Bias/Pooled_SD]
  
  Balance_object<-copy(final_results)[,c("Covariates","SD_Bias"),with=F]
  class(Balance_object)<-append(class(Balance_object),"Balance")
  #############Creating Diagnostic Graphs for PSM variables if indicated by user################
  ###Preparing PSM Summary as printed to console#####
  if(psm!=F){
    PSM_All<-tryCatch(copy(final_datapsm)[is.na(psmscore)==F,list(N=.N,Mean=mean(psmscore,na.rm=T),SD=sd(psmscore,na.rm=T),Median=median(psmscore,na.rm=T),Range=(paste0(format(round(range(psmscore,na.rm=T),4),nsmall = 4),collapse = "-")))],error=function(e) copy(final_datapsm)[is.na(eval(ll$psm))==F,list(N=.N,Mean=mean(eval(ll$psm),na.rm=T),SD=sd(eval(ll$psm),na.rm=T),Median=median(eval(ll$psm),na.rm=T),Range=(paste0(format(round(range(eval(ll$psm),na.rm=T),4),nsmall = 4),collapse = "-")))])
    
    
    PSM_Treat<-tryCatch(copy(final_datapsm)[is.na(psmscore)==F,list(N=.N,Mean=mean(psmscore,na.rm=T),SD=sd(psmscore,na.rm=T),Median=median(psmscore,na.rm=T),Range=(paste0(format(round(range(psmscore,na.rm=T),4),nsmall = 4),collapse = "-"))),by=Treatvar][,c(Treatvar,"Range"):=lapply(.SD,as.character),.SDcols=c(Treatvar,"Range")],error=function(e) copy(final_datapsm)[is.na(eval(ll$psm))==F,list(N=.N,Mean=mean(eval(ll$psm),na.rm=T),SD=sd(eval(ll$psm),na.rm=T),Median=median(eval(ll$psm),na.rm=T),Range=(paste0(format(round(range(eval(ll$psm),na.rm=T),4),nsmall = 4),collapse = "-"))),by=Treatvar][,c(Treatvar,"Range"):=lapply(.SD,as.character),.SDcols=c(Treatvar,"Range")])
    
    PSM_summary<-rbind(PSM_All,PSM_Treat,fill=T)[,colnames(PSM_Treat)[-match(c(Treatvar,"Range"),colnames(PSM_Treat))]:=lapply(.SD,function(i) round(i,4)),.SDcols=colnames(PSM_Treat)[-match(c(Treatvar,"Range"),colnames(PSM_Treat))]]
    
    ####Creating data and color format for graphs####
    
    PSM_graph_data<-tryCatch(copy(final_datapsm[is.na(psmscore)==F]),error=function(e) copy(final_datapsm[is.na(eval(ll$psm))==F]))
    plotcolors<-c("cyan","magenta")
    
    ####Histograms#######
    PSM_hist<-tryCatch(histogram(as.formula(paste("~","psmscore","|",Treatvar)),data=eval(final_datapsm[is.na(psmscore)==F]),group=eval(ll$Treatvar),type="count",main=paste("Predicted PSM Scores by",sQuote(ll$Treatvar)),xlab=c("PSM Scores"),col=plotcolors,strip=strip.custom(bg="lightgrey",par.strip.text=list(cex=0.7,fontface="bold")),scales=list(x="same",y="free",alternating=F,tck=c(1,0)),panel=function(x,col=col,...){
      panel.histogram(x,col=col[packet.number()],...)
    }),error=function(e) histogram(as.formula(paste("~",psm,"|",Treatvar)),data=eval(final_datapsm[is.na(eval(ll$psm))==F]),group=eval(ll$Treatvar),type="count",main=paste("Predicted PSM Scores by",sQuote(ll$Treatvar)),xlab=c("PSM Scores"),col=plotcolors,strip=strip.custom(bg="lightgrey",par.strip.text=list(cex=0.7,fontface="bold")),scales=list(x="same",y="free",alternating=F,tck=c(1,0)),panel=function(x,col=col,...){
      panel.histogram(x,col=col[packet.number()],...)
    }))
    
    ####Density Plots#######
    PSM_dens<-tryCatch(densityplot(as.formula(paste("~","psmscore")),data=PSM_graph_data,group = eval(ll$Treatvar),auto.key = list(border=T,background="lightgrey",fontface="bold",cex=0.7,corner=c(1,1)),main=paste("Density Plot of Predicted PSM Scores by",sQuote(ll$Treatvar)),xlab=c("PSM Scores"),scales = list(tck=c(1,0)),plot.points=plot.points),error=function(e) densityplot(as.formula(paste("~",psm)),data=PSM_graph_data,group = eval(ll$Treatvar),auto.key = list(border=T,background="lightgrey",fontface="bold",cex=0.7,corner=c(1,1)),main=paste("Density Plot of Predicted PSM Scores by",sQuote(ll$Treatvar)),xlab=c("PSM Scores"),scales = list(tck=c(1,0)),plot.points=plot.points))
    
    
    grid.arrange(PSM_hist,PSM_dens,ncol=1)
    
    #####Inverse Weights Summary######
    Inverse_dens<-densityplot(as.formula(paste("~","Inverse_weight")),data=PSM_graph_data,group = eval(ll$Treatvar),auto.key = list(border=T,background="lightgrey",fontface="bold",cex=0.7,corner=c(1,1)),main=paste("Density Plot of Inverse Weights by",sQuote(ll$Treatvar)),xlab=c("PSM Weights"),scales = list(tck=c(1,0)),plot.points=plot.points)
  }
  
  ##Allows for sloppy data.table syntax
  if(pretty_print==T){
    dataframe_name<-pairlist_converter(match.call()$Matchingcovs)
    tmpfile<-paste0("Descriptive Stats for ",dataframe_name,".xlsx")
    tmpsheet<-"DemoStats"
    tmpsheet2<-"BalanceBefore"
    if(add_new_file==T){
      Pretty_excel_print(copy(Group_means_sds)[,-grep("^i.",colnames(Group_means_sds)),with=F],var_col = "Covariates",grouped_col = 2,filename = tmpfile ,sheetname=tmpsheet,verbose=verbose)
    } else {
      Pretty_excel_print(copy(Group_means_sds)[,-grep("^i.",colnames(Group_means_sds)),with=F],var_col = "Covariates",grouped_col = 2,filename = tmpfile,sheetname=tmpsheet,append = T,verbose=verbose)
    }
    
    invisible(Pretty_excel_print(copy(final_results)[,c("Pooled_SD","SD_Bias"):=lapply(.SD,function(i) round(i,4)),.SDcols=c("Pooled_SD","SD_Bias")],var_col = "Covariates",grouped_col = 2,filename = paste0("Descriptive Stats for ",dataframe_name,".xlsx"),sheetname=tmpsheet2,append = T,verbose=verbose))
  }
  
  ###Assigning variables to global environment for user if desired#####
  if(global==T){
    dataname<-as.character(pairlist_converter(as.list(match.call())$Matchingcovs))
    
    if(!is.null(ref_vars)){
      assign(paste(dataname,"PSM",sep="_"),suppressWarnings(eval(copy(final_datapsm)[,-ref_vars,with=F])),envir = .GlobalEnv)
    } else {
      assign(paste(dataname,"PSM",sep="_"),suppressWarnings(eval(copy(final_datapsm)[,with=F])),envir = .GlobalEnv)
    }
    if(!is.null(ref_vars)){
      assign(paste(dataname,"ref",sep="_"),tryCatch({
        copy(final_datapsm)[,c(ID_vars,Treatvar,ref_vars),with=F]
      },error=function(e) print(paste0("Cannot create  dataframe for reference variables."))),envir = .GlobalEnv)
    } else print(paste0("Cannot create  dataframe for reference variables."))
    
    assign(paste(dataname,"outcome",sep="_"),tryCatch(Outcome_vars,error=function(e) {
      NULL
      if(verbose){
        print(paste0("Cannot create  dataframe for outcome variables because none indicated."))
      }}),envir = .GlobalEnv)
    
    assign(paste(dataname,"desc_stats",sep="_"),copy(Group_means_sds)[,which(colnames(Group_means_sds)%notin%c("Covariates","N_Comparison","N_Treatment")):=lapply(.SD,function(i) format(round(i,4),nsmall = 4)),.SDcols=which(colnames(Group_means_sds)%notin%c("Covariates","N_Comparison","N_Treatment"))],envir = .GlobalEnv)
    
    assign(paste(dataname,"PooledSD",sep="_"),eval(copy(final_results)[,c("Covariates","Pooled_SD"),with=F]),envir = .GlobalEnv)
    
    assign(paste(dataname,"SDmeans",sep="_"),Balance_object,envir = .GlobalEnv)
  }
  
  ################Printing and saving all results to list based on how PSM was entered by user########
  if(psm!=F){
    if(psm==T){
      cat("\n\t\tSummary of Propensity Scores\n------------------------------------------------------------------------------------\n")
      print(eval(PSM_summary[is.na(eval(ll$Treatvar)),eval(Treatvar):="All"]))
      cat("\n\tGoodness of Fit Statistics for Logistic Regression\n------------------------------------------------------------------------------------\n")
      print(pseudo_test)###Need to change to if else with no ref
      
      return(suppressWarnings(list(PSM_data=final_datapsm[,-ref_vars,with=F],Descriptive_stats=Group_means_sds[,which(colnames(Group_means_sds)%notin%c("Covariates","N_Comparison","N_Treatment")):=lapply(.SD,function(i) format(round(i,4),nsmall = 4)),.SDcols=which(colnames(Group_means_sds)%notin%c("Covariates","N_Comparison","N_Treatment"))],All_calculations=final_results,Pooled_SD=copy(final_results)[,c("Covariates","Pooled_SD"),with=F],SD_Means=Balance_object,PSM_summary=eval(PSM_summary[is.na(eval(ll$Treatvar)),eval(Treatvar):="All"]),Regression_output=glmtest,Inverse_weights=Inverse_dens,Outcome_vars=tryCatch(Outcome_vars,error=function(e) NULL),Ref_vars=tryCatch(copy(final_datapsm)[,c(ID_vars,Treatvar,ref_vars),with=F]))))
    } else {
      cat("\n\t\tSummary of Propensity Scores\n------------------------------------------------------------------------------------\n")
      print(eval(PSM_summary[is.na(eval(ll$Treatvar)),eval(Treatvar):="All"]))
      
      
      return(suppressWarnings(list(PSM_data=tryCatch(copy(final_datapsm)[,-ref_vars,with=F],error=function(e) copy(final_datapsm)),Descriptive_stats=Group_means_sds[,which(colnames(Group_means_sds)%notin%c("Covariates","N_Comparison","N_Treatment")):=lapply(.SD,function(i) format(round(i,4),nsmall = 4)),.SDcols=which(colnames(Group_means_sds)%notin%c("Covariates","N_Comparison","N_Treatment"))][Covariates=="psmscore",Covariates:=psm],All_calculations=final_results[Covariates=="psmscore",Covariates:=psm],Pooled_SD=copy(final_results)[,c("Covariates","Pooled_SD"),with=F],SD_Means=Balance_object,PSM_summary=eval(PSM_summary[is.na(eval(ll$Treatvar)),eval(Treatvar):="All"]),Inverse_weights=Inverse_dens,Outcome_vars=tryCatch(Outcome_vars,error=function(e) NULL),Ref_vars=tryCatch(copy(final_datapsm)[,c(ID_vars,Treatvar,ref_vars),with=F]))))
    }
  } else {
    return(list(PSM_data=final_datapsm,Descriptive_stats=Group_means_sds[,which(colnames(Group_means_sds)%notin%c("Covariates","N_Comparison","N_Treatment")):=lapply(.SD,function(i) format(round(i,4),nsmall = 4)),.SDcols=which(colnames(Group_means_sds)%notin%c("Covariates","N_Comparison","N_Treatment"))],All_calculations=final_results,Pooled_SD=copy(final_results)[,c("Covariates","Pooled_SD"),with=F],SD_Means=Balance_object,Outcome_vars=tryCatch(Outcome_vars,error=function(e) NULL),Ref_vars=tryCatch(copy(final_datapsm)[,c(ID_vars,Treatvar,ref_vars),with=F])))
  }
}
