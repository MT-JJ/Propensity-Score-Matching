###Consolidates functions from Optmatch and Matchit and provides options to produce multiple types of matches. User can also use nonglm generated psm scores if desired but defaults to glm distances.
############Nearest Neighbor Matching##################
Nearest_Match_fn<-function(Matchingcovs_PSM,Treatvar,ID_vars,include.psm=F,distance="logit",custom_dist=F,matching_model="glm",family="binomial",weights=NULL,ratio=1,m.order="random",caliper=NULL,exact_vars=NULL,ties=T,distance.tolerance=1e-05,global=T,...){
  '%notin%'<-Negate('%in%')
  require(data.table)
  dataframe_name<-pairlist_converter(as.list(match.call())$Matchingcovs_PSM)
  X_psm<-suppressWarnings(data.table(eval(Matchingcovs_PSM))[,-c("Inverse_weight"),with=F])
  if(sum(sapply(X_psm,function(i) sum(is.na(i))))){stop("You need to remove missing values first.  Use SD_mean_before function to run diagnostics and remove missing values.")}
  
  #Checking for numeric data
  if(sum(sapply(copy(X_psm)[,-c(Treatvar,ID_vars),with=F],is.numeric)==F)!=0){stop(paste("You need to change",paste(sQuote(names(which(sapply(copy(X_psm)[,-c(Treatvar,ID_vars),with=F],is.numeric)==F))),collapse=","),"to numeric data first!"))}
  
  ll<-lapply(as.list(match.call(expand.dots = F))[c("Treatvar","ID_vars")],as.name)
  ##Checking if the psmvariable is attached to dataframe
  if(include.psm==F){
    X<-suppressWarnings(copy(X_psm)[,-c("psmscore"),with=F])
    } else {X<-copy(X_psm)}
  Treatvar_class<-class(X[[Treatvar]])
  X[[Treatvar]]<-switch(Treatvar_class,
                        "factor"={print(paste0(sQuote(ll$Treatvar)," is a factor. Assigning 0 to ",sQuote(levels(X[[Treatvar]])[1])," and 1 to ",sQuote(levels(X[[Treatvar]])[2]),". If this is incorrect assign higher level of ",sQuote(ll$Treatvar)," to the treatment and rerun."));as.numeric(X[[Treatvar]])-1
                          },
                        "numeric"={print(paste0(sQuote(ll$Treatvar)," is numeric. Assigning 0 to ",sQuote(min(X[[Treatvar]]))," and 1 to ",sQuote(max(X[[Treatvar]])),". If this is incorrect assign higher value of ",sQuote(ll$Treatvar)," to the treatment and rerun."));plyr::mapvalues(X[[Treatvar]],sort(unique(X[[Treatvar]])),c(0,1))
                          },
                        "character"={tmp_fac<-factor(ll$Treatvar);print(paste0(sQuote(X[[Treatvar]])," is a character vector. Assigning 0 to ",sQuote(levels(tmp_fac)[1])," and 1 to ",sQuote(levels(tmp_fac)[2]),". If this is incorrect convert ",sQuote(ll$Treatvar)," to a factor or numeric variable and rerun."));as.numeric(tmp_fac)-1
                        },
                        stop(print(paste("Cannot match on",sQuote(ll$Treatvar),"because it is not a factor or numeric vector.  Check to see if you assigned correct variable to",sQuote("Treatvar"))))
                        )
  if(custom_dist==F){
    if(xtabs(paste("~",Treatvar),data = X)[1]<xtabs(paste("~",Treatvar),data=X)[2]&m.order=="random"){
      print("Fewer control than treated units preventing default random matching order. Nearest Neighbor without Replacement will match Treatment units in the order of largest to smallest. Specify 'smallest' in m.order option if you wish to change the matching behavior or add a caliper to prevent potential poor matches.")
      m.order="largest"
      }
    ##GLM to pass to matchit or optmatch
    if(matching_model=="glm"){
      matching_model<-substitute(glm(paste(Treatvar,"~.",paste0("-(",c(ID_vars,exact_vars),")")),data=X,family = family,weights=weights))
    } else matching_model<-substitute(matching_model)
    
    if(is.null(caliper)){caliper=0}
    
    Matching_nearest_results<-suppressWarnings(lapply(lapply(c(TRUE,FALSE),function(i) MatchIt::matchit(eval(matching_model),data=X,ratio=ratio,m.order=m.order,caliper=caliper,exact=exact_vars,replace=i)),"[","match.matrix"))
    
    All_match_results<-lapply(setattr(Matching_nearest_results,"names",c("Nearest_with_replacement","Nearest_without_replacement")),function(i) data.table(setNames(data.frame(i),paste0(rep("Comparison",ratio),c(1:ratio))),keep.rownames = "Treat"))
    
    Matched_sets<-lapply(All_match_results,function(i) melt(copy(i)[,colnames(i):=lapply(.SD,function(j) eval(X[[ID_vars]])[match(j,attr(X,"row.names"))]),.SDcols=colnames(i)][,Match_num:=1:.N],id.vars="Match_num",variable.name=Treatvar,value.name=ID_vars)[,Match_num:=as.numeric(Match_num)])
    
    Unmatched_treatment<-lapply(Matched_sets,function(i) copy(i)[Match_num%in%copy(i)[is.na(eval(ll$ID_vars))]$Match_num][eval(ll$Treatvar)=="Treat"][,c(ID_vars,"Match_num"),with=F])
    Function_ind<-"MatchIt"
    } else{
      Function_ind<-"Matching"
      if(distance%in%colnames(X)){
        Matching_nearest_results<-setattr(lapply(c(TRUE,FALSE), function(i) Matching::Match(Tr=X[[Treatvar]],X=X[[distance]],M=ratio,ties = ties,weights=weights,caliper = caliper,distance.tolerance = distance.tolerance,replace=i)),"names",c("Nearest_with_replacement","Nearest_without_replacement"))
        
        Matched_sets<-lapply(Matching_nearest_results,function(i) unique(melt(data.table(data.frame(i[match(c("index.treated","index.control"),names(i))]))[,Match_num:=rep(seq_along(rle(index.treated)$lengths),rle(index.treated)$lengths)],id.vars="Match_num"),by=c("Match_num","value"))[,-c("variable"),with=F])
        
        Unmatched_treatment<-lapply(Matching_nearest_results,function(i) unname(unlist(copy(X)[i[["index.dropped"]],ID_vars,with=F])))
      } else {stop(paste0("Name of your 'custom distance' must appear in " ,sQuote(dataframe_name)))}
      }
  if(Function_ind=="MatchIt"){
    Final_matched_sets<-setattr(lapply(seq_along(Matched_sets),function(i) copy(X_psm)[eval(copy(Matched_sets[[i]])[,-c(Treatvar),with=F]),on=ID_vars][Match_num%notin%Unmatched_treatment[[i]]$Match_num]),"names",names(Matched_sets))
    print("Matching conducted using the 'MatchIt' package.")
    } else {
      Final_matched_sets<-setattr(lapply(seq_along(Matched_sets),function(i) copy(X)[,value:=as.numeric(row.names(.SD))][Matched_sets[[i]],on="value"][,-c("value"),with=F]),"names",names(Matched_sets))
      print("Matching conducted using the 'Matching' package.")
      }
  if(is.null(exact_vars)){
    exact_data<-suppressWarnings(copy(X)[,-c(Treatvar,distance,ID_vars),with=F])
    exact_patterns<-unique(exact_data,by=colnames(exact_data))[,Match_num:=.I]
    exact_data_IDs<-copy(X)[,ROWID:=.I][,c("ROWID",ID_vars),with=F]
    Exact_all_test<-suppressWarnings(tryCatch(Matching::Match(Tr=X[[Treatvar]],X=exact_data,exact = T,caliper=NULL),error=function(e) 1))
    if(is.list(Exact_all_test)){
      Exact_all<-rbindlist(lapply(c("index.control","index.treated"),function(i) setnames(data.table(data.frame(Exact_all_test[[i]])),"ROWID")))
      Final_matched_sets<-append(Final_matched_sets,list(Exact=unique(exact_patterns[eval(copy(X)[exact_data_IDs,on=ID_vars][Exact_all,on="ROWID"][,-c("ROWID"),with=F]),on=colnames(exact_data)],by=ID_vars)))
      Unmatched_treatment<-append(Unmatched_treatment,list(data.table(data.frame(ID=unname(unlist(copy(X)[Exact_all_test[["index.dropped"]],ID_vars,with=F]))))))
      names(Unmatched_treatment)[3]<-"Exact"
    }
  }
  
  if(global){
    basenames<-lapply(seq_along(Final_matched_sets),function(i) paste(names(Final_matched_sets)[i],as.character(dataframe_name),sep="_"))
    if(exists(basenames[[1]],envir = .GlobalEnv)){
      basenames<-Global_name_gen(names(Final_matched_sets),dataframe_name)
      }
    lapply(seq_along(Final_matched_sets),function(i){
      assign(basenames[[i]],Final_matched_sets[[i]],envir = .GlobalEnv)
    })
  }
  
  ###Printing quick summary to console
  lapply(seq_along(Final_matched_sets),function(i) print(paste(nrow(Unmatched_treatment[[i]]),"treatment subjects removed and",nrow(unique(copy(Final_matched_sets[[i]]),by=ID_vars)),"total unique subjects in the",names(Final_matched_sets)[i],"matched set.")))
  return(list(Matched_Sets=Final_matched_sets,Unmatched_IDs=Unmatched_treatment,Raw_function_output=Matching_nearest_results))
}

##############Optimal Matching Function##########################
Optimal_match_fn<-function(Matchingcovs_PSM,Treatvar,ID_vars,ratio=1,fitted_glm=T,glm_family="binomial",glmweights=NULL,other_fitted_model=NULL,custom_dist=NULL,Match_type=c("both","Pair","Full"),exact_vars=NULL,include_exact_psm=F,caliper=NULL,Caliper_vars=NULL,Caliper_glm=T,dist_method=c("euclidean","mahalanobis","rank_mahalanobis","custom"),standardization.scale=sd,min.controls = 0, max.controls = Inf,omit.fraction=NULL,global=T,verbose=T){
  
  if(global){
    suppressPackageStartupMessages(require(optmatch))
    } else require(optmatch)
  
  dataframe_name<-pairlist_converter(as.list(match.call())$Matchingcovs_PSM)
  
  X<-suppressWarnings(data.table(eval(Matchingcovs_PSM))[,-c("Inverse_weight"),with=F])
  ll_calls<-as.list(match.call(expand.dots = F))
  ll<-lapply(ll_calls[c("Treatvar","ID_vars")],as.name)
  
  ###Processing Treatment variable by converting it to factor and then assigning numeric values (1 for treatment, 0 for comparison)
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
  
  orig_treat_levels<-setNames(levels(X[[Treatvar]]),c("Comparison","Treatment"))##Link between original levels and comparison and treatment.
  X<-copy(X)[,eval(Treatvar):=as.numeric(eval(ll$Treatvar))-1]
  
  X_optmatch<-suppressWarnings(X[,-c("psmscore"),with=F])
  
  #Match_on syntax
  ##1. This nested switch statament allows for all valid classes of "x" for match_on.
  if(is.null(other_fitted_model)){
    glm_method<-switch(dist_method[1],
                       "euclidean"={
                         if(!include_exact_psm){
                           substitute(glm(paste(eval(Treatvar),"~",paste0(colnames(X_optmatch)[-match(c(Treatvar,ID_vars,exact_vars),colnames(X_optmatch))],collapse = "+")),data=X_optmatch,family=glm_family,weights=glmweights))
                         } else {
                             substitute(glm(as.formula(paste(Treatvar,"~",paste0(colnames(X_optmatch)[-match(c(Treatvar,ID_vars),colnames(X_optmatch))],collapse = "+"))),data=X_optmatch,family=glm_family,weights=glmweights))
                           }
                         },
                       "mahalanobis"={
                         as.formula(paste(Treatvar,"~",paste0(colnames(X_optmatch)[-match(c(ID_vars,exact_vars,Treatvar),colnames(X_optmatch))],collapse="+")))
                         },
                       "rank_mahalanobis"={
as.formula(paste(Treatvar,"~",paste0(colnames(X_optmatch)[-match(c(ID_vars,exact_vars,Treatvar),colnames(X_optmatch))],collapse="+")))
                         },
                       "custom"={
                         if(!is.null(custom_dist)){
    list(x=setNames(X_optmatch[[custom_dist]],X_optmatch[[ID_vars]]),z=X_optmatch[[Treatvar]])
                         }
                         else {
                           stop(paste0("Entering custom implies you created your own distance measure outside this function and attached that measure to your data set.  Please include the name of this variable in ",dataframe_name,"."))
                          }
                         },

                          stop("You need to specifiy a valid dist-method in the function call!"))
    } else {
      glm_method<-other_fitted_model
      }
  
  ##2. Creates caliper specification
  if(!is.null(Caliper_vars)){
    if(Caliper_glm){
      Caliper_glm_model<-substitute(glm(paste(Treatvar,"~",paste0(Caliper_vars,collapse = "+")),data=X_optmatch,family=glm_family,weights=glmweights))
      
      Caliper_model<-suppressWarnings(caliper(eval(match_on(eval(Caliper_glm_model),data=X_optmatch,standardization.scale=standardization.scale)),width = caliper))
      } else{
        Caliper_model<-caliper(match_on(as.formula(paste(Treatvar,"~",paste0(Caliper_vars,collapse = "+"))),data=X_optmatch),width = caliper)}##Possibly add an else statement for mahalanobis all vars
    }
  
  ##3. Creates exact matches
  if(!is.null(exact_vars)){
    Within_model<-as.formula(paste(Treatvar,"~",paste0(exact_vars,collapse="+")))
    }
  ####This switch combines all possible forms of syntax options for the match_on matrix which will then be passed onto pair match or full match.
  exact_cal_var_ind<-as.character(interaction(is.null(exact_vars),is.null(Caliper_vars)))
  match_on_matrix<-switch(exact_cal_var_ind,
                          "TRUE.TRUE"={
                            if(is.null(custom_dist)) {
                              match_on(eval(glm_method),data=X_optmatch,caliper=caliper,standardization.scale=standardization.scale)
                            } else{
                                match_on(x=glm_method$x,data=X_optmatch,z=glm_method$z,caliper=caliper,standardization.scale=standardization.scale)
                              }
                            },
                          "TRUE.FALSE"={
                            if(is.null(custom_dist)) {
                            match_on(eval(glm_method),data=X_optmatch,standardization.scale=standardization.scale)+Caliper_model
                            } else {
                                match_on(x=glm_method$x,data=X_optmatch,z=glm_method$z,standardization.scale=standardization.scale)+Caliper_model}
                            },
                          "FALSE.TRUE"={
                            if(is.null(custom_dist)) {
                              match_on(eval(glm_method),data=X_optmatch,within = exactMatch(Within_model,data=X_optmatch),caliper=caliper,standardization.scale=standardization.scale)
                            } else {
                                match_on(x=glm_method$x,data=X_optmatch,z=glm_method$z,within = exactMatch(Within_model,data=X_optmatch),caliper=caliper,standardization.scale=standardization.scale)
                              }
                            },
                          "FALSE.FALSE"={
                            if(is.null(custom_dist)) {
                              match_on(eval(glm_method),data=X_optmatch,within=exactMatch(Within_model,data=X_optmatch),standardization.scale=standardization.scale)+Caliper_model
                            } else {
                                match_on(x=glm_method$x,data=X_optmatch,z=glm_method$z,within = exactMatch(Within_model,data=X_optmatch),standardization.scale=standardization.scale)+Caliper_model
                              }
                            })
  
  matches<-switch(Match_type[1],
                  "both"={
                    list(Full=full(match_on_matrix,data=X_optmatch,min.controls=min.controls,max.controls=max.controls,omit.fraction=omit.fraction),Pair=pair(match_on_matrix,controls=ratio,data=X_optmatch))
                    },
                  "Pair"={
                    list(Pair=pair(match_on_matrix,controls = ratio,data=X_optmatch))
                    },
                  "Full"={
                    list(Full=full(match_on_matrix,data=X_optmatch,min.controls=min.controls,max.controls=max.controls,omit.fraction=omit.fraction))
                    },
stop("Must enter a valid type of matching!"))
  
  final_matches_all<-tryCatch(
      {
        lapply(matches,function(i) data.table(data.frame(X,Match_num=as.numeric(plyr::mapvalues(i,levels(i),1:nlevels(i))),Match_status=matched(i)))[,eval(Treatvar):=factor(eval(ll$Treatvar),labels = c("Comparison","Treatment"))])
      },error=function(e) 1)
    
    final_matches<-lapply(final_matches_all,function(i) i[Match_status==T])
    Unmatched_treat<-lapply(final_matches_all,function(i) i[Match_status==F&eval(ll$Treatvar)=="Treatment",ID_vars,with=F])
    
    #Printing quick summary to console
    lapply(seq_along(final_matches),function(i) print(paste(nrow(Unmatched_treat[[i]]),"treatment subjects removed and",length(unique(final_matches[[i]][[ID_vars]])),"total unique subjects in the",names(final_matches)[i],"matched set.")))
  
  if(global){
    basenames<-lapply(seq_along(final_matches),function(i) paste(names(final_matches)[i],as.character(dataframe_name),sep="_"))
    
    if(exists(basenames[[1]],envir = .GlobalEnv)){
      basenames<-Global_name_gen(names(final_matches),dataframe_name)
      }
    lapply(seq_along(final_matches),function(i){
      assign(basenames[[i]],final_matches[[i]],envir = .GlobalEnv)
      })
    }
  
  return(list(Matched_sets=final_matches,UnmatchedIDs=Unmatched_treat))
}

