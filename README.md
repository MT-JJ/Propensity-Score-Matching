# Propensity-Score-Matching-in-R
These functions depend on using the Data Wrangling Functions in R

*These PSM functions are an attempt to consolidate and unify the current propensity score matching packages to provide a guide for practitioners interested in conducting propensity score matching using R.*  

**The following propensity score matching packages are required to use these functions:**  
1) MatchIt   
2) optmatch  
3) Matching  

**Additional dependencies include:**  
1) data.table  
2) plyr  
3) lattice  
4) latticeExtra  
5) gridExtra  
6) maptools  
7) DescTools  
8) openxlsx  
9) devtools  

For Windows users Rtools must be installed and make sure instructions are followed from this website: https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows
Specifically, note that PATH must be checked during the wizard setup even though it is not the default option.

*The logic behind the functions reflect three general phases  used during the process.*

**Phase 1:** Diagnostics phase where descriptive statistics for the treatment and comparison groups are analyzed before matching. During this phase logistic regression is conducted to obtain predicted values based on the covariates included in the model.

**Phase 2:** Matching phase where different techniques can be used to produce the best matches based on the covariates provided in Phase 1.

**Phase 3:** Assessment of the matched sets to determine if baseline equivalency is obtained. Note that outcome variables are unnecessary throughout these phases and discouraged to avoid cherry-picking the matches. Instead the determination of the method used should be based on the results of balance on the background characteristics.  Once your sample is chosen during Phase 3 you can then attach the desired outcome variables to this sample and conduct your outcome analysis. Ideally you will be able to run t-tests but perfect balance may not be possible depending on your population. Hence, controlling for some background characteristics in subsequent analysis may be necessary.


*Note: While these functions can each be utilized individually they are optimized to work together and the function arguments are similar across all functions. Hence, the Diagnostics function will prepare your data set to work seamlessly with the matching functions and then balance function.*


#Function details and explanations:

##Phase 1: PSM_diagnostics Function

###Default call:
*PSM_diagnostics<-function(Matchingcovs, Treatvar, ID_vars, convert_ind=T, ref_vars=NULL, missing_check=F, family="binomial", plot.points=T, pretty_print=F, add_new_file=T, psm=T, global=T, ...)*

**Matchingcovs** --- Your initial dataframe that includes your entire sample. Must be a data.frame type object that can be coerced to a data.table.

**Treatvar** --- A character vector of the name of your treatment variable which must be included in Matchingcovs.

**ID_vars** --- A character vector of the name of your ID variable which must be included in Matchingcovs. This is necessary for obtaining matching numbers in subsequent phases.

**convert_ind** --- A logical (T or F) indicating whether user wishes to create indicator variables from categorical variables using a design matrix. This function uses stats::model.matrix to create the initial design matrix. Note this means you do not need to create indicator variables before running this function. It defaults to TRUE.

**ref_vars** --- A character vector containing either the names of the indicator variables present in Matchingcovs that the user wishes to assign as the reference variable during logistic regression or one of the categories present in the categorical variable(s) if convert-ind=TRUE.  For simplicty, this defaults to NULL assuming the user wishes the function to choose the first alphabetical category in each variable as the reference variable.  If this is NULL and convert_ind=FALSE, then the reference variable will assumed to not be present in Matchingcovs and thus left out of the balance analysis.

**missing_check** --- Logical (T or F) indicating whether user wishes to impute missing values on matching covariates (default is FALSE). This follows Rosenbaum's (2010) imputation technique for propensity score matching where numeric variables with missing values are assigned a value just outside the range of the sample and then an additional indicator variable is created with missing values coded as 1.  The missing indicator variable is used as a matching covariate.  For categorical variables, the missing value is replaced with a string variable ("missing"). This is then simply added to the design matrix which is created using stats::model.matrix and will be included as one of the indicator variables.

**family** --- Refers to the type of logistic regression to use and defaults to binomial (i.e., logits), which is most common. See glm function in R for other specifications.

**plot_points** --- Logical (T or F) if points should  be added to the denisty plot graph that is displayed in the output.

**pretty_print** --- This logical option (T or F) prints the descriptive statistics and standarized mean differences to an excel sheet that is nicely formatted for convenience.

**add_new_file** --- This logical option (T or F) allows the user to overwrite an existing file by the same name if TRUE or create a new file if user wishes to save earlier output. This defaults to TRUE as normally user will run the diagnostic function once but could be helpful if user wishes to specify multiple logistic regressions to compare predicted PSM values.

**psm** --- This indicates whether the user wishes the diagnostic function to calculate the predicted propensity score. Initially this defaults to a logical (T) but possible values could be F or could be a column name in Matchingcovs if the predicted propensity score was calculated prior to using this function.

**global** --- This is a logical (T or F) specifiying if the user wishes the results of this function to be added as variables in the global environment. The returned results are stored in a list which can also be accessed. This defaults to TRUE.

##Phase 2: PSM Functions

*There are two different types of functions here that end up with the same output structure despite using different packages.*

###1) Nearest_Match_fn

####Default call:
*Nearest_Match_fn<-function(Matchingcovs_PSM, Treatvar, ID_vars, include.psm=F, distance="logit", custom_dist=F, matching_model="glm", family="binomial", weights=NULL, ratio=1, m.order="random", caliper=NULL, exact_vars=NULL, ties=T, distance.tolerance=1e-05, global=T, ...)*

**Matchingcovs_PSM** --- A data.frame object that includes the propensity score. If Diagnostic function with global=T was used, then this will be the dataframe created that is the name of the dataframe you entered in the previous function with _PSM at the end.

**Treatvar, ID_vars** --- Same as those used in PSM_diagnostics function.

**include.psm** --- This is a logical (T or F) indicating whether the user would like to include the psm variable in the logistic regression that the MatchIt package automatically conducts when using the nearest neighbor algorithm. This defaults to FALSE but this option provides flexibility.

**distance** --- This is for compatibility with the MatchIt function but defaults to logit as used in the diagnostic function. Options here include "mahalanobis",binomial link functions if matching model is "glm", generalized additive model link functions is matching model ="gam","nnet" for neural network model, and "rpart" for classification trees. See MatchIt package for more details.

**custom_distance** --- Allows user to enter a pre-specified distance for each observation if calculated outside these functions. Note it defaults to FALSE and simply allows for more flexibility beyond the MatchIt options.

**matching_model** --- Provides complete compatibility with MatchIt and can be either "glm" (generalized linear model) which is the default for basic logistic regression but can also be any model conducted outside of this function that is compatiable with MatchIt such as gam or other more complex logistic regression specifications.

**family** --- Relevant if matching_model="glm" and same as the PSM_diagnostic function and defaults to binomial for logits.

**weights** --- Relevant if matching_model="glm" and allows user to enter weights for each observation for weighted regression. Defaults to NULL and most common use would be to use the PSM scores or inverse weights calculated in the Diagnostic function to conduct a weighted logistic regression to generate the distances used for nearest neighbor matching.

**ratio** --- A numeric value that specifies the number of comparison subjects to match to each treatment unit. Defaults to 1:1 matching.

**m.order** --- can be either "random","largest",or "smallest" and represents which treatment units are matched first for nearest neighbor matching without replacement. Defaults to random but if distributions between groups are dissimilar better matches may be obtained if those treatment units with the largest propensity scores are matched first. Calipers may be more appropriate though.

**caliper** --- Maximum distance allowed between matched comparison and treatment units as measured by standard deviations on propensity scores with 0.2 as the one most commonly recommended in the literature. If no matches for a treatment unit within the caliper can be found then that treatment unit will not be matched and will be excluded from the sample.

**exact_vars** --- This is a character vector of names corresponding to those covariates where exact matching should be done. These names must be contained in the Matchingcovs_PSM dataframe. If NULL which is the default, this function will still try to do exact matching on all covariates and return those matches. If this is TRUE then the exact matching algorithm in MatchIt will be used but otherwise exact matching will be done by the Matching package.

**ties** --- This is for compatibility with the Matching package for exact matching and will select all exact matches rather than picking one at random. Defaults to TRUE and has no effect on Nearest Neighbor matching. 

**distance.tolerance** --- This is for compatibility with Matching package and is similar to caliper but this specification is helpful for those wishing to implement radius matching of 0.05 as described in Dehejia & Wahba (2002). The default is set close to 0, but if set to 0.05 this will include all comparison subjects within a 0.05 distance of each treatment group. Note that using this method means that comparison subjects could potentially be matched to multiple treatment units but paired t-tests should not be used in outcome analysis.

**global** --- A logical (T or F) specifiying if the user wishes the results of this function (matched datasets) to be added as variables in the global environment. The returned results are stored in a list which can also be accessed. This defaults to TRUE.


###2) Optimal_match_fn

####Default call:
*Optimal_match_fn<-function(Matchingcovs_PSM, Treatvar, ID_vars, ratio=1, fitted_glm=T, glm_family="binomial", glmweights=NULL, other_fitted_model=NULL, custom_dist=NULL, Match_type=c("both", "pair", "full"), exact_vars=NULL, include_exact_psm=F, caliper=NULL, Caliper_vars=NULL, Caliper_glm=T, dist_method=c("euclidean", "mahalanobis", "rank_mahalanobis", "custom"), standardization.scale=sd, min.controls = 0, max.controls = Inf, omit.fraction=NULL, global=T)*

***Matchingcovs_PSM, Treatvar, ID_vars, ratio, glm_family, glmweights, custom_dist, exact_vars, caliper, global*** --- (see Nearest_Match_fn)

**fitted_glm** --- Logical (T or F) that specifies if the distances used for optimal matching should be based on basic logistic regression model. If false then other_fitted_model should be specified similar to the matching_model option in the Nearest_Match_fn.

**other_fitted_model** --- Allows flexibility for user to include a different model as used in the matching_model option in Nearest_Match_fn.

**Match_type** --- Allows user to specify whether to use optimal pair matching or optimal full  matching or both (default). If "both" is specified and pair match is unattainable do to caliper then only full matching results will be returned.

**include_exact_psm** --- Logical (T or F) that allows user to specify if the logistic regression used to generate euclidean distances of predicted propensity scores are to be based on the variables where exact matching will be conducted.  If TRUE then you are choosing to use those exact matching variables twice in the matching process.

**Caliper_vars** --- This is a character vector of the column names that the specified caliper should be based on. Allows flexibility to assign a caliper to some rather than all matching covariates. This defaults to NULL but allows the user to create a type of subclassification that does not use the stricter form of exact.

**Caliper_glm** --- Logical (T or F) indicating whether a caliper distance created from the Caliper vars should be calculated using the basic logistic regression. Defaults to T.

**dist_method** --- This can be either "euclidean" (default), "mahalanobis", "rank_mahalanobis", or "custom." Euclidean is essentially the logit distance obtained from running the basic logistic regression. The others are provided for compatibility with the optmatch package. If "custom" is specified, then custom_dist cannot be NULL and should be generated outside the function.

**standardization.scale** --- Allows for compatibility with optmatch and is only relevant when dist_method="euclidean". The default, sd creates a scaled propensity score and is recommended by Rosenbaum (2010). The other options include "mad"" which is a standardized score that is less sensitive to outliers but not well-documented in the literature, and "unstandardized" which are the raw propensity scores.

**min.controls,max.controls** --- When full matching is conducted these specify the range of how many comparison subjects should be matched to each treatment.  The default is to match everyone but sometimes this results in having thousands of comparison students matched to one treatment subject while other treatment subjects are matched to one or not at all.

**omit.fraction** --- Used for compatibility with full matching in optmatch and should be a positive fraction less than one. This option allows the user to specify a fraction of the comparison group to leave unmatched. These will be the worst matches.



##Phase 3: Baseline_Balance_fun

###Default call:
*Baseline_equiv_fn<-function(MatchedSets, Treatvar, ID_vars, pre_pooled_sd, Match_num="Match_num", SD_mean_before=NULL, ref_var=NULL, other_psm=NULL, plot.points=T, outlier_dist=3, global=T)*

**MatchedSets** --- Should be the dataframe object created by either the Nearest_Match_fn or Optimal_match_fn. If global=T was specified in the function(s) used in Phase 2 then these will have the appropriate suffix.

**Treatvar, ID_vars** --- Same as in all prior functions.

**pre_pooled_sd** --- Should include the standardized mean differences obtained from the PSM_diagnostics and will have the suffix PooledSD if global=T was specified in earlier call.

**Match_num** --- Character vector referring to the name of the variable that contains the matched units and if matches were obtained using the functions from Phase 2 will default to "Match_num". This is only for flexbility if matching was done outside these functions and columnn name differs from "Match_num".

**SD_mean_before** --- Optional but must be specified if a comparison between before and after standardized mean differences are desired. Note if global=T was specified in the PSM_diagnostics function then it will have the suffix SDmeans.

**ref_var** --- If balance statistics on the reference variables are desired, then this should be the data set of reference variables for the matches. If global=T was specified in diagnostics function then this will have the suffix ref. It will automatically select only those subjects who were matched by their ID variables even though the object created in the diagnostics include all subjects.

**other_psm** --- Only needed if your psm variable was not created from the earlier functions and must refer to the variable name in your MatchedSets.

**plot.points** --- Refers to the same option as in PSM_diagnostics for density plots.

**outlier_dist** --- Used to identify any matches that may be large on the output graphs. It defaults to identifying the matched pairs whose propensity score differences exceed 3 standard deviations of all such differences. Can be any number but too small will make the graphs less reader friendly. The points on the graphs will be labeled by their match number ("Match_num").

**global** --- Same use as in all previous functions.

