#This script tests the correlation between fixel-weighted SC and FC at each of the 7021 edges (119x119 atlas) across participants of each subgroup. 

library(dplyr); 

#node labels
labels_schaefer119 = read.csv("labels_schaefer119.csv")

for (FBAmeas in c('fd','log_fc','fdc','dti'))
{
  
#load cohortwise matrices
SCmat=read.csv(paste0('matrices/cohort_matrices/cpl_SCmatlist_',FBAmeas,'.csv'), row.names = 1)
FCmat=read.csv('matrices/cohort_matrices/cpl_FCmatlist.csv', row.names = 1);

#remove edges that do not exist in the structural connectome
FCmat[is.na(SCmat)] <- NA

#load behavioural data
beh_data=read.csv('ADNI_behdata_cleaned.csv')
beh_data = beh_data[which(beh_data$Subject.ID %in% rownames(FCmat)),]
#binarize sex
beh_data$Sex = dplyr::recode(beh_data$Sex, 'F' = 1, 'M' = 0) 
ICVlist = read.table('ICV.txt', row.names = 1)

#removing ADNI's subjective memory complaint (SMC) group
SCmat=SCmat[which(beh_data$Diag!='SMC'),]
FCmat=FCmat[which(beh_data$Diag!='SMC'),]
beh_data=beh_data[which(beh_data$Diag!='SMC'),]

#removing cerebellum edges (track tips may still be connected to ROI despite cerebellum manually brushed off from the fixel template)
all_edges_labels=readRDS('edge_labels_schaefer119.rds')
SCmat[,grep(pattern='CEREBELLUM', all_edges_labels)] <- NA
FCmat[,grep(pattern='CEREBELLUM', all_edges_labels)] <- NA

#########################################################################Weighting fixel averages with streamline count

if (FBAmeas!='dti')  #not needed for non-fixel based SC
{
  
message("\nWeighting fixel edges to their streamline count...\n")

#turns streamline count (from tck2connectome) matrix to edge vector
ADNI_connectome=read.csv('matrices/cohort_matrices/ADNI_connectome.csv', header = F)
SC_Nstreamline=ADNI_connectome[upper.tri(ADNI_connectome, diag = FALSE)]

SCmatorig=SCmat
for (e in 1:nrow(SCmat)) {SCmat[e,]=SCmatorig[e,]*SC_Nstreamline}

}

#######################Unpermuted group-wise coupling#################

message("\nEstimating coupling with observed diagnoses...\n")

group_edge_rvals = data.frame(matrix(NA, nrow = 3, ncol = ncol(SCmat), dimnames = list(c("AD","MCI","CN"))))

for (diag in c('AD','MCI','CN')) {
  
  SCxFC_coupling=matrix(NA, nrow=2, ncol=ncol(SCmat))
  
  #identify subj from selected group
  groupidx = which(beh_data$Diag==diag)
  
  for (v in 1:ncol(SCmat))
  { #only existing edges
    if  (all(is.na(c(NA, NaN, SCmat[groupidx,v])))==F & 
         all(is.na(c(NA, NaN, FCmat[groupidx,v])))==F & 
         #at least 3 subj 
         length(which(!is.na(FCmat[groupidx,v]))) >= 3) 
    {
      cormodel=cor.test(as.numeric(SCmat[groupidx,v]), 
                        as.numeric(FCmat[groupidx,v]))
      SCxFC_coupling[1,v]=cormodel$estimate
    }
  }  
  
  group_edge_rvals[diag,]=SCxFC_coupling[1,]
}

group_edge_rvals = data.frame(matrix(NA, nrow = 3, ncol = ncol(SCmat), dimnames = list(c("AD","MCI","CN"))))

for (diag in c('AD','MCI','CN')) {
  
  SCxFC_coupling=matrix(NA, nrow=2, ncol=ncol(SCmat))
  
  #identify subj from selected group
  groupidx = which(beh_data$Diag==diag)
  
  for (v in 1:ncol(SCmat))
  { #only existing edges
    if  (all(is.na(c(NA, NaN, SCmat[groupidx,v])))==F & 
         all(is.na(c(NA, NaN, FCmat[groupidx,v])))==F & 
         #at least 3 subj 
         length(which(!is.na(FCmat[groupidx,v]))) >= 3) 
    {
      cormodel=cor.test(as.numeric(SCmat[groupidx,v]), 
                        as.numeric(FCmat[groupidx,v]))
      SCxFC_coupling[1,v]=cormodel$estimate
    }
  }  
  
  group_edge_rvals[diag,]=SCxFC_coupling[1,]
}


#unpermutted group differences
cpl_unperm=rbind(
ADvCN_unperm=group_edge_rvals['AD',]-group_edge_rvals['CN',],
ADvMCI_unperm=group_edge_rvals['AD',]-group_edge_rvals['MCI',],
MCIvCN_unperm=group_edge_rvals['MCI',]-group_edge_rvals['CN',])

#################Shuffling groups before coupling#######################

#permuted differences
nperm=5000; 
message(paste0("\nEstimating coupling with permuted diagnoses [",nperm, "] ...\n"))

#prepare multiple threads
nthread=5;
cl=parallel::makeCluster(nthread);
doParallel::registerDoParallel(cl);
`%dopar%` = foreach::`%dopar%`;

#progress bar
doSNOW::registerDoSNOW(cl);
pb=txtProgressBar(max = nperm, style = 3);
progress=function(n) setTxtProgressBar(pb, n);
opts=list(progress = progress);

#activate parallel processing
unregister_dopar = function() { 
.foreachGlobals <- utils::getFromNamespace(".foreachGlobals", "foreach");
env = .foreachGlobals; rm(list=ls(name=env), pos=env)}

##fitting permuted regression model and extracting t-stats in parallel streams
start=Sys.time()

cpl_allperms=foreach::foreach(perm=1:nperm, .combine= rbind, .options.snow = opts)  %dopar%
{
  
    #prepare data.frame to store group-wise coupling values 
    group_edge_rvals_perm = data.frame(matrix(NA, nrow = 3, ncol = ncol(SCmat), dimnames = list(c("AD","MCI","CN"))));
    
    #get SC-FC correlation across subj of each group separately
    for (comparisons in c('ADvCN','ADvMCI','MCIvCN')) 
    {
      #extract name of the two groups
      gnames=stringr::str_split(comparisons, pattern = 'v')[[1]][1:2];
      #only shuffle across the two groups compared
      diag_perm=beh_data$Diag;
      diag_perm[which(diag_perm == gnames[1]|diag_perm == gnames[2])] <- sample(diag_perm[which(diag_perm == gnames[1]|diag_perm == gnames[2])]);
      
      #compute coupling of the two groups
      for (diag in gnames) {
        #individual group coupling vector
        SCxFC_coupling=matrix(NA, nrow=1, ncol=ncol(SCmat));
        #identify subj from selected permuted group
        groupidx = which(diag_perm==diag);
        #correlate edges across participants from the group
        for (v in 1:ncol(SCmat))
        { #only existing edges
          if  (all(is.na(c(NA, NaN, SCmat[groupidx,v])))==F & 
               all(is.na(c(NA, NaN, FCmat[groupidx,v])))==F &
               #at least 3 subj 
               length(which(!is.na(SCmat[groupidx,v]))) >= 3) 
          {cormodel=cor.test(as.numeric(SCmat[groupidx,v]), 
                             as.numeric(FCmat[groupidx,v]))
            SCxFC_coupling[1,v]=cormodel$estimate}
        } 
        group_edge_rvals_perm[diag,]=SCxFC_coupling[1,];
      }
      
       #calculate coupling difference 
       groupdiff=as.numeric(group_edge_rvals_perm[gnames[1],])-as.numeric(group_edge_rvals_perm[gnames[2],]);
       #add meaningful name for cpl_perm 
       assign(paste0(comparisons,'_perm'), groupdiff);
      
      }    
    
    #append comparison results based on permutation
    cpl_perm=rbind(ADvCN_perm,ADvMCI_perm,MCIvCN_perm);
    return(cpl_perm)
}

#end bar
end=Sys.time()
message(paste("\nCompleted in ",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))





###############################################################
########################SUMMARY AND REPORT#####################


#threshold edges from the observed differences to extract significant edges
for (mod in 1:nrow(cpl_unperm))
{
  #get only permutations from the one pair of groups
  grpdiff=cpl_unperm[mod,]
  #extract name of group pair to identify their permuted values
  gnames=stringr::str_split(row.names(cpl_unperm)[[mod]], 
                            pattern = '_')[[1]][1]
  grpindex=which(grepl(gnames, rownames(cpl_allperms))==T)
  
  #prepare index of edges that pass the p<.0025 threshold
  sigedges=matrix(NA, nrow=1, ncol=ncol(cpl_unperm))
  for (e in 1:ncol(cpl_unperm))
  {
    # Calculate the one-sided p-value for each edge
    if (!is.na(grpdiff[e])) {
      #edges either above or under, not tested in absolute differences
      if (grpdiff[e] > 0) 
      { extreme_edges <- sum(as.numeric(cpl_allperms[grpindex,e]) >= as.numeric(grpdiff[e]), na.rm = T) } 
      else 
      { extreme_edges <- sum(as.numeric(cpl_allperms[grpindex,e]) <= as.numeric(grpdiff[e]), na.rm = T) }  
      
      p_value <- extreme_edges / nrow(cpl_allperms[grpindex,])
      if (p_value < (0.01/4))
      {sigedges[e]=as.numeric(grpdiff[e])}
    }
    
  }
  
  
  #print results, distinguishing group names/directions
  gnames_split=stringr::str_split(row.names(cpl_unperm)[[mod]], 
                            pattern = '_|v')[[1]][1:2]
  print(paste('For', gnames_split[1], 'compared to', gnames_split[2], ',', length(which(sigedges>0)), FBAmeas, 'and FC edges had significantly increased coupling, while', length(which(sigedges<0)), 'significantly decreased.'))
  
  assign(paste0(gnames,'_sigedges'), sigedges)
}

#merging sigedges for each comparison
sigedges=rbind(ADvCN_sigedges,ADvMCI_sigedges,MCIvCN_sigedges)
row.names(sigedges)=lapply(stringr::str_split(row.names(cpl_unperm), pattern = '_'), function(x) x[1])

#To quickly get the names of the significant edges
#needs the row name to be kept, hence the  drop = FALSE
source('jobs/statistics/edge_identifier.r')
ADvCN_edge_names=edge_identify(sigedges[1, ,drop = FALSE],
                               group_edge_rvals)
ADvMCI_edge_names=edge_identify(sigedges[2, ,drop = FALSE],
                                group_edge_rvals)
MCIvCN_edge_names=edge_identify(sigedges[3, ,drop = FALSE],
                                group_edge_rvals)

#save general outcome
saveRDS(list(
  group_edge_rvals=group_edge_rvals, 
  cpl_allperms=cpl_allperms, 
  cpl_unperm=cpl_unperm, 
  sigedges=sigedges, 
  edge_names= list(ADvCN_edge_names=ADvCN_edge_names,
                   ADvMCI_edge_names=ADvMCI_edge_names,
                   MCIvCN_edge_names=MCIvCN_edge_names)), 
  file=paste0('output/edge_wise_output/p-value-based/',FBAmeas,'xFC_sigedges.rds'))

}