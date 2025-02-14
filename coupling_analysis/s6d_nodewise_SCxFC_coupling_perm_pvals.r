#This script tests "node coupling" between fixel-weighted SC matrices and FC matrices, correlating SC and FC nodal strength values, for each node, across participants of each group separately


library(igraph);
library(dplyr); 

#node labels
labels_schaefer119 = read.csv("labels_schaefer119.csv")

for (FBAmeas in c('fd','log_fc','fdc','dti'))
{

#load cohortwise matrices [named 2D here because they will be put back into the matrix format]
SCmat_2d=read.csv(paste0('matrices/cohort_matrices/cpl_SCmatlist_',FBAmeas,'.csv'), row.names = 1)
FCmat_2d=read.csv('matrices/cohort_matrices/cpl_FCmatlist.csv', row.names = 1);

#load behavioural data
beh_data=read.csv('ADNI_behdata_cleaned.csv')
beh_data = beh_data[which(beh_data$Subject.ID %in% rownames(FCmat_2d)),]
#binarize sex
beh_data$Sex = dplyr::recode(beh_data$Sex, 'F' = 1, 'M' = 0) 
ICVlist = read.table('ICV.txt', row.names = 1)

#removing ADNI's subjective memory complaint (SMC) group
SCmat_2d=SCmat_2d[which(beh_data$Diag!='SMC'),]
FCmat_2d=FCmat_2d[which(beh_data$Diag!='SMC'),]
beh_data=beh_data[which(beh_data$Diag!='SMC'),]

#########################################################################
#CONTROLLING FOR STREAMLINE COUNT

#This counts the streamlines that were produced with tck2connectome

#tl,dr: This matters because average FD/FDC/FC says nothing about the actual "shape" of the tract. It could be spurious and ridiculously small. 1 streamline can have the same average FD as 1000, yet have a different effect.

#Dhollander:
#"So for example, if you mask had just 4 fixels, where 1 sits in a voxel alone and the 3 other ones sit together in a single voxel, the result of mrstats -output mean is the mean of the FD of all 4 fixels, equally weighted. There’s no “separate” voxel-wise averaging that happens first. This might be important depending on what you actually want to compute, especially for FD very specifically (because it has some properties of a density, even though it’s far from a proper density)! If your initial ROI(s) is/are just WM “blobs” without a specific tract hypothesis, then in lots of scenarios you may not want this behaviour" [he laters talks about voxel-ROI mask averaging which can't apply for connectome2tck protocols as it only has tck, not voxels]
#https://community.mrtrix.org/t/calculating-average-fba-metrics-of-specific-tracts/1805/15

if (FBAmeas!='dti') #not needed for non-fixel based SC
{
  
message("\nWeighting fixel edges to their streamline count...\n")

#turns streamline count (from tck2connectome) matrix to edge vector
ADNI_connectome=read.csv('matrices/cohort_matrices/ADNI_connectome.csv', header = F)
SC_Nstreamline=ADNI_connectome[upper.tri(ADNI_connectome, diag = FALSE)]

SCmatorig=SCmat_2d
for (e in 1:nrow(SCmat_2d)) {SCmat_2d[e,]=SCmatorig[e,]*SC_Nstreamline}

}
###########################################################################
#################Computing weighted node negree (strength) in each modality


message("\nComputing weighted node negree (strength) in SC and FC...\n")

#table where degree vals will be saved cohort wise
SCnodedegree=matrix(NA, nrow=nrow(SCmat_2d), ncol=119,  
              dimnames = list(row.names(SCmat_2d),NULL))
FCnodedegree=SCnodedegree

for (subj in 1:nrow(SCnodedegree)) #for each subject
{
  #Turning each subject level edges back to matrix format
  
  #Get matrix out of individual data
  FCmat = matrix(0, nrow = 119, ncol = 119)
  FCmat[upper.tri(FCmat, diag = F)] = as.numeric(FCmat_2d[subj,1:7021])
  #mirror lower triangle
  tm <- t(FCmat)[,-nrow(FCmat)] 
  FCmat[lower.tri(FCmat)] <- tm[lower.tri(tm, diag = F)]
  FCmat = as.data.frame(FCmat)
  colnames(FCmat) = labels_schaefer119$oldlabels
  row.names(FCmat) = labels_schaefer119$oldlabels
  
  #Get matrix out of individual data
  SCmat = matrix(0, nrow = 119, ncol = 119)
  SCmat[upper.tri(SCmat, diag = F)] = as.numeric(SCmat_2d[subj,1:7021])
  #mirror lower triangle
  tm <- t(SCmat)[,-nrow(SCmat)] 
  SCmat[lower.tri(SCmat)] <- tm[lower.tri(tm, diag = F)]
  SCmat = as.data.frame(SCmat)
  colnames(SCmat) = labels_schaefer119$oldlabels
  row.names(SCmat) = labels_schaefer119$oldlabels
  SCmat[is.na(SCmat)] <- 0
  
  #create graph objects
  SCgraph = graph_from_adjacency_matrix(as.matrix(SCmat), 
                                        mode = "undirected",
                                        weighted = T, diag = F)
  FCgraph = graph_from_adjacency_matrix(as.matrix(FCmat), 
                                        mode = "undirected",
                                        weighted = T, diag = F)
  
  #weighted node degree (strength)
  SCdegree=strength(SCgraph, vids = V(SCgraph), loops = TRUE)
  FCdegree=strength(FCgraph, vids = V(FCgraph), loops = TRUE)
  
  for (node in 1:length(SCdegree))
  { 
   SCnodedegree[subj,node]=SCdegree[node]
   FCnodedegree[subj,node]=FCdegree[node]
  }

}


################################################################
################Node degree coupling in each group#############


#unpermuted observations
message("\nCoupling nodes across groups...\n")

#prepare node coupling degree tables
group_node_rvals = data.frame(matrix(NA, nrow = 3, ncol = ncol(SCmat), dimnames = list(c("AD","MCI","CN"),labels_schaefer119$oldlabels)))

#in each group separately
for (diag in c('AD','MCI','CN')) 
{
  #Selecting group of interest
  SCdegree=SCnodedegree[which(beh_data$Diag == diag),]
  FCdegree=FCnodedegree[which(beh_data$Diag == diag),]
  
  #Get one r FC node degree x SC node degree for each node (across diag subjects)
  for (n in 1:ncol(SCdegree))
  {
   cormodel=cor.test(SCdegree[,n], FCdegree[,n])
   group_node_rvals[diag,n]=cormodel$estimate
  }
}

#unpermutted group differences
cpl_unperm=rbind(
  ADvCN_unperm=group_node_rvals['AD',]-group_node_rvals['CN',],
  ADvMCI_unperm=group_node_rvals['AD',]-group_node_rvals['MCI',],
  MCIvCN_unperm=group_node_rvals['MCI',]-group_node_rvals['CN',])

################################################################
##########PERMUTING Node degree coupling in each group########
  
#permuted differences
nperm=5000; 
message(paste0("\nEstimating node coupling with permuted diagnoses [",nperm, "] ...\n"))


#prepare multiple threads
nthread=6;
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
  group_node_rvals_perm = data.frame(matrix(NA, nrow = 3, ncol = ncol(SCmat), dimnames = list(c("AD","MCI","CN"),  labels_schaefer119$oldlabels)))
    
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
        
        #Selecting group of interest
        SCdegree=SCnodedegree[groupidx,]
        FCdegree=FCnodedegree[groupidx,]
        
        #Get one r FC node degree x SC node degree for each node (across diag subjects)
        for (n in 1:ncol(SCdegree))
        {
          cormodel=cor.test(SCdegree[,n], FCdegree[,n])
          group_node_rvals_perm[diag,n]=cormodel$estimate
        }
       } 

      #calculate coupling difference 
      groupdiff=as.numeric(group_node_rvals_perm[gnames[1],])-as.numeric(group_node_rvals_perm[gnames[2],]);
      #add meaningful name for cpl_perm 
      assign(paste0(comparisons,'_perm'), groupdiff);
          
    
      
    }
  
    #append comparison results based on permutations
    cpl_perm=rbind(ADvCN_perm,ADvMCI_perm,MCIvCN_perm);
    return(cpl_perm)
}  
  
  
  #end bar
  end=Sys.time()
  message(paste("\nCompleted in ",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))
  
  
  
  
###############################################################
#######################Extract significant results#############

  #threshold nodes from the observed differences to extract significant nodes
  for (mod in 1:nrow(cpl_unperm))
  {
    #get only permutations from the said pair of groups
    grpdiff=cpl_unperm[mod,]
    #extract name of group pair to identify their permuted values
    gnames=stringr::str_split(row.names(cpl_unperm)[[mod]], 
                              pattern = '_')[[1]][1]
    grpindex=which(grepl(gnames, rownames(cpl_allperms))==T)
    
    #prepare index of nodes that pass the p<.0025 threshold
    signodes=matrix(NA, nrow=1, ncol=ncol(cpl_unperm))
    for (e in 1:ncol(cpl_unperm))
    {
      # Calculate the one-sided p-value for each node
      if (!is.na(grpdiff[e])) {
        #nodes either above or under, not tested in absolute differences
        if (grpdiff[e] > 0) 
        { extreme_nodes <- sum(as.numeric(cpl_allperms[grpindex,e]) >= as.numeric(grpdiff[e]), na.rm = T) } 
        else 
        { extreme_nodes <- sum(as.numeric(cpl_allperms[grpindex,e]) <= as.numeric(grpdiff[e]), na.rm = T) }  
        
        p_value <- extreme_nodes / nrow(cpl_allperms[grpindex,])
        if (p_value < (0.01/4))
        {signodes[e]=as.numeric(grpdiff[e])}
      }
      
    }
    
    #append node label
    colnames(signodes)=labels_schaefer119$oldlabels
    
    #print results, distinguishing group names/directions
    gnames_split=stringr::str_split(row.names(cpl_unperm)[[mod]], 
                                    pattern = '_|v')[[1]][1:2]
    print(paste0('[',FBAmeas,']'))
    print(paste('For', gnames_split[1], 'compared to', gnames_split[2], ':', length(which(signodes>0)), 'nodes had significantly increased coupling.'))
    
    if (length(which(signodes>0)) != 0) 
    {print(colnames(signodes)[which(signodes>0)])}
    
    cat(paste(length(which(signodes<0)), 'nodes had significantly decreased coupling.'))
    
    if(length(which(signodes<0)) != 0) 
      {print(colnames(signodes)[which(signodes<0)])}
   
    assign(paste0(gnames,'_signodes'), signodes)
  }
  
  
signodes=rbind(ADvCN_signodes,ADvMCI_signodes,MCIvCN_signodes)
row.names(signodes)=lapply(stringr::str_split(row.names(cpl_unperm), pattern = '_'), function(x) x[1])


#save general outcome
saveRDS(list(
  group_node_rvals=group_node_rvals, 
  cpl_allperms=cpl_allperms,
  cpl_unperm=cpl_unperm, 
  signodes=signodes), 
  file=paste0('output/node_wise_output/p-value-based/',FBAmeas,'xFC_signodes.rds'))

}