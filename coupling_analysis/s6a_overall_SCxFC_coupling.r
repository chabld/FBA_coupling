#This script tests the correlation between overall SC and  FC weights (correlation between all matching edges (119x119 atlas)), and at subgroup-level (correlating group-wise averaged edge-weights). 

#load cohortwise matrices for each corresponding SC matrix
#these are matrices with N rows per participants * E weights per edges (corresponding to the 119*119 regions, without the diagonal of self-connected regions and reverse directions)
for (FBAmeas in c('fd','log_fc','fdc','dti'))
{
SCmat=read.csv(paste0('matrices/cohort_matrices/cpl_SCmatlist_',FBAmeas,'.csv'), row.names = 1)
FCmat=read.csv('matrices/cohort_matrices/cpl_FCmatlist.csv', row.names = 1);

#remove FC edges that do not exist in the structural connectome
FCmat[is.na(SCmat)] <- NA

#load behavioural data
beh_data=read.csv('ADNI_behdata_cleaned.csv') 
beh_data = beh_data[which(beh_data$Subject.ID %in% rownames(FCmat)),]
beh_data$Sex = dplyr::recode(beh_data$Sex, 'F' = 1, 'M' = 0) #binarize sex
ICVlist = read.table('ICV.txt', row.names = 1)

#removing ADNI's subjective memory complaint (SMC) group
SCmat=SCmat[which(beh_data$Diag!='SMC'),]
FCmat=FCmat[which(beh_data$Diag!='SMC'),]
beh_data=beh_data[which(beh_data$Diag!='SMC'),]

#removing cerebellum edges (track tips may still be connected to ROI despite cerebellum manually brushed off from the fixel template)
all_edges_labels=readRDS('edge_labels_schaefer119.rds')
SCmat[,grep(pattern='CEREBELLUM', all_edges_labels)] <- NA
FCmat[,grep(pattern='CEREBELLUM', all_edges_labels)] <- NA

########################################################################
#Weighting fixel averages with streamline count
#This counts the streamlines that were produced with tck2connectome

#tl,dr: This matters because average FD/FDC/FC says nothing about the actual "shape" of the tract. It could be spurious and ridiculously small. 1 streamline can have the same average FD as 1000, yet have a different importance

#Dhollander:
#"So for example, if you mask had just 4 fixels, where 1 sits in a voxel alone and the 3 other ones sit together in a single voxel, the result of mrstats -output mean is the mean of the FD of all 4 fixels, equally weighted. There’s no “separate” voxel-wise averaging that happens first. This might be important depending on what you actually want to compute, especially for FD very specifically (because it has some properties of a density, even though it’s far from a proper density)! If your initial ROI(s) is/are just WM “blobs” without a specific tract hypothesis, then in lots of scenarios you may not want this behaviour" [he laters talks about voxel-ROI mask averaging which can't apply for connectome2tck protocols as it only has tck, not voxels]
#https://community.mrtrix.org/t/calculating-average-fba-metrics-of-specific-tracts/1805/15

if (FBAmeas!='dti') #not needed for non-fixel based SC
{
  message("\nWeighting fixel edges to their streamline count...\n")
  
  #turns streamline matrix to edge vector
  ADNI_connectome=read.csv('matrices/cohort_matrices/ADNI_connectome.csv', header = F)
  SC_Nstreamline=ADNI_connectome[upper.tri(ADNI_connectome, diag = FALSE)]
  SCmatorig=SCmat
  for (e in 1:nrow(SCmat)) {SCmat[e,]=SCmatorig[e,]*SC_Nstreamline}
}

########################################################################
#overall coupling as correlation between all matching edges of SC and FC per subject

all_r=matrix(NA,nrow=nrow(SCmat))
row.names(all_r)=beh_data$Subject.ID
for (v in 1:nrow(SCmat))
{
  #correlation test needs at least 3 observations
  if(length(SCmat[v,!is.na(SCmat[v,])]) >= 3 & 
     length(FCmat[v,!is.na(FCmat[v,])]) >= 3)
 {  
    vector_cor=cor.test(as.numeric(SCmat[v,]), as.numeric(FCmat[v,]),na.rm=TRUE)
    coupling_val=round(vector_cor$estimate,2)
    all_r[v,]=coupling_val
  }
}

print(paste0('Across the  whole dataset, overall connectome edges, SC (', FBAmeas,') and FC correlate at a range from r = ', min(all_r, na.rm = T),', to r = ', max(all_r, na.rm = T)))

#ANOVA to test effect of diagnosis group on individual SC-FC correlation
dataset= data.frame(all_r=all_r, Diag=beh_data$Diag, 
                    Sex=beh_data$Sex, Age=beh_data$Age)
aov_diag = aov(all_r ~ Diag + Sex + Age, data = dataset)
print(summary(aov_diag))

#if significant effect of group, marginal means contrasts post hoc test
if (summary(aov_diag)[[1]][1,5] < .05)
{
  library(emmeans)
  emm=emmeans(aov_diag, ~ Diag);
  con=pairs(emm, adjust='none');
  print(con)
}

#save those for later plotting
assign(paste0('all_r_',FBAmeas), all_r)

########################################################################
##############################group average and comparisons#############
#overall coupling as correlation between all matching edges of SC and FC, after the SC/FC weights were averaged across group

#To store whole models for later plotting
overall_group_r=matrix(nrow=1, ncol=3)
colnames(overall_group_r)=c('overallSC','overallFC','diag')

group_r=list()
for (diag in c('AD','MCI','CN')) 
{
  overallSC=colMeans(SCmat[which(beh_data$Diag==diag),],na.rm=TRUE)
  overallFC=colMeans(FCmat[which(beh_data$Diag==diag),],na.rm=TRUE) 
  vector_cor=cor.test(overallSC,overallFC)
  coupling_val=round(vector_cor$estimate,3)
  group_r[[diag]]=coupling_val
  
  print(paste0('In ', diag, ' overall (edge-wise SC/FC averaged across groups), SC (', FBAmeas,') and FC correlate at r = ', coupling_val,', at p = ', format(vector_cor$p.value, scientific=F, digits=4)))
  
 
  overall_group_r=rbind(overall_group_r,data.frame(overallSC,overallFC,diag))
  
}
overall_group_r=overall_group_r[-1,]
assign(paste0(FBAmeas, '_overall_group_r'), overall_group_r)  

#test difference between average group coupling using Fisher-Z transform
results1=cocor::cocor.indep.groups(r1.jk=as.numeric(group_r[[1]]),
  n1=sum(beh_data$Diag=='AD'), r2.hm=as.numeric(group_r[[2]]),
  n2=sum(beh_data$Diag=='MCI'));
results2=cocor::cocor.indep.groups(r1.jk=as.numeric(group_r[[1]]),
  n1=sum(beh_data$Diag=='AD'), r2.hm=as.numeric(group_r[[3]]),
  n2=sum(beh_data$Diag=='CN'));
results3=cocor::cocor.indep.groups(r1.jk=as.numeric(group_r[[2]]),
 n1=sum(beh_data$Diag=='MCI'), r2.hm=as.numeric(group_r[[3]]),
 n2=sum(beh_data$Diag=='CN'));
print(paste0('AD vs MCI, AD vs CN, MCI vs CN:'))
print(c(results1@fisher1925$p.value, results2@fisher1925$p.value, results3@fisher1925$p.value))
}

overallmatrix=cbind(all_r_fd, all_r_log_fc, all_r_fdc, all_r_dti)
colnames(overallmatrix)=c('all_r_fd', 'all_r_log_fc', 'all_r_fdc', 'all_r_dti')
saveRDS(overallmatrix, file='output/overall_output/Whole_connectome_overall_coupling_all_r.rds')