#This script is the same as 6a, except it restricts its analysis to within-DMN edges instead of edges of the whole connectome
#The script tests the correlation between overall fixel-weighted SC and overall FC (correlation between edges of the default mode network), and at subgroup-level (correlated group-wise averaged edge-weights. 

stattable=matrix(NA, nrow=3, ncol=8) #for paper
row.names(stattable)=c('AD','MCI','CN')
colnames(stattable)=c('fd.r','fd.p','log_fc.r','log_fc.p','fdc.r','fdc.p','dti.r','dti.p')

#load cohortwise matrices
for (FBAmeas in c('fd','log_fc','fdc','dti'))
{
SCmat=read.csv(paste0('matrices/cohort_matrices/cpl_SCmatlist_',FBAmeas,'.csv'), row.names = 1)
FCmat=read.csv('matrices/cohort_matrices/cpl_FCmatlist.csv', row.names = 1);

#remove edges that do not exist in the structural connectome
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

#removing cerebellum edges (track tips connected to ROI despite cerebellum brushed off the fixel template)
all_edges_labels=readRDS('edge_labels_schaefer119.rds')
SCmat[,grep(pattern='CEREBELLUM', all_edges_labels)] <- NA
FCmat[,grep(pattern='CEREBELLUM', all_edges_labels)] <- NA

########################################################################
######Selecting only within-DMN edges

#node labels
labels_schaefer119 = read.csv("labels_schaefer119.csv")

#Default mode network nodes based on the Yeo classification  
DMN_node_pattern='Default'

#read edge labels to identify within DMN edges
edge_labels_schaefer119 <- readRDS("edge_labels_schaefer119.rds")

#count internal DMN edges (if 1, only connecting 1 DMN node, if 2, two, so within DMN)
library(stringr)
DMN_edges <- str_count(edge_labels_schaefer119,paste(DMN_node_pattern, collapse = "|"))

#remove non-within-DMN edges
SCmat[,which(DMN_edges!=2)] <- NA
FCmat[,which(DMN_edges!=2)] <- NA

#########################################################################Weighting fixel averages with streamline count

if (FBAmeas!='dti') #not needed for non-fixel based SC
{
  
message("\nWeighting fixel edges to their streamline count...\n")

#turns streamline matrix to edge vector
ADNI_connectome=read.csv('matrices/cohort_matrices/ADNI_connectome.csv',
                         header = F)
SC_Nstreamline=ADNI_connectome[upper.tri(ADNI_connectome, diag = FALSE)]

SCmatorig=SCmat
for (e in 1:nrow(SCmat)) {SCmat[e,]=SCmatorig[e,]*SC_Nstreamline}

}

#########################################################################overall coupling as correlation between all (within-DMN) matching edges of SC and FC per subject

all_r=matrix(NA,nrow=nrow(SCmat))
row.names(all_r)=beh_data$Subject.ID
for (v in 1:nrow(SCmat))
{
  #needs at least 3 observations
  if(length(SCmat[v,!is.na(SCmat[v,])]) >= 3 & 
     length(FCmat[v,!is.na(FCmat[v,])]) >= 3)
 {  
    vector_cor=cor.test(as.numeric(SCmat[v,]), as.numeric(FCmat[v,]),na.rm=TRUE)
    coupling_val=round(vector_cor$estimate,2)
    all_r[v,]=coupling_val
  }
}

print(paste0('Across the  whole dataset, overall DMN (edge-wise SC/FC per subject), SC (', FBAmeas,') and FC correlate at a range from r = ', min(all_r, na.rm = T),', to r = ', max(all_r, na.rm = T)))

#test effect of Diag on individual SC-FC correlation
dataset= data.frame(all_r=all_r, Diag=beh_data$Diag, 
                    Sex=beh_data$Sex, Age=beh_data$Age)
aov_diag = aov(all_r ~ Diag + Sex + Age, data = dataset)
print(summary(aov_diag))

#if significant diag effect, marginal means contrasts post hoc test
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
#overall coupling as correlation between all (within-DMN) matching edges of SC and FC, after the SC/FC weights were averaged across group


group_r=list()
for (diag in c('AD','MCI','CN')) 
{
  overallSC=colMeans(SCmat[which(beh_data$Diag==diag),],na.rm=TRUE)
  overallFC=colMeans(FCmat[which(beh_data$Diag==diag),],na.rm=TRUE) 
  vector_cor=cor.test(overallSC,overallFC)
  coupling_val=round(vector_cor$estimate,3)
  group_r[[diag]]=coupling_val
  
  print(paste0('In ', diag, ' overall DMN edges, SC (', FBAmeas,') and FC correlate at r = ', coupling_val,', at p = ', format(vector_cor$p.value, scientific=F, digits=4)))
  
  
  stattable[diag,paste0(FBAmeas,'.r')]=coupling_val
  stattable[diag,paste0(FBAmeas,'.p')]=format(vector_cor$p.value, scientific=F, digits=4)
}

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

overallDMNmatrix=cbind(all_r_fd, all_r_log_fc, all_r_fdc, all_r_dti)
colnames(overallDMNmatrix)=c('all_r_fd', 'all_r_log_fc', 'all_r_fdc', 'all_r_dti')
saveRDS(overallDMNmatrix, 'output/overall_output/DMN_overall_coupling_all_r.rds')
