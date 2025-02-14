#This script uses the permuted node-wise group differences to estimate a p-value and correct them for FDR afterward

library(dplyr); 
labels_schaefer119 = read.csv("labels_schaefer119.csv")


#compute p values of all node-wise differences, across all pairs of groups and SC metrics
for (FBAmeas in c('fd','fdc','log_fc', 'dti'))
{

node_cpl_summary=readRDS(paste0('output/node_wise_output/p-value-based/',FBAmeas,'xFC_signodes.rds'))
cpl_unperm=node_cpl_summary$cpl_unperm
cpl_allperms=node_cpl_summary$cpl_allperms

p_values <- matrix(NA, nrow = nrow(cpl_unperm), ncol = ncol(cpl_unperm), dimnames = list(rep("", nrow(cpl_unperm)), rep("", ncol(cpl_unperm))))
                              
for (mod in 1:nrow(cpl_unperm))
  {
  #extract name of group pair to identify their permuted values
  gnames=stringr::str_split(row.names(cpl_unperm)[[mod]], 
                            pattern = '_')[[1]][1]
  grpindex=which(grepl(gnames, rownames(cpl_allperms))==T)
  row.names(p_values)[mod]=gnames
  #get only permutations from the said pair of groups
  grpdiff=cpl_unperm[mod,]
  
  # Calculate the 2-sided p-value for each node
  for (node in 1:ncol(cpl_unperm)) {
  if (!is.na(grpdiff[node])) {
    #nodes either above or under, not tested in absolute differences
    if (grpdiff[node] > 0) 
      { extreme_nodes <- sum(as.numeric(cpl_allperms[grpindex,node]) >= as.numeric(grpdiff[node]), na.rm = T) } 
    else 
    { extreme_nodes <- sum(as.numeric(cpl_allperms[grpindex,node]) <= as.numeric(grpdiff[node]), na.rm = T) }  
    
    p_values[mod, node] <- extreme_nodes / nrow(cpl_allperms[grpindex,])
   }
   } 
  }

row.names(p_values)=paste0(row.names(p_values),'_',FBAmeas)
if (!exists('all_ps')){all_ps=p_values} else {all_ps=rbind(all_ps, p_values)}
}

#save all p values
saveRDS(all_ps, './output/node_wise_output/all_ps.rds')
# Adjust p-values for all ps
adjusted_vector <- p.adjust(as.vector(all_ps), method = "BH")
# Reshape back into the original matrix dimensions
all_ps_adj <- matrix(adjusted_vector, nrow=nrow(all_ps), ncol = ncol(all_ps), dimnames = list(rownames(all_ps), colnames(all_ps)))

#######################report sig nodes#########################

for (FBAmeas in c('fd','fdc','log_fc', 'dti'))
{
#get original coupling r values
node_cpl_summary=readRDS(paste0('output/node_wise_output/p-value-based/',FBAmeas,'xFC_signodes.rds'))
group_node_rvals=node_cpl_summary$group_node_rvals
cpl_unperm=node_cpl_summary$cpl_unperm

selected_rows = all_ps_adj[grep(paste0(FBAmeas, "\\b"), rownames(all_ps_adj)), ]

for (mod in 1:nrow(selected_rows))
{

  #print results, distinguishing group names/directions
  gnames=stringr::str_split(row.names(cpl_unperm)[[mod]], 
                            pattern = '_')[[1]][1]
  gnames_split=stringr::str_split(row.names(cpl_unperm)[[mod]], 
                                  pattern = '_|v')[[1]][1:2]
  print(paste('For', gnames_split[1], 'compared to', gnames_split[2], ',', length(which(selected_rows[mod,]< 0.05)), FBAmeas, 'and FC node(s) remained significant after FDR.'))
  
  if (length(which(selected_rows[mod,]< 0.05))!=0) 
  {
  for (sign in which(selected_rows[mod,]< 0.05))
  {
    if (cpl_unperm[mod,sign] > 0)
    {print(paste(colnames(cpl_unperm[sign]),'node(s) had significantly increased coupling.'))} else 
    {print(paste(colnames(cpl_unperm[sign]),'node(s) had significantly decreased coupling.'))}
  }}
  
  #save signodes in a matrix
  signodes=matrix(NA, nrow=1, ncol=ncol(selected_rows))
  colnames(signodes)=labels_schaefer119$oldlabels
  signodes[which(selected_rows[mod,]< 0.05)]=as.numeric(cpl_unperm[mod,which(selected_rows[mod,]< 0.05)])
  row.names(signodes)=stringr::str_split(row.names(selected_rows)[[mod]],pattern = '_')[[1]][1]
  
  assign(paste0(gnames,'_signodes'), signodes)
}

#merging signodes for each comparison
signodes=rbind(ADvCN_signodes,ADvMCI_signodes,MCIvCN_signodes)
row.names(signodes)=lapply(stringr::str_split(row.names(cpl_unperm), pattern = '_'), function(x) x[1])

node_cpl_summary$signodes=signodes
saveRDS(node_cpl_summary, paste0('output/node_wise_output/p-value-based/',FBAmeas,'xFC_signodes_FDR.rds'))
}

