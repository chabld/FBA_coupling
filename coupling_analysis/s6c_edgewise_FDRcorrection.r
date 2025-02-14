#This script uses the permuted edge-wise group differences to estimate a p-value and correct them for FDR afterward

library(dplyr); 

#compute p values of all edge-wise differences, across all pairs of groups and SC metrics
for (FBAmeas in c('fd','fdc','log_fc', 'dti'))
{

edge_cpl_summary=readRDS(paste0('output/edge_wise_output/p-value-based/',FBAmeas,'xFC_sigedges.rds'))
cpl_unperm=edge_cpl_summary$cpl_unperm
cpl_allperms=edge_cpl_summary$cpl_allperms

p_values <- matrix(NA, nrow = nrow(cpl_unperm), ncol = ncol(cpl_unperm), dimnames = list(rep("", nrow(cpl_unperm)), rep("", ncol(cpl_unperm))))
                              
for (mod in 1:nrow(cpl_unperm))
{
  #get only permutations from the one pair of groups
  grpdiff=cpl_unperm[mod,]
  #extract name of group pair to identify their permuted values
  gnames=stringr::str_split(row.names(cpl_unperm)[[mod]], 
                            pattern = '_')[[1]][1]
  grpindex=which(grepl(gnames, rownames(cpl_allperms))==T)
  row.names(p_values)[mod]=gnames

  for (e in 1:ncol(cpl_unperm)) 
  {
    # Calculate the one-sided p-value for each edge
    if (!is.na(grpdiff[e])) {
      #edges either above or under, not tested in absolute differences
      if (grpdiff[e] > 0) 
      { extreme_edges <- sum(as.numeric(cpl_allperms[grpindex,e]) >= as.numeric(grpdiff[e]), na.rm = T) } 
      else 
      { extreme_edges <- sum(as.numeric(cpl_allperms[grpindex,e]) <= as.numeric(grpdiff[e]), na.rm = T) }  
      
      p_values[mod, e] <- extreme_edges / nrow(cpl_allperms[grpindex,])
    }
  } 
}

row.names(p_values)=paste0(row.names(p_values),'_',FBAmeas)
if (!exists('all_ps')){all_ps=p_values} else {all_ps=rbind(all_ps, p_values)}
}

#save all p values
saveRDS(all_ps, './output/edge_wise_output/all_ps.rds')
# Adjust p-values for all ps
adjusted_vector <- p.adjust(as.vector(all_ps), method = "BH")
# Reshape back into the original matrix dimensions
all_ps_adj <- matrix(adjusted_vector, nrow=nrow(all_ps), ncol = ncol(all_ps), dimnames = list(rownames(all_ps), colnames(all_ps)))

#########################report sig edges########################

for (FBAmeas in c('fd','fdc','log_fc', 'dti'))
{
#get original coupling r values
edge_cpl_summary=readRDS(paste0('output/edge_wise_output/p-value-based/',FBAmeas,'xFC_sigedges.rds'))
group_edge_rvals=edge_cpl_summary$group_edge_rvals
cpl_unperm=edge_cpl_summary$cpl_unperm

selected_rows = all_ps_adj[grep(paste0(FBAmeas, "\\b"), rownames(all_ps_adj)), ]

for (mod in 1:nrow(selected_rows))
{
#print results, distinguishing group names/directions
gnames_split=stringr::str_split(row.names(selected_rows)[[mod]], 
                                pattern = '_|v')[[1]][1:2]
print(paste('For', gnames_split[1], 'compared to', gnames_split[2], ',', length(which(selected_rows[mod,]< 0.05)), FBAmeas, 'and FC edges remained significant after FDR.'))

#To quickly get the names of the significant edges
source('jobs/statistics/edge_identifier.r')
sigedges=matrix(NA, nrow=1, ncol=ncol(selected_rows))
sigedges[which(selected_rows[mod,]< 0.05)]=as.numeric(cpl_unperm[mod,which(selected_rows[mod,]< 0.05)])
row.names(sigedges)=stringr::str_split(row.names(selected_rows)[[mod]],pattern = '_')[[1]][1]
assign(paste0(row.names(sigedges),'_sigedges'), sigedges)
#needs the row name to be kept, hence the  drop = FALSE
edge_names=edge_identify(sigedges[1, ,drop = FALSE],
                               group_edge_rvals)
print(edge_names)
}

#merging sigedges for each comparison
sigedges=rbind(ADvCN_sigedges,ADvMCI_sigedges,MCIvCN_sigedges)
row.names(sigedges)=lapply(stringr::str_split(row.names(cpl_unperm), pattern = '_'), function(x) x[1])

edge_cpl_summary$sigedges=sigedges
saveRDS(edge_cpl_summary, paste0('output/edge_wise_output/p-value-based/',FBAmeas,'xFC_sigedges_FDR.rds'))
}

