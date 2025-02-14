#script dedicated to a function that identifies and labels the edges from the index of significant edges, and the associated coupling r value per group, provided in the s6c edge-wise comparison script

#sigedges is a vector containing all edges, 0 if NULL, and a difference value if significant
#group_edge_rvals a vector containing all edges, each with a r value, being the group-wise coupling value as per s6c

edge_identify = function(sigedges, group_edge_rvals)
{
  
edgeidx=which(sigedges!=0)
  
if (length(edgeidx) == 0) 
  { return(print("There are no significant edges.")) }
#get labels of the atlas
labels = read.csv('labels_schaefer119.csv')

#Make dummy matrix
matrix_sc = matrix(0, nrow = 119, ncol = 119)

#add label names to rows and cols
colnames(matrix_sc) = labels$oldlabels
rownames(matrix_sc) = labels$oldlabels

#identify the sig edge's corresponding row and col numbers in the matrix upper triangle
rowcol= which(upper.tri(matrix_sc, diag=F), arr.ind = TRUE)[edgeidx,]

#get group names sigedges
gnames_split=stringr::str_split(row.names(sigedges), 
                                pattern = 'v')[[1]][1:2]

#extract corresponding colnames/rownames to name the connected ROIs
pos_edgename = data.frame(Edges=NA, Group1_r=NA, Group2_r=NA)
colnames(pos_edgename)[2]=gnames_split[1]
colnames(pos_edgename)[3]=gnames_split[2]
neg_edgename = pos_edgename

if (is.matrix(rowcol)==F) #won't be a matrix if only 1 edge
{
    if (sigedges[edgeidx]>0)
    {pos_edgename = data.table::rbindlist(list(
     pos_edgename, data.frame(
     paste0(colnames(matrix_sc)[rowcol[1]], ' x ',
                      row.names(matrix_sc)[rowcol[2]]),    
     group_edge_rvals[gnames_split[1],edgeidx], 
     group_edge_rvals[gnames_split[2],edgeidx])        
    ), use.names=FALSE)
    }
    else if (sigedges[edgeidx]<0) {
     neg_edgename= data.table::rbindlist(list(
     neg_edgename, data.frame(
     paste0(colnames(matrix_sc)[rowcol[1]], ' x ',
             row.names(matrix_sc)[rowcol[2]]),    
     group_edge_rvals[gnames_split[1],edgeidx], 
     group_edge_rvals[gnames_split[2],edgeidx])        
     ), use.names=FALSE)
     }
  
  #remove 1st line
  pos_edgename=pos_edgename[-1,]; 
  neg_edgename=neg_edgename[-1,];
  return(list(pos_edges=pos_edgename,neg_edges=neg_edgename))
  
  
} else 
{
  
    for (e in 1:nrow(rowcol)) 
    { 
      #print negative and positive separately
      if (sigedges[edgeidx[e]]>0)
      { pos_edgename = data.table::rbindlist(list(
            pos_edgename,
            data.frame(Edges=paste0(colnames(matrix_sc)[rowcol[e,1]], ' x ',row.names(matrix_sc)[rowcol[e,2]]),
            group_edge_rvals[gnames_split[1],edgeidx[e]], 
            group_edge_rvals[gnames_split[2],edgeidx[e]])
            ), use.names=FALSE)
      }else if (sigedges[edgeidx[e]]<0)
      {
        neg_edgename = data.table::rbindlist(list(
            neg_edgename,
            data.frame(Edges=paste0(colnames(matrix_sc)[rowcol[e,1]], ' x ',row.names(matrix_sc)[rowcol[e,2]]),
            group_edge_rvals[gnames_split[1],edgeidx[e]], 
            group_edge_rvals[gnames_split[2],edgeidx[e]])
            ), use.names=FALSE)
      }
    
    }
  
  #remove 1st line
  pos_edgename=pos_edgename[-1,]; 
  neg_edgename=neg_edgename[-1,];
  return(list(pos_edges=pos_edgename,neg_edges=neg_edgename))
}

}