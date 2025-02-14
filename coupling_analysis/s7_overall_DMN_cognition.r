#This script aims at testing the relevance of overall DMN coupling, for nodes and within-DMN edges, with regards to cognitive measures (memory composite score computed from the RAVLT, and MMSE)

library(dplyr); 
library(igraph);

#preparing concatenated statistics tables
DMN_stattable=list(); DMN_confinttable=list();

#select cognitive measures of interest:
varofinterest_names=c('MEMORY_COMPOSITE','MMSCORE')
cogn=0
for (varofinterest_name in varofinterest_names)
{
  cogn=cogn+1
#######################################################################################Load base data######################
  #node labels
  labels_schaefer119 = read.csv("labels_schaefer119.csv")

  #load FC matrix
  FCmat=read.csv('matrices/cohort_matrices/cpl_FCmatlist.csv', row.names = 1);

  #load behavioural data (ADNI data could not be made available as it requires request for access)
  beh_data=read.csv('ADNI_behdata_cleaned.csv')
  beh_data = beh_data[which(beh_data$Subject.ID %in% rownames(FCmat)),]
  #binarize sex
  beh_data$Sex = dplyr::recode(beh_data$Sex, 'F' = 1, 'M' = 0) 
  #load ICV
  ICVlist = read.table('ICV.txt', row.names = 1)

##########################################################################make MEMORY composite score with CFA

if (varofinterest_name=='MEMORY_COMPOSITE')
{
  library(lavaan)
  MEMO_MODEL <- ' MEMORY_COMPOSITE =~ RAVLT_learning + RAVLT_immediate + AVDEL30MIN + AVDELTOT '
  fit <- suppressWarnings(cfa(model = MEMO_MODEL, data = beh_data[which(beh_data$Diag!='SMC'),],  estimator = "MLR" , missing="fiml")) 
  SumMEMO_MODEL = summary(fit, fit.measures=TRUE,standardized=TRUE)
  
  #normalise it (0 to 1)  
  normalize <- function(x) {(x - min(x,na.rm = T)) / (max(x,na.rm = T) - min(x,na.rm = T))}
  MEMORY_COMPOSITE <- normalize(lavPredict(fit, type = "lv"))
  
  #append to beh_data
  beh_data$MEMORY_COMPOSITE=NA
  beh_data[which(beh_data$Diag!='SMC'),'MEMORY_COMPOSITE']=MEMORY_COMPOSITE
}
  
################################################################
#load SC mat depending on SC metric
for (FBAmeas in c('fd','fdc','log_fc','dti'))
{

  #load cohortwise matrices
  SCmat=read.csv(paste0('matrices/cohort_matrices/cpl_SCmatlist_',FBAmeas,'.csv'), row.names = 1)
  
  #removing ADNI's subjective memory complaint (SMC) group
  SCmat=SCmat[which(beh_data$Diag!='SMC'),]
  FCmat=FCmat[which(beh_data$Diag!='SMC'),]
  beh_data=beh_data[which(beh_data$Diag!='SMC'),]

#########################################################################Weighting fixel averages with streamline count
  
  if (FBAmeas!='dti') #not needed for non-fixel based SC
  {
  message("\nWeighting fixel edges to their streamline count...\n")
  
  #turns streamline matrix to edge vector
  ADNI_connectome=read.csv('matrices/cohort_matrices/ADNI_connectome.csv', header = F)
  SC_Nstreamline=ADNI_connectome[upper.tri(ADNI_connectome, diag = FALSE)]
  
  SCmatorig=SCmat
  for (e in 1:nrow(SCmat)) {SCmat[e,]=SCmatorig[e,]*SC_Nstreamline}
  }
  

################################################################
###################Extract within-DMN edges#####################


  #Default mode network node index based on the Yeo-7 classification
  DMN_nodeidx=c(grep('Default',labels_schaefer119$oldlabels))
    
  #read edge labels to identify within DMN edges
  edge_labels_schaefer119 <- readRDS("edge_labels_schaefer119.rds")
  
  #count specifically DMN edges (if 1, at least connecting 1 DMN node, if 2, two, so within DMN)
  library(stringr)
  DMN_edges <- str_count(edge_labels_schaefer119,paste(labels_schaefer119$oldlabels[DMN_nodeidx], collapse = "|"))
  sigedges=which(DMN_edges==2)
    #For each subject, compute coupling across all these edges
        #individual level
        all_r_edges=matrix(NA,nrow=nrow(SCmat))
        row.names(all_r_edges)=beh_data$Subject.ID
        for (v in 1:nrow(SCmat))
        {
          #needs at least 3 observations
          if(length(sigedges) >= 3
             & length(which(!is.na(SCmat[v,sigedges]))) >= 3)
          {  
            vector_cor=cor.test(
              as.numeric(SCmat[v,sigedges][which(!is.na(SCmat[v,sigedges]))]), 
              as.numeric(FCmat[v,sigedges][which(!is.na(SCmat[v,sigedges]))]))
            coupling_val=vector_cor$estimate
            all_r_edges[v,]=coupling_val
          }
        }
        

################################################################
##Computing weighted node negree (nodal strength) in each SC type for DMN nodes
        
        #get index of DMN nodes
        signodes=matrix(NA, nrow=1, ncol=119)
        colnames(signodes)=labels_schaefer119$oldlabels
        signodes=which(colnames(signodes) %in% labels_schaefer119$oldlabels[DMN_nodeidx])
        
        message("\nComputing weighted node negree (strength) in SC and FC...\n")
        #distinguish for clarity (2d are matrices with N row per participant * E weight per edge, not adjacency matrices)
        FCmat_2d=FCmat;SCmat_2d=SCmat;
        #table where degree vals will be saved cohort wise
        SCnodedegree=matrix(NA, nrow=nrow(SCmat), ncol=119,  
                            dimnames = list(row.names(SCmat),NULL))
        FCnodedegree=SCnodedegree
        
        for (subj in 1:nrow(SCnodedegree)) #for each subject
        {
          #Get FC adjacency matrix out of individual data
          FCmatadj = matrix(0, nrow = 119, ncol = 119)
          FCmatadj[upper.tri(FCmatadj, diag = F)] = as.numeric(FCmat_2d[subj,1:7021])
          #mirror lower triangle
          tm <- t(FCmatadj)[,-nrow(FCmatadj)] 
          FCmatadj[lower.tri(FCmatadj)] <- tm[lower.tri(tm, diag = F)]
          FCmatadj = as.data.frame(FCmatadj)
          colnames(FCmatadj) = labels_schaefer119$oldlabels
          row.names(FCmatadj) = labels_schaefer119$oldlabels
          
          #Get SC adjacency matrix out of individual data
          SCmatadj = matrix(0, nrow = 119, ncol = 119)
          SCmatadj[upper.tri(FCmatadj, diag = F)] = as.numeric(SCmat_2d[subj,1:7021])
          #mirror lower triangle
          tm <- t(SCmatadj)[,-nrow(SCmatadj)] 
          SCmatadj[lower.tri(SCmatadj)] <- tm[lower.tri(tm, diag = F)]
          SCmatadj = as.data.frame(SCmatadj)
          colnames(SCmatadj) = labels_schaefer119$oldlabels
          row.names(SCmatadj) = labels_schaefer119$oldlabels
          SCmatadj[is.na(SCmatadj)] <- 0
          
          #create graph objects
          SCgraph = graph_from_adjacency_matrix(as.matrix(SCmatadj), mode = "undirected", weighted = T, diag = F)
          FCgraph = graph_from_adjacency_matrix(as.matrix(FCmatadj), mode = "undirected", weighted = T, diag = F)
          
          #compute weighted node degree (strength)
          SCdegree=strength(SCgraph, vids = V(SCgraph), loops = TRUE)
          FCdegree=strength(FCgraph, vids = V(FCgraph), loops = TRUE)
          
          for (node in signodes)
          { 
            SCnodedegree[subj,node]=SCdegree[node]
            FCnodedegree[subj,node]=FCdegree[node]
          }
          
        }        
        
        #Get one Pearson r correlating FC node strength to SC node strength for each node (across subjects of a group)
        message("\nCoupling nodes across groups...\n")
        
        #prepare node coupling degree tables
        all_r_nodes=matrix(NA,nrow=nrow(SCmat))
        row.names(all_r_nodes)=beh_data$Subject.ID
        for (n in 1:nrow(SCmat))
        {
          cormodel=cor.test(SCnodedegree[n,signodes], FCnodedegree[n,signodes])
          all_r_nodes[n,]=cormodel$estimate
        } 
          
################################################################
#####Modelling DMN edges/nodes across the whole sample##########
        
        #define response variable and only keep participants with available neuropsychological data
        varofinterest= grep(varofinterest_name, colnames(beh_data));
        beh_data_pred = beh_data[is.na(beh_data[,varofinterest])==FALSE,];
        overallcoupling_edges = all_r_edges[which(rownames(all_r_edges) %in% beh_data_pred$Subject.ID),];
        overallcoupling_nodes=all_r_nodes[which(rownames(all_r_nodes) %in% beh_data_pred$Subject.ID),];
        ICV_pred = ICVlist[which(rownames(ICVlist) %in% beh_data_pred$Subject.ID),]
        
        dataset_pred = cbind(beh_data_pred[,c('Age', 
                                              'Sex', 
                                              'Diag')],
                             overallcoupling_edges,
                             overallcoupling_nodes,
                             ICV_pred);
        row.names(dataset_pred)=beh_data_pred$Subject.ID
        
#-------------------------------RUN MODEL--------------------
        
        stattable=matrix(NA, nrow=2) #will store statistics
        confinttable=matrix(NA, nrow=2) #will stores CIs of adjusted beta values
        
        for (diag2 in c('All','AD','MCI','CN')) 
        {
          
          if (diag2=='All')
          {
            dataset = dataset_pred #whole dataset  
          } else { 
            #only specific group
            dataset = dataset_pred[which(dataset_pred$Diag==diag2),]
          } 
          
          ##define dependent variable
          y = beh_data_pred[which(beh_data_pred$Subject.ID %in% row.names(dataset)),varofinterest];
          #define independent variables
          SCFCcoupling_edges = dataset$overallcoupling_edges;
          SCFCcoupling_nodes = dataset$overallcoupling_nodes;
          age = dataset[,'Age']; sex = dataset[,'Sex']; 
          ICV2 = dataset$ICV_pred
          
          #run linear model
          model <- lm(y ~ SCFCcoupling_edges + SCFCcoupling_nodes + age + sex + ICV2, data=dataset)
          
#---------------------ASSUMPTIONS CHECKS-------------------
          
          library(car)
          #check multicollinearity: 
          #VIF values: below 10 is good, see Myers, R.H. (1990) Classical and Modern Regression Application (2nd edn). Boston: Duxbury Press.
          #"tolerance" value (1/VIF): below 0.02 is a problem, cf Menard, S. (1995). Applied logistic regression analysis. Sage university paper series on quantitative applications in the social sciences, 07-106. Thousand Oaks, CA: Sage.
          #Intercorrelation above 0.9 : Mayers, A. (2013). Introduction to statistics and SPSS in psychology. Edinburgh Gate, Harlow: Pearson Education Limited; Field, A. (2013). Discovering statistics using IBM SPSS statistics.).
          
          if(!all(vif(model) <= 5) | !all(1/vif(model) >= 0.2) |
             cor.test(SCFCcoupling_edges,SCFCcoupling_nodes)$estimate > 0.8) {
            cat('Multicollinearity detected:\n')
            cat(which(vif(model) >= 5))
            cat('\n Correlation between node-wise and edge-wise coupling:\n')
            cat(cor.test(SCFCcoupling_edges,SCFCcoupling_nodes)$estimate)}
          
          #check independence of errors (Durbin-Watson statistics between 1 and 3 is good according to Mayers, A. (2013))
          if (durbinWatsonTest(model)[2] < 1 | durbinWatsonTest(model)[2] > 3) 
            { cat('Residuals are autocorrelated')}
          
          
#---------------------RESULTS SAVING AND COLLATING-------------
          
          
          print(paste('Effects of', FBAmeas,'x FC coupling in DMN on', colnames(beh_data)[varofinterest], 'in', diag2)) 
          
          #appends adjusted beta estimates
          print(cbind(summary(model)[[4]][2:3,1:4], adj.beta=QuantPsyc::lm.beta(model)[1:2]))
          #extracts only adj.beta and p columns:
           library(QuantPsyc)
           stattable=cbind(stattable, as.matrix(cbind(adj.beta=lm.beta(model)[1:2], summary(model)[[4]][2:3,1:4]))[,c(1,5)])
        #save confidence intervals: this needs standard errors because you want the confidence intervals in the forest plots to be based on adjusted beta, not the original beta
       std_errors <- coef(summary(model))[2:3, "Std. Error"]
       sd_x <- apply(na.omit(dataset[,c("overallcoupling_edges",
                                 "overallcoupling_nodes")]), 2, sd)
       sd_y <- sd(y)
       std_errors_beta <- std_errors * (sd_x / sd_y)
       confintbeta <- cbind(Lower = lm.beta(model)[1:2] - 1.96 * std_errors_beta, Upper = lm.beta(model)[1:2] + 1.96 * std_errors_beta)
       confinttable=cbind(confinttable,as.matrix(confintbeta))
       }
        
        assign(paste0('stattable_',FBAmeas), stattable[,-1])
        assign(paste0('confinttable_',FBAmeas), confinttable[,-1])

}

cog_stattable=rbind(fd=stattable_fd, 
                    log_fc=stattable_log_fc, 
                    fdc=stattable_fdc, 
                    dti=stattable_dti)
cog_confinttable=rbind(fd=confinttable_fd, 
                    log_fc=confinttable_log_fc, 
                    fdc=confinttable_fdc, 
                    dti=confinttable_dti)

DMN_stattable=append(DMN_stattable, list(cog_stattable));
DMN_confinttable=append(DMN_confinttable,list(cog_confinttable));
}

################################################################
######################Plotting forest plots ###################
library(forestplot)
#This was programmed to report results for 'All', not individual groups 

#FDR correction of p-values
corrected_p=p.adjust(unlist(lapply(DMN_stattable, function(x) x[, 2])), method='fdr')

#for each cognitive measure, forest plot of its adjusted beta, CI, p
plotlist=list()
for (i in 1:cogn)
{

#will select the p that correspond to the measure's respective table
start_idx <- floor((i - 1) * length(corrected_p) / cogn) + 1 
end_idx <- min(floor(i * length(corrected_p) / cogn))

forestdata <- data.frame(
  beta  = DMN_stattable[[i]][c(1:8),1],
  lower = DMN_confinttable[[i]][c(1:8),1],
  upper = DMN_confinttable[[i]][c(1:8),2],
  level = c("FD–FC edges", "FD–FC nodes", 
            "FbC–FC edges", "FbC–FC nodes", 
            "FDC–FC edges", "FDC–FC nodes", 
            "Streamline SC–FC edges", "Streamline SC–FC nodes"),
  pvalue = corrected_p[start_idx:end_idx])

#check significance after fdr
add_asterisks <- function(p) { if (p < 0.001) { return("***") } else if (p < 0.01) { return("**") } else if (p < 0.05) { return("*") } else { return("") } }
forestdata$asterisks <- sapply(forestdata$pvalue, add_asterisks)

if(varofinterest_names[i]=='MMSCORE')
{title="Effects of coupling on MMSE score"}
else if (varofinterest_names[i]=='MEMORY_COMPOSITE')
{title="Effects of coupling on memory performance"}
else {title=paste( "Effects of coupling on", varofinterest_names[i])}

level_colors = c("FD–FC edges"="#648FFF", "FD–FC nodes"="#648FFF", 
  "FbC–FC edges"= "#785EF0", "FbC–FC nodes"= "#785EF0", 
  "FDC–FC edges"= "#0072B2", "FDC–FC nodes"= "#0072B2", 
  "Streamline SC–FC edges"= "#009E73", "Streamline SC–FC nodes"= "#009E73")
# Create labeltext matrix 
labeltext_matrix <- cbind( forestdata$level, forestdata$asterisks ) 
# Assign colors for boxes, lines, and labels 
box_colors <- sapply(forestdata$level, function(level)
  level_colors[level]);
label_gpar <- lapply(1:nrow(forestdata), function(i) { list(
  level = gpar(fontsize = 10, fontface = "bold", col = level_colors[forestdata$level[i]]), 
  asterisks = gpar(fontsize = 15, fontface = "bold", col = level_colors[forestdata$level[i]]))})
# Update your plot code
plot <- forestplot( 
  labeltext = labeltext_matrix, 
  mean = forestdata$beta, 
  lower = forestdata$lower, 
  upper = forestdata$upper, 
  fn.ci_norm = fpDrawCircleCI, #turns boxes to circles
  vertices = TRUE, #adds whiskers
  col = fpColors( box = box_colors, line = box_colors, zero = "gray50" ), 
  shapes_gp = fpShapesGp( default = gpar(pch = 16, col = "black"), 
                          box = lapply(1:nrow(forestdata), function(i) gpar(col = box_colors[i], fill = box_colors[i])),
                          line = lapply(1:nrow(forestdata), function(i) gpar(col = box_colors[i], fill = box_colors[i])) ), 
  txt_gp = fpTxtGp( label= lapply(1:nrow(forestdata), function(i) list( gpar(fontsize = 10, fontface = "bold", col = level_colors[forestdata$level[i]]), gpar(fontsize = 10) )), 
                    ticks = gpar(fontsize = 15), 
                    xlab = gpar(fontsize = 15, fontface = "bold", col = "black"), 
                    title = gpar(fontsize = 15, fontface = "bold", col = "black") ), 
  xlab = "Adjusted beta", 
  title = title)
                    
plotlist=append(plotlist, list(grid::grid.grabExpr(print(plot))))
}


library(grid)
library(ggplot2)
library(gridExtra)

allplots=do.call(grid.arrange, c(plotlist, ncol = 1))

#saves plot 
ggsave("figures/DMN_coupling_cognition.png", plot = allplots, width = 2000, height = 1990, units = "px", dpi = 300)
