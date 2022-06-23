library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(grid)

#Generating barplot of all OTUs from all samples

#Read in data
data <- read.csv(file = "All_OTUs.txt", header=T, sep="\t")

#convert data frame from a "wide" format to a "long" format
melted_data = melt(data, id = c("Sample"))

#Keep samples ordered the way they're ordered in the data
melted_data$Sample <- factor(melted_data$Sample,levels=unique(melted_data$Sample))

#Making a nice palette with distinctive colours
library(randomcoloR)
n <- 750
palette <- distinctColorPalette(n)

#Generating plot
barplot = ggplot(melted_data, aes(x = Sample, fill = forcats::fct_rev(variable), y = value)) + 
    geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1), 
    axis.title.y = element_text(size = 10), legend.title = element_text(size = 16), 
    axis.text.y = element_text(colour = "black", size = 10, face = "bold")) + 
    scale_y_continuous(expand = c(0,0)) + 
    theme(legend.position = "none") +
    scale_fill_manual(values = rep(palette,10))+
    labs(x = "", y = "Relative Abundance (%)", fill = "OTU")

barplot

#Generating barplot of conserved OTUs from all samples
#Read in data
conserved_data <- read.csv(file = "Conserved_OTUs.txt", header=T, sep="\t")

#convert data frame from a "wide" format to a "long" format
conserved_melted_data = melt(conserved_data, id = c("Sample"))

#Keep samples ordered the way they're ordered in the data
conserved_melted_data$Sample <- factor(conserved_melted_data$Sample,levels=unique(conserved_melted_data$Sample))

#Making a nice palette with distinctive colours
library(randomcoloR)
n <- 7
palette <- distinctColorPalette(n)

#Generating plot
conserved_barplot = ggplot(conserved_melted_data, aes(x = Sample, fill = forcats::fct_rev(variable), y = value)) + 
    geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1), 
    axis.title.y = element_text(size = 10), legend.title = element_text(size = 10), 
    axis.text.y = element_text(colour = "black", size = 10, face = "bold")) + 
    scale_y_continuous(expand = c(0,0)) + 
    scale_fill_manual(values = palette)+
    labs(x = "", y = "Relative Abundance (%)", fill = "OTU")

conserved_barplot

#Generate 3D NMDS plot from all OTUs in 2019
abund = read.csv("NMDS_data.txt", header= TRUE, sep="\t")

#make community matrix - extract columns with abundance information
nums = abund[,4:ncol(abund)]

#turn abundance data frame into a matrix
nums_matrix = as.matrix(nums)

#Now to create our 3D plot
set.seed(123)
nmds = metaMDS(nums_matrix, k=3, distance = "bray")

#extract NMDS scores (x,y and z coordinates)
data.scores = as.data.frame(scores(nmds))

#add columns to data frame 
data.scores$Sample = abund$Sample
data.scores$Time = abund$Time
data.scores$Location = abund$Location
 
head(data.scores)

#Load up libraries
library(plotly)

#Generate plot - points coloured by location and shaped by time
plot_3D <- plot_ly(data.scores, x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, color = ~Location, symbol = ~Time, symbols = c('circle','diamond','square'), marker = list(size = 5), text = ~paste('<br>Sample:', Sample, '<br>Location:', Location)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'NMDS1'),
                     yaxis = list(title = 'NMDS2'),
                     zaxis = list(title = 'NMDS3')))
#Show plot
plot_3D

#Generate plot - points coloured by time and shaped by location
plot_3D2 <- plot_ly(data.scores, x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, color = ~Time, symbol = ~Location, symbols = c('circle','diamond','square','cross'), marker = list(size = 5), text = ~paste('<br>Sample:', Sample, '<br>Location:', Location)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'NMDS1'),
                     yaxis = list(title = 'NMDS2'),
                     zaxis = list(title = 'NMDS3')))
#Show plot
plot_3D2

#Calculate pairwise ANOSIM scores, grouped by location

#Find unique groups in data
groups = unique(abund[c("Location")])
group_list = as.list(groups$Location)

#Calculate pairwise ANOSIM statistics
for (i in group_list) {
    for (j in group_list) {
    if (i==j) next
    sub_df1 = subset(abund, subset=(Location==i))
    sub_df2 = subset(abund, subset=(Location==j))
    sub_df = rbind(sub_df1, sub_df2)
    nums = sub_df[,4:ncol(sub_df)]
    nums_matrix = as.matrix(nums)
    set.seed(123)
    ano = anosim(nums_matrix, sub_df$Location, distance = "bray", permutations = 9999)
    Rvalue = ano$statistic
    pvalue = ano$signif
    cat(i,"\t",j ,"\tR-value:\t", Rvalue, "\tP-value:\t", pvalue, '\n')

    }
}

#Calculate pairwise ANOSIM scores, grouped by time

#Find unique groups in data
groups = unique(abund[c("Time")])
group_list = as.list(groups$Time)

#Calculate pairwise ANOSIM statistics
for (i in group_list) {
    for (j in group_list) {
    if (i==j) next
    sub_df1 = subset(abund, subset=(Time==i))
    sub_df2 = subset(abund, subset=(Time==j))
    sub_df = rbind(sub_df1, sub_df2)
    nums = sub_df[,4:ncol(sub_df)]
    nums_matrix = as.matrix(nums)
    set.seed(123)
    ano = anosim(nums_matrix, sub_df$Time, distance = "bray", permutations = 9999)
    Rvalue = ano$statistic
    pvalue = ano$signif
    cat(i,"\t",j ,"\tR-value:\t", Rvalue, "\tP-value:\t", pvalue, '\n')

    }
}

#Generate diversity boxplots faceted by time

#Generate 3D NMDS plot from all OTUs in 2019
abund = read.csv("NMDS_data.txt", header= TRUE, sep="\t")

#make community matrix - extract columns with abundance information
nums = abund[,4:ncol(abund)]

#turn abundance data frame into a matrix
nums_matrix = as.matrix(nums)

#Calculate shannon diversity
shan = diversity(nums_matrix, index="shannon")
shan.scores = as.data.frame(scores(shan))
colnames(shan.scores) = c('Shannon')

#Calculate simpson diversity
simp = diversity(nums_matrix, index="simpson")
simp.scores = as.data.frame(scores(simp))
colnames(simp.scores) = c('Simpson')

#Create main dataframe with metadata
scores_df <- abund[c(1,2,3)]

#Add in shannon and simpson diversity scores
scores_df$Shannon = shan.scores$Shannon
scores_df$Simpson = simp.scores$Simpson


box_data = melt(scores_df, id = "Location")
sub_box = subset(box_data,variable !='Sample')
sub_box$value = as.numeric(as.character(sub_box$value))

#Plot Shannon diversity box plots
shan_box_plot = ggplot(scores_df, aes(x = Location, y = Shannon, fill = Time)) + 
     geom_boxplot(colour = "black", position = position_dodge(1)) +
     facet_wrap(~Location,ncol = 4, scales = "free_x") 

shan_box_plot

#Plot Simpson diversity box plots
simp_box_plot = ggplot(scores_df, aes(x = Location, y = Simpson, fill = Time)) + 
     geom_boxplot(colour = "black", position = position_dodge(1)) +
     facet_wrap(~Location,ncol = 4, scales = "free_x") 

simp_box_plot

#Fitting significant taxa to NMDS plot

#Code taken from https://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html and adapted to suit what we needed
 
bv.step <- function(fix.mat, var.mat, 
                    fix.dist.method="bray", var.dist.method="euclidean", correlation.method="spearman",
                    scale.fix=FALSE, scale.var=TRUE,
                    max.rho=0.95,
                    min.delta.rho=0.001,
                    random.selection=TRUE,
                    prop.selected.var=0.2,
                    num.restarts=10,
                    var.always.include=NULL,
                    var.exclude=NULL,
                    output.best=10
){
 
  if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
  if(sum(var.always.include %in% var.exclude) > 0){stop("var.always.include and var.exclude share a variable")}
  require(vegan)
 
  if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
  if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}
 
  fix.dist <- vegdist(as.matrix(fix.mat), method=fix.dist.method)
 
  #an initial removal phase
  var.dist.full <- vegdist(as.matrix(var.mat), method=var.dist.method)
  full.cor <- suppressWarnings(cor.test(fix.dist, var.dist.full, method=correlation.method))$estimate
  var.comb <- combn(1:ncol(var.mat), ncol(var.mat)-1)
  RES <- data.frame(var.excl=rep(NA,ncol(var.comb)), n.var=ncol(var.mat)-1, rho=NA)
  for(i in 1:dim(var.comb)[2]){
    var.dist <- vegdist(as.matrix(var.mat[,var.comb[,i]]), method=var.dist.method)
    temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
    RES$var.excl[i] <- c(1:ncol(var.mat))[-var.comb[,i]]
    RES$rho[i] <- temp$estimate
  }
  delta.rho <- RES$rho - full.cor
  exclude <- sort(unique(c(RES$var.excl[which(abs(delta.rho) < min.delta.rho)], var.exclude)))
 
  if(random.selection){
    num.restarts=num.restarts
    prop.selected.var=prop.selected.var
    prob<-rep(1,ncol(var.mat))
    if(prop.selected.var< 1){
      prob[exclude]<-0
    }
    n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
  } else {
    num.restarts=1
    prop.selected.var=1  
    prob<-rep(1,ncol(var.mat))
    n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
  }
 
  RES_TOT <- c()
  for(i in 1:num.restarts){
    step=1
    RES <- data.frame(step=step, step.dir="F", var.incl=NA, n.var=0, rho=0)
    attr(RES$step.dir, "levels") <- c("F","B")
    best.comb <- which.max(RES$rho)
    best.rho <- RES$rho[best.comb]
    delta.rho <- Inf
    selected.var <- sort(unique(c(sample(1:dim(var.mat)[2], n.selected.var, prob=prob), var.always.include)))
    while(best.rho < max.rho & delta.rho > min.delta.rho & RES$n.var[best.comb] < length(selected.var)){
      #forward step
      step.dir="F"
      step=step+1
      var.comb <- combn(selected.var, RES$n.var[best.comb]+1, simplify=FALSE)
      if(RES$n.var[best.comb] == 0){
        var.comb.incl<-1:length(var.comb)
      } else {
        var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
        temp <- NA*1:length(var.comb)
        for(j in 1:length(temp)){
          temp[j] <- all(var.keep %in% var.comb[[j]]) 
        }
        var.comb.incl <- which(temp==1)
      }
 
      RES.f <- data.frame(step=rep(step, length(var.comb.incl)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]+1, rho=NA)
      for(f in 1:length(var.comb.incl)){
        var.incl <- var.comb[[var.comb.incl[f]]]
        var.incl <- var.incl[order(var.incl)]
        var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
        temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
        RES.f$var.incl[f] <- paste(var.incl, collapse=",")
        RES.f$rho[f] <- temp$estimate
      }
 
      last.F <- max(which(RES$step.dir=="F"))
      RES <- rbind(RES, RES.f[which.max(RES.f$rho),])
      best.comb <- which.max(RES$rho)
      delta.rho <- RES$rho[best.comb] - best.rho 
      best.rho <- RES$rho[best.comb]
 
      if(best.comb == step){
        while(best.comb == step & RES$n.var[best.comb] > 1){
          #backward step
          step.dir="B"
          step <- step+1
          var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
          var.comb <- combn(var.keep, RES$n.var[best.comb]-1, simplify=FALSE)
          RES.b <- data.frame(step=rep(step, length(var.comb)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]-1, rho=NA)
          for(b in 1:length(var.comb)){
            var.incl <- var.comb[[b]]
            var.incl <- var.incl[order(var.incl)]
            var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
            temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
            RES.b$var.incl[b] <- paste(var.incl, collapse=",")
            RES.b$rho[b] <- temp$estimate
          }
          RES <- rbind(RES, RES.b[which.max(RES.b$rho),])
          best.comb <- which.max(RES$rho)
          best.rho<- RES$rho[best.comb]
        }
      } else {
        break()
      }
 
    }
 
    RES_TOT <- rbind(RES_TOT, RES[2:dim(RES)[1],])
    print(paste(round((i/num.restarts)*100,3), "% finished"))
  }
 
  RES_TOT <- unique(RES_TOT[,3:5])
 
 
  if(dim(RES_TOT)[1] > output.best){
    order.by.best <- RES_TOT[order(RES_TOT$rho, decreasing=TRUE)[1:output.best],]
  } else {
    order.by.best <-  RES_TOT[order(RES_TOT$rho, decreasing=TRUE), ]
  }
  rownames(order.by.best)<-NULL
 
  order.by.i.comb <- c()
  for(i in 1:length(selected.var)){
    f1 <- which(RES_TOT$n.var==i)
    f2 <- which.max(RES_TOT$rho[f1])
    order.by.i.comb <- rbind(order.by.i.comb, RES_TOT[f1[f2],])
  }
  rownames(order.by.i.comb)<-NULL
 
  if(length(exclude)<1){var.exclude=NULL} else {var.exclude=exclude}
  out <- list(
    order.by.best=order.by.best,
    order.by.i.comb=order.by.i.comb,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(order.by.best$var.incl[1], ",")))], collapse=","),
    best.model.rho=order.by.best$rho[1],
    var.always.include=var.always.include,
    var.exclude=var.exclude
  )
  out
 
}


abund_table = read.csv("Family_spe.txt", row.names=1, header= TRUE, sep="\t")
#Transpose the data to have sample names on rows
abund_table<-t(abund_table)
head(abund_table)
#Get grouping information
grouping_info = read.csv("Grouping_info.txt", row.names=1, header= TRUE, sep="\t")
head(grouping_info)
 
#Parameters
cmethod<-"pearson" #Correlation method to use: pearson, pearman, kendall
fmethod<-"bray" #Fixed distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
vmethod<-"bray" #Variable distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
nmethod<-"bray" #NMDS distance method:  euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao

res.bv.step.biobio <- bv.step(wisconsin(abund_table), wisconsin(abund_table), 
                              fix.dist.method=fmethod, var.dist.method=vmethod,correlation.method=cmethod,
                              scale.fix=FALSE, scale.var=FALSE, 
                              max.rho=0.95, min.delta.rho=0.001,
                              random.selection=TRUE,
                              prop.selected.var=0.3,
                              num.restarts=10,
                              output.best=10,
                              var.always.include=NULL) 
 
#Get the 10 best subset of taxa
taxaNames<-colnames(abund_table)
bestTaxaFit<-""
for(i in (1:length(res.bv.step.biobio$order.by.best$var.incl)))
{
  bestTaxaFit[i]<-paste(paste(taxaNames[as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[i], split=",")))],collapse=' + '), " = ",res.bv.step.biobio$order.by.best$rho[i],sep="")
}
bestTaxaFit<-data.frame(bestTaxaFit)
colnames(bestTaxaFit)<-"Best combination of taxa with similarity score"
 
bestTaxaFit

#Generate NMDS plot
MDS_res=metaMDS(abund_table, distance = nmethod, k = 2, trymax = 50)
 
bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[1], ",")))
bio.fit <- envfit(MDS_res, abund_table[,bio.keep,drop=F], perm = 999)
 
bio.fit
 
#Get site information
df<-scores(MDS_res,display=c("sites"))
 
#Add grouping information
df<-data.frame(df,Type=grouping_info[rownames(df),1:2])
 
#Get the vectors for bioenv.fit
df_biofit<-scores(bio.fit,display=c("vectors"))
df_biofit<-df_biofit*vegan:::ordiArrowMul(df_biofit)
df_biofit<-as.data.frame(df_biofit)

#Manually fix water rows that started with numbers (confuses R)
df[c("Type", "Type.Time")][is.na(df[c("Type", "Type.Time")])] <- "2021 Apr"
df[c("Type", "Type.Location")][is.na(df[c("Type", "Type.Location")])] <- "Water"


#Draw samples
fitted_NMDS<-ggplot()+
geom_point(data=df,aes(NMDS1,NMDS2,colour=Type.Time,shape=Type.Location))+
geom_segment(data=df_biofit, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)+
geom_text(data=as.data.frame(df_biofit*1.1),aes(NMDS1, NMDS2, label = rownames(df_biofit)),color="#808080",alpha=0.5)+
theme_bw()

fitted_NMDS


