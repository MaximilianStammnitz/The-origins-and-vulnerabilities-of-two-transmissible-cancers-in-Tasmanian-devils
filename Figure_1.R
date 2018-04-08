# Figure 1: Origins of DFT1 and DFT2

# Max Stammnitz, Transmissible Cancer Group, 
# mrs72@cam.ac.uk

# Please cite:
# "The origins and vulnerabilities of two transmissible cancers in Tasmanian devils"
# Cancer Cell 33(4), 607-619 (2018).

# Libaries
# install.packages('dendextend')
library(dendextend)

# Colors
DFT1.col <- 'cornflowerblue'
DFT2.col <- 'red'
heat.cols <- c("white", "gray80", "black")

# 1. Load genotype table from Mendeley Data:
# ###

# 2. Only chose one SNP per (linked) RadSeq locus
first.SNP.over.linkage <- function(x){
  
  # a. Pre-format the data table
  x[,-c(1:7)] <- apply(x[,-c(1:7)], 2, 
                       function(y) 
                       {
                         y[y=='1\\1'] = 0; 
                         y[y=='1\\2' | y=='2\\1'] = 0.5; # 1\2\2 calls counted as HETs
                         y[y=='2\\2'] = 1;
                         y[y=='1\\2\\2'] = NA;
                         y})
  
  
  # b. Make 'iterations' random combinations of linked SNP loci
  
  # Isolate constant body of unique RadSeq-SNPs
  const.loci <- names(table(x[,'Locus No.'])[as.numeric(table(x[,'Locus No.']))==1])
  const.loci.snps <- x[match(as.numeric(const.loci),x[,'Locus No.']),]
  
  # Isolate variable body of internally linked RadSeq-SNPs
  linked.loci <- names(table(x[,'Locus No.'])[as.numeric(table(x[,'Locus No.']))>1])
  linked.loci.snps <- vector(mode='list', length=length(unique(linked.loci)))
  names(linked.loci.snps) <- linked.loci
  for (i in 1:length(linked.loci.snps)){
    linked.loci.snps[[i]] <- x[which(x[,'Locus No.']==as.numeric(linked.loci[[i]])),]
  }
  
  # Take first elements from each body of linked SNPs
  linked.loci.snps.firsts <- as.character(sapply(linked.loci.snps, function(x) {rownames(x[1,])}))
  linked.loci.snps.firsts <- match(linked.loci.snps.firsts,rownames(x))
  
  # Obtain variable matrizes (vectorised)
  x_firsts <- x[linked.loci.snps.firsts,]
  
  # Reunite variable and constant bodies (vectorised)
  y <- rbind(const.loci.snps,x_firsts)
  
  # Sort complete lists of SNPs  (vectorised)
  y <- y[order(y[,'CHR'],y[,'POS']),]
#   rownames(y)<- paste(y[,'CHR'],
#                       y[,'Locus No.'],
#                       sep=":Loc:")
  y <- y[,-c(1:7)]
  y <- as.matrix(y)
  
  # c. Output
  return(y)
}
genotypes <- first.SNP.over.linkage(genotypes)

# 3. Impute NA-values (see Methods)
impute.linked.snp <- function(x){
  
  # a. Pre-format the data table
  x <- apply(x, 2, function(y) {y[y=='1\\1'] = 0; 
                                y[y=='1\\2' | y=='2\\1'] = 0.5; # 1\2\2 calls counted as HETs
                                y[y=='2\\2'] = 1;
                                y[y=='1\\2\\2'] = NA;
                                y})
  
  # b. Create euclidian SNP distance matrix (yet excluding NAs from calculations)
  dist.SNPs <- as.matrix(dist(x))
  colnames(dist.SNPs) <- rownames(dist.SNPs)
  
  # c. Iterate over each sample and impute NA values
  for (i in 1:ncol(x)){
    
    # Find NAs in sample
    # cat('\n', i, '/', ncol(x))
    impute <- which(is.na(x[,i])==T)
    
    if(length(impute)!=0){
      
      # Impute each NA value based on shortest SNP distance across all samples (i.e. closest linkage)
      for (j in 1:length(impute)){
        
        # besides 'self-distance of 0'
        linkage.degrees <- names(sort(dist.SNPs[which(rownames(x)[impute[j]]==rownames(dist.SNPs)),])[-1])
        
        # Fill in imputed value from this particular sample i
        for (k in 1:length(linkage.degrees)){
          
          # take closest linked SNP genotype that is not NA
          if(is.na(x[which(rownames(x)==linkage.degrees[k]),i])==F){
            x[as.numeric(impute[j]),i] <- x[which(rownames(x)==linkage.degrees[k]),i]
            break
          }
        }
      }  
    }
  }
  
  # d. Output
  return(x)
}
genotypes <- impute.linked.snp(genotypes)
class(genotypes) <- 'numeric'

# 4. Rename samples
rename.samples <- function(x){
  
  # a. Location
  colnames(x) <- gsub('AR_', 'Arthur River ', colnames(x))
  colnames(x) <- gsub('WN_', 'Woolnorth ', colnames(x))
  colnames(x) <- gsub('NP_', 'Narawntapu ', colnames(x))
  colnames(x) <- gsub('MW_', 'Mount William ', colnames(x))
  colnames(x) <- gsub('FN_', 'Freycinet ', colnames(x))
  colnames(x) <- gsub('FO_', 'Forestier ', colnames(x))
  
  # b. Year
  for (i in 1:length(c(1999:2013))){
    year <- c(1999:2013)[i]
    colnames(x) <- gsub(paste0(substr(year, 3, 4),'_'), paste0(year,''), colnames(x))
  }
  
  # c. Clip off gender/sample number
  colnames(x)[c(1:398)] <- substr(colnames(x)[c(1:398)],1,nchar(colnames(x)[c(1:398)])-4)
  
  # d. Rename Channel samples
  colnames(x) <- gsub('202H1', 'Channel 2014', colnames(x))
  colnames(x) <- gsub('203H', 'Channel 2014', colnames(x))
  
  # e. Output
  return(x)
}
genotypes.orig <- genotypes
genotypes <- rename.samples(genotypes)

# 5. Plot clustering map
samples <- genotypes %>% t %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", value = 'black', k = 1) %>% 
  set("branches_lwd", 0.5)

snps <- genotypes %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", value = 'black', k = 1)

pdf('Figure_1.pdf',
    width=50, height=50)
mar.default <- c(5,5,5,5) + 0.5
par(mar = mar.default)
genotypes.plot <- genotypes
heatmap(x = genotypes.plot, 
        Rowv=snps, 
        Colv=samples, 
        col=heat.cols,
        na.rm=T,
        margins = c(10, 15))
dev.off()