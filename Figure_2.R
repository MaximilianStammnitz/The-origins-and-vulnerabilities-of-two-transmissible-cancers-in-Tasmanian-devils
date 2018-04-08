# Figure 2: Single-Nucleotide Variants and Indels in DFT1 and DFT2

# Max Stammnitz, Transmissible Cancer Group, University of Cambridge, UK
# mrs72@cam.ac.uk

# Please cite:
# "The origins and vulnerabilities of two transmissible cancers in Tasmanian devils"
# Cancer Cell 33(4), 607-619 (2018).

# Libaries
library(stringr)
library(scales)
library(sigfit)
library(coda)
library(rstan)

# select colors
DFT1.col <- 'cornflowerblue'
DFT2.col <- 'red'
heat.cols <- c("white", "gray80", "black")

# 1. Load SNV data and Mutational Signatures sets
# DFT1.somatic.SNVs <- 
# DFT2.somatic.SNVs <- 

# 2. Generate spectra
MUT.spectrum.96 <- function(x){
  
  # a. Isolate triplets form Platypus VCF
  context <- rapply(strsplit(as.character(x[,"INFO"]),";"), function(x) x[12])
  context <- str_split_fixed(context,"=",2)[,2]
  context <- rapply(strsplit(context,""), function(x) paste(x[10:12],collapse=""))
  
  # b. Isolate base changes
  context.changes <- matrix(NA,nrow=length(context),ncol=3)
  colnames(context.changes) <- c("REF", "ALT-5'", "ALT-3'")
  context.changes[,1] <- context
  context.changes[,2:3] <- str_split_fixed(context.changes[,1],"",3)[,c(1,3)]
  context.changes[,2] <- paste(context.changes[,2],as.character(x[,'ALT']),context.changes[,3],sep="")
  context.changes <- context.changes[,c(1,2)]
  colnames(context.changes) <- c("REF", "ALT")
  
  # Pyrimidine-context substitutions
  equivalent.mut <- matrix(c('ACA>AAA', 'TGT>TTT', # AC - AA, or TG - TT ###  C>A or G>T
                             'ACC>AAC', 'GGT>GTT',
                             'ACG>AAG', 'CGT>CTT',
                             'ACT>AAT', 'AGT>ATT',
                             'CCA>CAA', 'TGG>TTG', # CC - CA, or GG - GT
                             'CCC>CAC', 'GGG>GTG',
                             'CCG>CAG', 'CGG>CTG',
                             'CCT>CAT', 'AGG>ATG',
                             'GCA>GAA', 'TGC>TTC', # GC - GA, or CG - CT
                             'GCC>GAC', 'GGC>GTC',
                             'GCG>GAG', 'CGC>CTC',
                             'GCT>GAT', 'AGC>ATC',
                             'TCA>TAA', 'TGA>TTA', # TC - TA, or AG - AT
                             'TCC>TAC', 'GGA>GTA',
                             'TCG>TAG', 'CGA>CTA',
                             'TCT>TAT', 'AGA>ATA',
                             'ACA>AGA', 'TGT>TCT', # AC - AG, or TG - TC ###  C>G or G>C
                             'ACC>AGC', 'GGT>GCT',
                             'ACG>AGG', 'CGT>CCT',
                             'ACT>AGT', 'AGT>ACT',
                             'CCA>CGA', 'TGG>TCG', # CC - CG, or GG - GC
                             'CCC>CGC', 'GGG>GCG',
                             'CCG>CGG', 'CGG>CCG',
                             'CCT>CGT', 'AGG>ACG',
                             'GCA>GGA', 'TGC>TCC', # GC - GG, or CG - CC
                             'GCC>GGC', 'GGC>GCC',
                             'GCG>GGG', 'CGC>CCC',
                             'GCT>GGT', 'AGC>ACC',
                             'TCA>TGA', 'TGA>TCA', # TC - TG, or AG - AC
                             'TCC>TGC', 'GGA>GCA',
                             'TCG>TGG', 'CGA>CCA',
                             'TCT>TGT', 'AGA>ACA',
                             'ACA>ATA', 'TGT>TAT', # AC - AT, or TG - TA ###  C>T or G>A
                             'ACC>ATC', 'GGT>GAT',
                             'ACG>ATG', 'CGT>CAT',
                             'ACT>ATT', 'AGT>AAT',
                             'CCA>CTA', 'TGG>TAG', # CC - CT, or GG - GA
                             'CCC>CTC', 'GGG>GAG',
                             'CCG>CTG', 'CGG>CAG',
                             'CCT>CTT', 'AGG>AAG',
                             'GCA>GTA', 'TGC>TAC', # GC - GT, or CG - CA
                             'GCC>GTC', 'GGC>GAC',
                             'GCG>GTG', 'CGC>CAC',
                             'GCT>GTT', 'AGC>AAC',
                             'TCA>TTA', 'TGA>TAA', # TC - TT, or AG - AA
                             'TCC>TTC', 'GGA>GAA',
                             'TCG>TTG', 'CGA>CAA',
                             'TCT>TTT', 'AGA>AAA',
                             'ATA>AAA', 'TAT>TTT', # AT - AA, or TA - TT ###  T>A or A>T
                             'ATC>AAC', 'GAT>GTT',
                             'ATG>AAG', 'CAT>CTT',
                             'ATT>AAT', 'AAT>ATT',
                             'CTA>CAA', 'TAG>TTG', # CT - CA, or GA - GT
                             'CTC>CAC', 'GAG>GTG',
                             'CTG>CAG', 'CAG>CTG',
                             'CTT>CAT', 'AAG>ATG',
                             'GTA>GAA', 'TAC>TTC', # GT - GA, or CA - CT
                             'GTC>GAC', 'GAC>GTC',
                             'GTG>GAG', 'CAC>CTC',
                             'GTT>GAT', 'AAC>ATC',
                             'TTA>TAA', 'TAA>TTA', # TT - TA, or AA - AT
                             'TTC>TAC', 'GAA>GTA',
                             'TTG>TAG', 'CAA>CTA',
                             'TTT>TAT', 'AAA>ATA',
                             'ATA>ACA', 'TAT>TGT', # AT - AC, or TA - TG ### T>C or A>G
                             'ATC>ACC', 'GAT>GGT',
                             'ATG>ACG', 'CAT>CGT',
                             'ATT>ACT', 'AAT>AGT',
                             'CTA>CCA', 'TAG>TGG', # CT - CC, or GA - GG
                             'CTC>CCC', 'GAG>GGG',
                             'CTG>CCG', 'CAG>CGG',
                             'CTT>CCT', 'AAG>AGG',
                             'GTA>GCA', 'TAC>TGC', # GT - GC, or CA - CG
                             'GTC>GCC', 'GAC>GGC',
                             'GTG>GCG', 'CAC>CGC',
                             'GTT>GCT', 'AAC>AGC',
                             'TTA>TCA', 'TAA>TGA', # TT - TC, or AA - AG
                             'TTC>TCC', 'GAA>GGA',
                             'TTG>TCG', 'CAA>CGA',
                             'TTT>TCT', 'AAA>AGA',
                             'ATA>AGA', 'TAT>TCT', # AT - AG, or TA - TC ### T>G or A>C
                             'ATC>AGC', 'GAT>GCT',
                             'ATG>AGG', 'CAT>CCT',
                             'ATT>AGT', 'AAT>ACT',
                             'CTA>CGA', 'TAG>TCG', # CT - CG, or GA - GC
                             'CTC>CGC', 'GAG>GCG',
                             'CTG>CGG', 'CAG>CCG',
                             'CTT>CGT', 'AAG>ACG',
                             'GTA>GGA', 'TAC>TCC', # GT - GG, or CA - CC
                             'GTC>GGC', 'GAC>GCC',
                             'GTG>GGG', 'CAC>CCC',
                             'GTT>GGT', 'AAC>ACC',
                             'TTA>TGA', 'TAA>TCA', # TT - TG, or AA - AC
                             'TTC>TGC', 'GAA>GCA',
                             'TTG>TGG', 'CAA>CCA',
                             'TTT>TGT', 'AAA>ACA'),
                           nrow=96, ncol=2, byrow=TRUE)
  
  # c. Define variants from input table and convert them to pyrimidine context
  triplett_snvs <- paste0(context.changes[,1], '>', context.changes[,2])
  for(mut in 1:nrow(equivalent.mut)){
    triplett_snvs <- gsub(equivalent.mut[mut,2], equivalent.mut[mut,1], triplett_snvs)
  }
  
  # d. Count triplets
  consensus.mut <- equivalent.mut[,1]
  counts <- table(triplett_snvs)
  if(length(counts)<96){
    add.names <- consensus.mut[which(is.na(match(consensus.mut,names(counts)))==T)]
    names.takeover <- names(counts)
    counts <- c(counts,rep(0,length(add.names)))
    names(counts) <- c(names.takeover,add.names)
  }
  counts <- counts[match(consensus.mut,names(counts))]
  
  # e. Output
  names(counts) <- consensus.mut
  return(counts)
}
DFT1.somatic.spectrum <- MUT.spectrum.96(DFT1.somatic.SNVs)
DFT2.somatic.spectrum <- MUT.spectrum.96(DFT2.somatic.SNVs)

# 3. Load Tasmanian devil base-triplet frequencies

# 4. Plot mutational spectra of DFT1 and DFT2
pdf("Figure_2a.pdf", height=12, width=20)
par(mar=c(5, 17, 4, 1), mgp=c(3,-0.5,-2.5))
colors = c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60")

DFT1.col <- 'cornflowerblue'
DFT2.col <- 'red'
mut <- c('C > A','C > G','C > T','T > A','T > C','T > G')

# Background colours
y.top = 90
borders = c(0, 306.9*1/6, 306.9*2/6, 306.9*3/6, 306.9*4/6, 306.9*5/6, 306.9)
alphas = c(rep(0.2,3), rep(0.3,3))
plot(1, type="n", ylim=c(0,y.top), xlim=c(0, 306.9), xlab="", ylab="", axes=F)
for (i in 1:6) {
  rect(xleft=borders[i], xright=borders[i+1], ybottom=0, 
       ytop=y.top, col=alpha(colors[i], alphas[i]), border="white")
  rect(xleft=borders[i], xright=borders[i+1], ybottom=y.top-7.5, 
       ytop=y.top, col=colors[i], border="white")
  text(x=(borders[i]+borders[i+1])/2, y=y.top-3.75, labels=mut[i], cex=3, col="white")
}

# Normalise to [%] mutations
devil.tripl.div <- c(rep(as.numeric(devil.7.1.genome.triplet.freq.32[,5])[1:16],3),
                     rep(as.numeric(devil.7.1.genome.triplet.freq.32[,5])[17:32],3))
DFT1.somatic.all.norm <- c(c(DFT1.somatic.spectrum/sum(DFT1.somatic.spectrum))/devil.tripl.div)*
  sum(devil.tripl.div)
DFT2.somatic.all.norm <- c(c(DFT2.somatic.spectrum/sum(DFT2.somatic.spectrum))/devil.tripl.div)*
  sum(devil.tripl.div)

# DFT1 spectrum
barplot(DFT1.somatic.all.norm, 
        col = rep(DFT1.col,96),
        width = 2,
        border = NA,
        main = '', 
        cex.main = 2,
        cex.names = 1.4,
        ylab = "", 
        yaxt = 'n', 
        ylim = c(0, 45),
        las = 2,
        font = 1,
        family = 'mono',
        xaxs = 'i',
        add = T,
        xlim = c(0, 306.9),
        space = c(rep(c(rep(0.593,15),0.6),5),
                  rep(0.593,16)),
        offset = 37,
        names.arg = NA)

# DFT2 spectrum
barplot(DFT2.somatic.all.norm, 
        col = rep(DFT2.col,96),
        width = 2,
        border = NA,
        main = '', 
        ylab = '', 
        yaxt = 'n',
        xaxt = 'n',
        ylim = c(0, 40),
        xaxs = 'i',
        add = T,
        xlim = c(0, 279*1.1),
        space = c(rep(c(rep(0.593,15),0.6),5),
                 rep(0.593,16)),
        names.arg = NA)

# Axes and Labels
axis(side=2, las=2, cex.axis=3, col=DFT2.col, col.axis=DFT2.col,
     at = c(0, 10), labels = c('0', ''), lwd = 5, pos=-2, hadj = 2.1)
axis(side=2, las=2, cex.axis=3, col=DFT2.col, col.axis=DFT2.col,
     at = c(10, 20, 30), labels = c('10', '20', '30'), lwd = 5, pos=-2, hadj = 1.6)
axis(side=2, las=2, cex.axis=3, col=DFT1.col, col.axis=DFT1.col,
     at = c(37, 47), labels = c('0', ''), lwd = 5, pos=-2, hadj = 2.1)
axis(side=2, las=2, cex.axis=3, col=DFT1.col, col.axis=DFT1.col,
     at = c(47, 57, 67, 77), labels = c('10', '20', '30', '40'), lwd = 5, pos=-2, hadj = 1.6)
title(ylab="SNV Proportion [%]", line=6, cex.lab=5)
dev.off()

# 5. Fitting species-agnostic COSMIC mutational signatures to DFT1 and DFT2 somatic

# Import human (hg19) genome triplet frequencies
# hg19.genome.triplet.freq.32 <-
  
# Import devil genome triplet frequencies (masked genome with regions filters)
# devil.7.1.genome.triplet.freq.32 <- 

# Calculate devil triplet to human triplet frequency ratios
normaliser <- c(rep(as.numeric(devil.7.1.genome.triplet.freq.32[,5])[1:16],3),
                rep(as.numeric(devil.7.1.genome.triplet.freq.32[,5])[17:32],3)) / 
  c(rep(as.numeric(hg19.genome.triplet.freq.32[,5])[1:16],3),
    rep(as.numeric(hg19.genome.triplet.freq.32[,5])[17:32],3))

# Normalise 30 human cosmic signatures by devil-human triplet frequency ratio
signature.types.devil <- apply(signature.types.human, 2, function(sig) {
  new.sig = sig * normaliser
  new.sig / sum(new.sig)
})
signature.types.devil <- signature.types.devil + 1e-9

## Unfiltered SNV/SNP counts
MCMC.Counts <- matrix(nrow = 96, ncol = 2, data = NA)
colnames(MCMC.Counts) <- c('DFT1-somatic', 'DFT2-somatic')
MCMC.Counts[,1] <- DFT1.somatic.spectrum
MCMC.Counts[,2] <- DFT2.somatic.spectrum

DFT1.somatic.fit <- sigfit::run_sampling(counts = MCMC.Counts[,'DFT1-somatic'], 
                                         signatures = signature.types.devil[,c(1,5)],
                                         iter = 100000, warmup = 1000)
DFT2.somatic.fit <- sigfit::run_sampling(counts = MCMC.Counts[,'DFT2-somatic'], 
                                         signatures = signature.types.devil[,c(1,5)],
                                         iter = 100000, warmup = 1000)

## Extract means and HDP-intervals
DFT1.somatic.fit <- extract(DFT1.somatic.fit)
means.DFT1 <- apply(DFT1.somatic.fit$exposures, 2, mean)
lower.DFT1 <- HPDinterval(as.mcmc(DFT1.somatic.fit$exposures), prob=0.95)[,1]
upper.DFT1 <- HPDinterval(as.mcmc(DFT1.somatic.fit$exposures), prob=0.95)[,2]
names(means.DFT1) <- names(lower.DFT1) <- names(upper.DFT1) <- c('1H', '5H')
DFT1.somatic.fit <- cbind(means.DFT1,lower.DFT1,upper.DFT1)

DFT2.somatic.fit <- extract(DFT2.somatic.fit)
means.DFT2 <- apply(DFT2.somatic.fit$exposures, 2, mean)
lower.DFT2 <- HPDinterval(as.mcmc(DFT2.somatic.fit$exposures), prob=0.95)[,1]
upper.DFT2 <- HPDinterval(as.mcmc(DFT2.somatic.fit$exposures), prob=0.95)[,2]
names(means.DFT2) <- names(lower.DFT2) <- names(upper.DFT2) <-c('1H', '5H')
DFT2.somatic.fit <- cbind(means.DFT2,lower.DFT2,upper.DFT2)

## Plot
pdf('Figure_2b.pdf',
    width=12, height=12)
par(mar=c(5, 17, 5, 4), mgp=c(3,-0.5,-2.5))
DFT1.col <- 'cornflowerblue'
DFT2.col <- 'red'
fits <- barplot(rbind(DFT1.somatic.fit[,1],
                      DFT2.somatic.fit[,1]), 
                beside = T,
                col = c(DFT1.col, DFT2.col),
                cex.lab = 2,
                cex.axis = 1.5,
                ylab = '',
                border='white',
                cex.names = 2.5,
                names.arg = c('Signature 1', 'Signature 5'),
                yaxt = 'n',
                ylim = c(0,1),
                xaxt = 'n')
arrows(fits, rbind(DFT1.somatic.fit[,2],DFT2.somatic.fit[,2]), 
       fits, rbind(DFT1.somatic.fit[,3],DFT2.somatic.fit[,3]), 
       angle=90, code=1, length=0.2, lwd=5, col='black')
arrows(fits, rbind(DFT1.somatic.fit[,3],DFT2.somatic.fit[,3]),
       fits, rbind(DFT1.somatic.fit[,1],DFT2.somatic.fit[,1]), 
       angle=90, code=1, length=0.2, lwd=5, col='black')
legend('topleft',
       legend = c('DFT1', 'DFT2'),
       col = c(DFT1.col, DFT2.col),
       pch = 15,
       cex = 5, 
       bty = 'n',
       text.col = c(DFT1.col, DFT2.col))
axis(side = 1, at = c(2,5), tick = FALSE, 
     labels = c('Signature 1', 'Signature 5'), cex.axis = 4,
     padj = 1.5)
axis(2, at = c(0,0.25), labels = c('0',''),
     cex.axis = 3, las = 2, lwd = 5, pos=0.5, hadj = 2)
axis(2, at = c(0.75,1), labels = c('','100'),
     cex.axis = 3, las = 2, lwd = 5, pos=0.5, hadj = 1.38)
axis(2, at = seq(f=0, t=1, length.out=5), 
     labels = c('', '25', '50', '75', ''),
     cex.axis = 3, las = 2, lwd = 5, pos=0.5, hadj = 1.6)
title(ylab = "Signature Proportion [%]", line = 10, cex.lab = 5)
dev.off()