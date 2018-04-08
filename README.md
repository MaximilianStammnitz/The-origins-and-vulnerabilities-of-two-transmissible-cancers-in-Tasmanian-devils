# The origins and vulnerabilities of two transmissible cancers in Tasmanian devils.
R scripts utilised for the analysis of Tasmanian devil genomic analyses used in our study "The origins and vulnerabilities of two transmissible cancers in Tasmanian devils". Cancer Cell 33(4), 607-619 (2018).

# Data
1. Genomic data (BAM files) can be accessed via the European Nucleotide Archive: https://www.ebi.ac.uk/ena/data/view/PRJEB21902
2. Additional tables, FISH, genome assembly, IHC and mutation data is available via Mendeley data: http://dx.doi.org/10.17632/znfphvhmbv.1 

# Figure 1:
- integrate and process SNP data from 398 Tasmanian devils by Brüniche-Olsen et al., 2016
- integrate genotype calls of our samples
  - 86T, 88T (DFT1)
  - 202T2, 203T3 (DFT2)
  - 91H, 202H1, 203H (Normals)
- draw hierarchical clustering map

# Figure 2:
- integrate somatic DFT1 and DFT2 single-nucleotide variants (SNVs)
- draw somatic mutational spectrum of DFT1 (86T-unique & 88T-unique)
- draw somatic mutational spectrum of DFT2 (202T2-unique & 203T3-unique)
- normalise 30 COSMIC signatures to Tasmanian devil 7.1 reference genome base-triplet frequency
- fit DFT1 and DFT2 mutational spectrum to different combinations of COSMIC signatures

# Figure 3
- integrate DFT1 and DFT2 structural variant (SV) calls
- draw circos plots by sample
- draw barplots of SV repair classes by sample

# Figure 4
- integrate copy-number calls
- integrate gene dosage alteration calls
- draw Venn diagrams
- draw PDGFRA, PDGFRB and B2M Campbellgrams

# Figure 5 (soon available as a shiny app)
- integrate drug screening IC50 data from DFT1, DFT2 and human cancer cell lines (https://www.cancerrxgene.org/, Yang et al., 2014)
- process samples
- draw hierarchical clustering map
- draw cross-sample IC50 boxplots with sample-wise 'beeswarming', for drugs Afatinib, Axitinib, Sorafenib, Dasatinib, Talazoparib and AZD7762
