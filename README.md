# circRNA translation

### Info

Data and Rscripts used in analysis of RiboSeq and MS data for [Hansen, TB, "Signal and noise in circRNA translation", methods, 2021](https://www.sciencedirect.com/science/article/pii/S104620232100044X)


Use **R/riboseq_reads.R** to plot the RiboSeq reads from high and low quality RiboSeq datasets (RiboSeq quality scores are found in data/RiboSeq_qual.txt, and computed as described in [Stagsted et al, "Noncoding AUG circRNAs constitute an abundant and conserved subclass of circles", LSA, 2019](https://www.life-science-alliance.org/content/2/3/e201900398))

**R/riboseq_correlation.R** plots the association between RiboSeq coverage and quality for mRNAs, lncRNAs, and circRNAs

**R/MS.R** plots Mass-spec analysis on data from [Doll et al](https://www.nature.com/articles/s41467-017-01747-2) using annotated proteins from UniProt with custom putative circRNA-derived peptides and the SARS-CoV2 proteome.


### License

Copyright (C) 2021 ncRNALab.  