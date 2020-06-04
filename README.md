# CNV_detection

## Background:
You are developing an NGS​-based assay for a recently characterized disorder called Deldupemia, an autosomal recessive disease caused by mutations in the CNSL gene. Due to high sequence homology in the CNSL region, deletions and duplications are common (thus, rather than all patients having two copies of DNA across the region, some have one, and others have three). The deletion/duplication breakpoints can vary from sample to sample (see table further below), and it has been hypothesized that the breakpoints correspond with ethnicity.
With help from someone on the molecular biology team, you developed probes for the CNSL region and performed a retrospective analysis of 10,000 samples spanning different ethnicities. The data in “cnsl_data.csv.gz” catalogs the depth of NGS reads at 100 different hybrid capture probe locations, 50 in the CNSL region, and 50 outside of the region. The depth at a probe is expected to be linearly proportional to the copy number of the DNA at that site (i.e., having three CNSL copies should give roughly 3x the depth as having one copy).

## Assumptions:
● Probes in the “non​CNSL” region are expected to have CN=2.

● A deletion or duplication is any contiguous stretch of at least four well behaved probes
that have copy number of ~1 or ~3, respectively.

● Due to variability of extraction efficiency in the lab and error in the quantification of DNA
libraries, each sample has a slightly different average NGS read depth across all probes.

● Each probe captures DNA with different efficiency relative to other probes, but you can
assume that a single probe is equally efficient across all samples.


## Objective

Characterize the deletion and duplication frequencies and breakpoint positions on a per-​ethnicity basis.


## Algorithm


● Step 1: filter the prob with less than 10 depth in at least 10 sample (CNSL_probe_5 is filtered in
this case)

● Step 2: fit the depth to a statistical model descripted on the paper named "Development and validation of a 36-gene sequencing assay for hereditary cancer risk assessment". In brief, the likelihood of observing a given number of mapped reads di,j at a given genomic position i for sample j with copy number ci,j is modeled as following:
p(di,j |ci,j ) = NegBinom(di,j |μ = ci,j μiμj ,r = ri)
where μi is the average depth for that targeted location across samples, μj is the average depth for that particular sample across targeted positions, ri is the dispersion parameter.

● Step 3: Assign the copy number with largest probability to each sequencing depth.

● Step 4: A deletion or duplication is any contiguous stretch of at least four well behaved probes that have copy number of ~1 or ~3, respectively.
