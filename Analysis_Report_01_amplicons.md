Analysis Report 1: Human Skin Microbiota Sexually Induced Differences
================
Katelyn Kukar
November 01, 2017

Introduction
============

Methodology to estimate microbial diversity has significantly increased with an effort to understand the distribution and diversity of microbes in their natural environments. Cultivation-dependent methodology works to mimic the natural environment of microbes in a laboratory setting, with enrichment cultures being particularly important for visualizing species that favor the selected media of growth (Nichols, 2007). Different ecological metrics have been used to categorize these cultivated colonies in order to deduce community data about known bacterial species, and potentially deduce information about unknown species (Birtel *et al.*, 2015). Understanding the cultivation methodology to successfully grow microbial species is imperative for gathering information abou bio diverse communities that can be beneficial in studying niche environments, like the human gut microbiota which interact with pathogens and probiotics (Wu *et al.*, 2012). Beyond healthcare, general dispersal of bacterial communities increases biodiversity in water sources, soils, and all micro environments that effect the well being of the larger ecosystem; with this assumption also holding true for bacterial organisms which reside on a human host or human host environment (Fuhrman, 1999).

Within the paper, *Forensic identification using skin bacterial communities*, published in 2010 by Noah Fierer the diversity in human skin bacterial communities is explained to be “far higher than previously recognized, with a high degree of inter-individual variability in the composition of bacterial communities.” This discovery leads the paper to discuss that these individualized communities can be utilized as personal identifiers for criminal cases as forensic evidence. Fierer et al. (2010) claims that “these skin bacteria may persist on touched surfaces for prolonged periods because many are highly resistant to environmental- mental stresses, including moisture, temperature, and UV radiation,” meaning the community samples can be collected long after the host has left (Fierer *et al.*, 2010). Though the implications for forensic research are still under scrutiny, the data can be more broadly analyzed to deduce if there are specific differences in bacterial communities between female and male participants. By analytically evaluating these communities separately -- to find patterns in abundance, richness, and diversity -- we deduce information about the specific community composition, which can yield information about the overall biodiversity. The inference that human hosts hold their own unique sets of bacteria helps foster the notion of unique micro environments that endosymbiotically regulate the human body.

By explicitly looking at community differences between male and female participants, information can be gathered about the overall health and well-being of each participant sex in regard to the surrounding environment. In a study published in Science in 2013, it was discovered that the gut microbiota are extremely relevant in determining autoimmune disease susceptibility (Markle *et al.*, 2013). Apparently, when male cecal contents were transferred to female mice these mice received a higher level of protection against “pancreatic islet inflammation, autoantibody production, and the development of diabetes,” meaning this difference within the sexes microbiota was translational to inhibit the spread of disease (Markle *et al.*, 2013). Therefore, understanding the differences in the skin microbiota could also potentially pose interesting medical questions that could yield insight on diseases current to specific, sexually biased traits.

Overall, divulging community data analysis from female and male participants, from the sequenced data sets from the Fierer et. al (2010) paper, can provide useful information about the resident skin microbiota of sexually diverse humans. By using sex as the means of a comparison further hypotheses can be made that explain differences based on hormonal differences and social differences since the presiding environment remains constant for both sample sets. We can thereby infer that because of these differences a community the female and male sex of humans house a community of bacteria that are significantly different in regard to diversity, richness, and abundance on human hands. We test this hypothesis utilizing computational amplicon sequencing analyses that parse data particular to the bacterial communities found within Fierer et al (2010).

Methods
=======

Sample origin and sequencing
----------------------------

### Sample collection

All analysis was done using Fierer et al. as a reference source for sequence data and sample origins. Within the paper it is stated that samples were taken from nine healthy adults (four female, five male) who worked within the same building at the University of Colorado between ages 18-40. The samples taken included computer mouse swabs touched by the owner 12h before swabbing and the palm of the individual's dominant hand used to operate the mouse. Individuals were required to maintain normal hygiene habits prior to the swab as to not increase variation. The swabs were taken using autoclave cotton-tipped swabs that were autoclave and pre-moistened with a sterile solution. The swabbing took place for 10s per sample taken with all samples being kept at -80C for less than a week before DNA extraction (Fierer *et al.*, 2010).

### DNA and Pyrosequencing

Fierer et. al 2010 explains that for each sample the 16s rRNA genes were amplified using the MO BIO PowerSoil DNA Isolation kit with the broken, frozen cotton swabs. These tubes were then horizontally shaken and kit procedures were followed for extraction. PCR reactions were carried out in triplicate repeats using HotMAsterMic with thermal cycling at 94C for 3min followed by 35 cycles of denaturation, annealing for 30s, and extension and final extension procedures. Replicate amplicons were pooled in agarose gel using SYBR Safe DNA gel stain from Invitrogen. Amplicon DNA concentrations were then measured with the final pool of DNA being precipitated, centrifuged, and centrifuged to create a pellet that was re-suspended in nuclease-free water. Pyrosequencing was then carried out on a 454 Life Sciences Genome Sequencer FLX instrument at the University of South Carolina (an NGS system used to sequence DNA samples) (Fierer *et al.*, 2010).

Computational
-------------

To create a meta data set with all of the vectors of interest we initially utilized the raw data sets collected from the 454 Sequencer and published by Fierer et al. This was done using a Dada pipeline within R. We ordered the samples first, then extracted the sample names from their fastq format, which was initially done using a remote server and fastq-dump to download the list of files in the run table to the raw directory. QC reports were then created for each of the runs utilizing fastqc and outputted as HTML to be readable. We then plotted the quality of each of the twenty samples of interest into a readable format in order to deduce the length to trim. It was found that quality is reduced after 200 bases, so the maximum acceptable length was made to be 225 bases. Using the dada pipeline the sequences were filtered, trimmed, into a new output folder and allowed to have 3 expected errors ((McMurdie and Holmes, 2013). A table of read counts was formatted to visualize the reads before and after filters and then again to visualize error trends. Duplicated sequences were then removed and dada was run on the reads based on 454 data recommendations. The sequences were aligned to craft a site by species matrix and a histogram representation of sequence lengths from all samples. Chimeras were removed and a singular table to give all pipeline information of the sequences trimmed and edited was created. A taxa code was initiated to yield a table with the taxonomy of each individual sequence to create a phylogeny that expresses the overall relatedness of each sequence. All of these tables were crafted in the dada pipeline through the sex of both male and female sample sets to visualize relatedness among all samples tested.

Once the data set was compiled by relatedness phyloseq organized all aspects of the data into a merged meta data set. This data set was parsed to remove any non-applicable data regarding the sex of participants (i.e. all samples which included swabs from electronic devices). The data was then melted to include all taxonomically related data sets that were coded for under the taxa file, separate from the metadata that was read in. This allowed for analysis via tables, figures, and ggplot graphs.

Dplyr and ggplot packages were used to analyze the data through representative figures and tables for the melted figure which was a compilation of the taxonomical data and pyrosequenced data. The pyrosequenced data set was also pruned in order to use pyroseq plot format on samples which include only males and females. The figures used look at abundance, richness, diversity metrics (applied through Shannon diversity), a phylogeny of the entire community set, and an ordination all of which are separated by the respective sex to divulge sex based community analysis on the bacterial species present.

Results
=======

Parsing Data
------------

To divulge significant evidence towards our hypothesis on the differences between sexually unique bacterial communities computational metrics were applied to the 2010 Fierer et al study in order to parse the data to uniquely address male and female samples alone.

The primary step for this analysis included analyzing the 454 sequenced samples to be trimmed in order to reduce errors and then compiled into an easy to read spreadsheet where sequence identifiers, sequence length, read errors, taxonomy, and sex-based sample type were listed for imperative analytic usage.

The primary graphs represent each of 20 samples original sequence error margins.

``` r
# Plots the quality profiles of all twenty samples
plotQualityProfile(filenames_forward_reads[1:20])
```

![](Analysis_Report_01_amplicons_files/figure-markdown_github-ascii_identifiers/check-quality-plots-1.png)

We can see from the quality profiles that most reads tend to get pretty bad in quality after around 200 bases. Therefore, we decided to set a maximum acceptable sequence length of 225 bases.

|                  |  Reads In|  Reads Out|
|------------------|---------:|----------:|
| ERR1942280.fastq |       404|        350|
| ERR1942281.fastq |       422|        194|
| ERR1942282.fastq |       412|         31|
| ERR1942283.fastq |       791|        426|
| ERR1942284.fastq |       677|        525|
| ERR1942285.fastq |       443|         72|
| ERR1942286.fastq |       667|        617|
| ERR1942287.fastq |       590|        541|
| ERR1942288.fastq |       908|        877|
| ERR1942289.fastq |       372|        147|
| ERR1942290.fastq |       468|        249|
| ERR1942291.fastq |       933|        819|
| ERR1942292.fastq |       724|        709|
| ERR1942293.fastq |       811|        470|
| ERR1942294.fastq |       938|        552|
| ERR1942295.fastq |       705|        620|
| ERR1942296.fastq |       754|        441|
| ERR1942297.fastq |       275|        246|
| ERR1942298.fastq |       562|        389|
| ERR1942299.fastq |      1025|        852|

    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Initializing error rates to maximum possible estimate.
    ## Sample 1 - 350 reads in 72 unique sequences.
    ## Sample 2 - 194 reads in 163 unique sequences.
    ## Sample 3 - 31 reads in 25 unique sequences.
    ## Sample 4 - 426 reads in 176 unique sequences.
    ## Sample 5 - 525 reads in 134 unique sequences.
    ## Sample 6 - 72 reads in 65 unique sequences.
    ## Sample 7 - 617 reads in 178 unique sequences.
    ## Sample 8 - 541 reads in 135 unique sequences.
    ## Sample 9 - 877 reads in 201 unique sequences.
    ## Sample 10 - 147 reads in 107 unique sequences.
    ## Sample 11 - 249 reads in 181 unique sequences.
    ## Sample 12 - 819 reads in 212 unique sequences.
    ## Sample 13 - 709 reads in 128 unique sequences.
    ## Sample 14 - 470 reads in 171 unique sequences.
    ## Sample 15 - 552 reads in 250 unique sequences.
    ## Sample 16 - 620 reads in 141 unique sequences.
    ## Sample 17 - 441 reads in 186 unique sequences.
    ## Sample 18 - 246 reads in 88 unique sequences.
    ## Sample 19 - 389 reads in 332 unique sequences.
    ## Sample 20 - 852 reads in 239 unique sequences.
    ##    selfConsist step 2 
    ##    selfConsist step 3 
    ## 
    ## 
    ## Convergence after  3  rounds.
    ## Total reads used:  9127

This table and corresponding error models were built to showcase the read errors from the sequencing data and the corresponding trimming profile and alignment that needed to be accounted for to allow for reproducible and unbiased analysis (removing repeats, trimming, accounting for error bias, and chimeric regions with all sequences being parsed with 454 acceptable protocols through a dada pipeline).

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](Analysis_Report_01_amplicons_files/figure-markdown_github-ascii_identifiers/visualize-errors-with-plots-1.png)

    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.

    ## Sample 1 - 350 reads in 72 unique sequences.
    ## Sample 2 - 194 reads in 163 unique sequences.
    ## Sample 3 - 31 reads in 25 unique sequences.
    ## Sample 4 - 426 reads in 176 unique sequences.
    ## Sample 5 - 525 reads in 134 unique sequences.
    ## Sample 6 - 72 reads in 65 unique sequences.
    ## Sample 7 - 617 reads in 178 unique sequences.
    ## Sample 8 - 541 reads in 135 unique sequences.
    ## Sample 9 - 877 reads in 201 unique sequences.
    ## Sample 10 - 147 reads in 107 unique sequences.
    ## Sample 11 - 249 reads in 181 unique sequences.
    ## Sample 12 - 819 reads in 212 unique sequences.
    ## Sample 13 - 709 reads in 128 unique sequences.
    ## Sample 14 - 470 reads in 171 unique sequences.
    ## Sample 15 - 552 reads in 250 unique sequences.
    ## Sample 16 - 620 reads in 141 unique sequences.
    ## Sample 17 - 441 reads in 186 unique sequences.
    ## Sample 18 - 246 reads in 88 unique sequences.
    ## Sample 19 - 389 reads in 332 unique sequences.
    ## Sample 20 - 852 reads in 239 unique sequences.

    ## The sequences being tabled vary in length.

**Figure 1**: Sequence Length Histogram ![](Analysis_Report_01_amplicons_files/figure-markdown_github-ascii_identifiers/histogram-of-sequence-lengths-1.png)

This table explicitly represents the sequence lengths of the trimmed data set to show they are all correctly analyzed beneath the given range of 225 bases.

After removing chimeras, we were left with 99.65% of our cleaned reads.

|            |  Input|  Filtered|  Denoised|  Sequence Table|  Non-chimeric|
|------------|------:|---------:|---------:|---------------:|-------------:|
| ERR1942280 |    404|       350|       350|             350|           350|
| ERR1942281 |    422|       194|       194|             194|           194|
| ERR1942282 |    412|        31|        31|              31|            31|
| ERR1942283 |    791|       426|       426|             426|           426|
| ERR1942284 |    677|       525|       525|             525|           525|
| ERR1942285 |    443|        72|        72|              72|            72|
| ERR1942286 |    667|       617|       617|             617|           585|
| ERR1942287 |    590|       541|       541|             541|           541|
| ERR1942288 |    908|       877|       877|             877|           877|
| ERR1942289 |    372|       147|       147|             147|           147|
| ERR1942290 |    468|       249|       249|             249|           249|
| ERR1942291 |    933|       819|       819|             819|           819|
| ERR1942292 |    724|       709|       709|             709|           709|
| ERR1942293 |    811|       470|       470|             470|           470|
| ERR1942294 |    938|       552|       552|             552|           552|
| ERR1942295 |    705|       620|       620|             620|           620|
| ERR1942296 |    754|       441|       441|             441|           441|
| ERR1942297 |    275|       246|       246|             246|           246|
| ERR1942298 |    562|       389|       389|             389|           389|
| ERR1942299 |   1025|       852|       852|             852|           852|

This markdown table showcases the edited data and the reads given out after each step of filtering for the corrected data set to be used in figure based computational analysis.

``` r
# show the results of the taxonomy assignment
unname(taxa)
```

Taxonomy was assigned for each sequence to describe which bacterial community was present within the samples tested for.

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 4.418541e-05 
    ## Run 1 stress 0 
    ## ... New best solution
    ## ... Procrustes: rmse 0.2503605  max resid 0.3937541 
    ## Run 2 stress 9.4848e-05 
    ## ... Procrustes: rmse 0.274344  max resid 0.4919324 
    ## Run 3 stress 0 
    ## ... Procrustes: rmse 0.2467099  max resid 0.3615862 
    ## Run 4 stress 6.496638e-05 
    ## ... Procrustes: rmse 0.2715979  max resid 0.3993061 
    ## Run 5 stress 3.408458e-05 
    ## ... Procrustes: rmse 0.2661741  max resid 0.4694951 
    ## Run 6 stress 0 
    ## ... Procrustes: rmse 0.1874121  max resid 0.4149212 
    ## Run 7 stress 8.56982e-05 
    ## ... Procrustes: rmse 0.3092958  max resid 0.5516642 
    ## Run 8 stress 9.388572e-05 
    ## ... Procrustes: rmse 0.2254944  max resid 0.3150192 
    ## Run 9 stress 4.835216e-05 
    ## ... Procrustes: rmse 0.2571429  max resid 0.4897982 
    ## Run 10 stress 0 
    ## ... Procrustes: rmse 0.17153  max resid 0.3092111 
    ## Run 11 stress 0 
    ## ... Procrustes: rmse 0.2727335  max resid 0.5467948 
    ## Run 12 stress 9.477793e-05 
    ## ... Procrustes: rmse 0.2812978  max resid 0.6012645 
    ## Run 13 stress 0 
    ## ... Procrustes: rmse 0.2969349  max resid 0.5152367 
    ## Run 14 stress 5.968682e-05 
    ## ... Procrustes: rmse 0.3078262  max resid 0.512252 
    ## Run 15 stress 0 
    ## ... Procrustes: rmse 0.2872692  max resid 0.4966087 
    ## Run 16 stress 9.764269e-05 
    ## ... Procrustes: rmse 0.2680452  max resid 0.4627489 
    ## Run 17 stress 9.766709e-05 
    ## ... Procrustes: rmse 0.2450739  max resid 0.4410873 
    ## Run 18 stress 7.244177e-05 
    ## ... Procrustes: rmse 0.2741797  max resid 0.4901149 
    ## Run 19 stress 8.827905e-05 
    ## ... Procrustes: rmse 0.2997358  max resid 0.4773969 
    ## Run 20 stress 0 
    ## ... Procrustes: rmse 0.2970872  max resid 0.4142021 
    ## *** No convergence -- monoMDS stopping criteria:
    ##     20: stress < smin

    ## Warning in metaMDS(veganifyOTU(physeq), distance, ...): Stress is (nearly)
    ## zero - you may have insufficient data

    ## Warning in postMDS(out$points, dis, plot = max(0, plot - 1), ...): skipping
    ## half-change scaling: too few points below threshold

Phyloseq is used to create a compiled data set with sequenced MetaData and taxonomy. This is parsed to remove Non-applicable sex based data (i.e. from electronic device swabs.). The data is also melted together to form a larger data set applicable for computational analyses of sexual data characteristics and taxonomical characteristics.

After parsing and filtering, the "melted" data set is now applicable for analytical research regarding the community analysis of male and female participants, and their corresponding bacterial communities.

Analyitics
----------

**Figure 2**: Abundance Measure for Male versus Female sample types

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1660 rows containing non-finite values (stat_boxplot).

![](Analysis_Report_01_amplicons_files/figure-markdown_github-ascii_identifiers/Abundance-1.png)

**Figure 3**: Richness Measure for Male versus Female sample types

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](Analysis_Report_01_amplicons_files/figure-markdown_github-ascii_identifiers/Richness-plot-1.png)

**Figure 4**: Alpha diversity measures of the two sample types

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](Analysis_Report_01_amplicons_files/figure-markdown_github-ascii_identifiers/Shannon%20Diversity%20Plot-1.png)

**Figure 5**: Male versus Female bar graph showcasing the phylum present within each bacterial community ![](Analysis_Report_01_amplicons_files/figure-markdown_github-ascii_identifiers/Abundance%20Taxonomy%20Graph-1.png)

**Figure 6**: Inferred phylogeny of sequences, with points on tips representing samples within which each particular taxa occurred. Tree represents maximum likelihood phylogeny inferred using RAxML.

    ## Warning: Removed 1 rows containing missing values (geom_segment).

![](Analysis_Report_01_amplicons_files/figure-markdown_github-ascii_identifiers/example-phyloseq-plot-1.png)

**Figure 7**: Ordination to order species along presumed ecological gradient to identify patterns of species distribution and abundance for male versus female samples

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 4.418541e-05 
    ## Run 1 stress 0 
    ## ... New best solution
    ## ... Procrustes: rmse 0.2577542  max resid 0.4821088 
    ## Run 2 stress 0 
    ## ... Procrustes: rmse 0.243252  max resid 0.4383912 
    ## Run 3 stress 0 
    ## ... Procrustes: rmse 0.2660546  max resid 0.5502146 
    ## Run 4 stress 0 
    ## ... Procrustes: rmse 0.2981806  max resid 0.4718184 
    ## Run 5 stress 0 
    ## ... Procrustes: rmse 0.2404484  max resid 0.4947174 
    ## Run 6 stress 0 
    ## ... Procrustes: rmse 0.2556273  max resid 0.5807703 
    ## Run 7 stress 0 
    ## ... Procrustes: rmse 0.17041  max resid 0.2740696 
    ## Run 8 stress 4.483392e-05 
    ## ... Procrustes: rmse 0.2574179  max resid 0.3840677 
    ## Run 9 stress 0 
    ## ... Procrustes: rmse 0.2560501  max resid 0.4840124 
    ## Run 10 stress 0 
    ## ... Procrustes: rmse 0.1789803  max resid 0.3300431 
    ## Run 11 stress 2.896177e-06 
    ## ... Procrustes: rmse 0.2721875  max resid 0.5789352 
    ## Run 12 stress 3.438036e-05 
    ## ... Procrustes: rmse 0.2769226  max resid 0.6033152 
    ## Run 13 stress 7.496779e-05 
    ## ... Procrustes: rmse 0.287852  max resid 0.5677947 
    ## Run 14 stress 9.759861e-05 
    ## ... Procrustes: rmse 0.2402742  max resid 0.4537296 
    ## Run 15 stress 0 
    ## ... Procrustes: rmse 0.255512  max resid 0.4982651 
    ## Run 16 stress 9.935633e-05 
    ## ... Procrustes: rmse 0.2958957  max resid 0.6594391 
    ## Run 17 stress 7.021784e-05 
    ## ... Procrustes: rmse 0.2122595  max resid 0.4071978 
    ## Run 18 stress 9.785171e-05 
    ## ... Procrustes: rmse 0.2138188  max resid 0.4404271 
    ## Run 19 stress 1.530938e-05 
    ## ... Procrustes: rmse 0.282998  max resid 0.5321475 
    ## Run 20 stress 7.746608e-05 
    ## ... Procrustes: rmse 0.3019012  max resid 0.5618569 
    ## *** No convergence -- monoMDS stopping criteria:
    ##     20: stress < smin

    ## Warning in metaMDS(veganifyOTU(physeq), distance, ...): Stress is (nearly)
    ## zero - you may have insufficient data

    ## Warning in postMDS(out$points, dis, plot = max(0, plot - 1), ...): skipping
    ## half-change scaling: too few points below threshold

![](Analysis_Report_01_amplicons_files/figure-markdown_github-ascii_identifiers/Ordination%20Male%20vs.%20Female-1.png)

Discussion
==========

Understanding the microbiota which live within our unique niche environments can yield information regarding medical advancements. Specifically, understanding how these microbiota differ within female and male populations can increase our understanding of hormone differences and gland production that can yield evidence towards how these microbes symbiotically interact with human host systems.

Here we prove that differences between male and female bacterial communities within the hand are present, in regards to; abundance, diversity, and richness. We further our analysis to divulge if the similar environments constitute a similar makeup of bacteria present through taxonomical representation.

### Computational

Computational analyses include preparing the sequencing data in order to have trimmed and aligned sequences that appropriately underscore the important information needed to be extrapolated for further analyses. The primary quality plot data showcases the region of the sequences that needs to be trimmed to account for sequence errors; primarily after 220 base pairs. The corresponding errors graphic, after trimming, shows the statistical margin of error that occurs when longer sequences are read through an NGS system like 454. The error graph shows a direct relationship between the number of errors and the lengths of the sequences. Once trimmed, Figure 1 portrays the new sequence lengths using a histogram to graphically underscore that all sequences are now under the 225 base pair limit. The next markdown table showcases how many reads were present before and after trimming, alignment, removing chimeras and removing duplicates to again underscore how the data was manipulated, in order to, increase statistical rational towards our hypothesis and decrease extraneous error in tests we wanted to run. Lastly, a taxonomy table that is compiled with ribosomal database information highlights the comparisson of sequence data to corresponding bacterial types.

### Abundance

Abundance is a measure of the total number of individuals residing within the same niche environment to statistically observe how consecrated the region is. If a high abundance occurs in one environment it can be observed as having more resources, or less environmentally de-stabling associations. Abundance is the relative representation of a species in a particular ecological niche, which in this case is the human hand skin for both male and female participants.

Evaluation of Figure 2 expresses an abundance of differing amounts between the sample sexes with male samples having approximately 3000 individual bacterial colony reads. On female hand samples there was only approximately 1500 colony reads, meaning male hands carried twice as much bacterial on their hands as females. This means male participants overall have a larger bacterial community residing on their hands than females.

However, when abundance levels are overlapped with taxonomical representation there is a similarity between the two bacterial communities represented; both communities have more proteobacteria than any other bacteria (Figure 5). We can conclude that the similar environmental pressures faced by both sample sets, due to living in Colorado, could account for this.

### Richness

Richness is a measure of different morphotypes within a bacterial community. Calculations of richness within a community are imperative to showcase if only one species is present, or many within the same niche environment. Biodiversity is measured through both richness and diversity factors. By observing biodiversity within communities, conservation of certain species can occur through the creation of more favorable niches. In the observation between male and female resident sample populations, like richness is observed in the overlap of some varying morphotypes. However, due to the ability for potential mutation or nutritional needs of the bacteria, some morphotypes favored the males over the females environment and vice versa.

Evaluation of Figure 3 showcases the differences within richness metrics between both female and male bacterial communities. Unlike abundance, where males showcased more colony reads, we visualize female samples to have an increased level of richness to deduce that a large variety of bacterial species are present on the female hand. Though male samples have some diversity, female samples show a richness level double that of male samples (9:4). We can conclude that male samples then have more colonies of a singular bacterial type, where female samples have fewer overall colonies but of more diverse types.

### Diversity

Diversity is a measure of the relative variability among organisms from all sources or morphologies. Shannon Diversity offers a statistical means of categorizing the diversity of individuals within a community by accounting for abundance and evenness. We apply Shannon Diversity as a metric to compare the male and female bacterial communities farther.

Within Figure 4 we see an increased alpha level for female samples than that of males, which is comparable to the richness measure given previously in Figure 3. This measure is almost double that of the male samples. We can conclude that for the region in which the bacteria inhabit (human hand), the female communities offer a more diverse community. High diversity is favored over low diversity since this ensures natural sustainability to the micro-community which populates the hand.

#### Phylogeny

Utilizing the Ribosomal Database Project (RDP via DADA2) we created a taxonomy based on the sequences found. This taxonomy allowed for a direct comparison of the bacterial communities found, and their sex-based origin.

Within Figure 6, we visualize this phylogeny that is colored over to represent male and female differences within each branch. We visualize more dispersion over all of the branches for female samples, showing multiple phylums, genus, and species all stemming from multiple ancestor lineages. We can conclude that this spread within the phylogeny directly compares to the level of richness represented in prior figures. However in comparison, male sample data, which has a decreased spread, still has taxonomies that overlap with female samples. This proves that though many of the microbiota are different between male and female samples, overall most communities present in male samples are also simultaneously presented in female samples. This means that male populations are more alike to female populations than female populations are to males.

### Ordination

Ordination can be used as a exploratory measure for data clustering of similar and different multivariate types. When looking at both male and female sample types we visualize more clustering in the male sample types, concluding that male populations show a higher level of relatedness over all tested variables (Figure 7). This relatedness counters the larger levels of dispersion seen within female samples.

Plotting an ordination is a measure for quantifying richness and diversity simultaneously over multiple axes, and here we visualize these traits are more closely related for male species. However, when looking at the male samples alone we do see some dispersion between the two predominant clusters, which we can deduce as two more diverse populations. Potentially, these clusters represent species found on a different region of the hand where there are more glands, more contact with the environment, or constitute more rough or raw skin. When looking at female samples, the populations are not closely clustered at all, and show a dramatic spread over the plot to accentuate the diversity between all groups found. Potentially, these groups have similar nutritional needs, but function differently within the skin microbiota population.

### Conclusion

Overall, we utilize measures of abundance, richness, and diversity to prove male and female hands have different bacterial communities. We see that males house a community of more abundant bacteria that are more closely related overall, while females have less bacteria but of more diverse morphologies. We also explain that though the communities are overall different they do represent similarities in the phylum of bacteria seen for both sex samples (Proteobacteria), which is most likely from the overall similar environment the study participants live in (University of Colorado at Boulder). We can postulate that the microbiota on the male hand are less diverse because of the skin pH, which is nominally different between the two sexes due to sweat glands which produce a more acidic environment unfavorable to many microbiota and are stimulated by hormonal differences between the sexes {Fierer *et al.* (2008)}. These conditions would also account for differences observed in abundance and richness. Testosterone production can effect sweat glands on the human skin, which are known for maintaining and regulated body temperature and expelling toxins to support immune responses (KUWAHARA *et al.*, 2006).

Further studies need to be conducted to directly prove this association with pH, and simultaneously disprove these conclusions for the different bacterial communities are not from hygiene differences or social constructs which imply females using more cosmetic products. More information also needs to be found regarding the production of testosterone and it's correlation to sweat gland activity versus estrogen. By understanding these differences between sexes, we can develop a better understanding of how these differences can be responsible for changes in the health and well-being of individuals studied.

Sources Cited
=============

Birtel,J. *et al.* (2015) Estimating bacterial diversity for ecological studies: Methods, metrics, and assumptions. *PloS one*, **10**, e0125356.

Fierer,N. *et al.* (2008) The influence of sex, handedness, and washing on the diversity of hand surface bacteria. *Proceedings of the National Academy of Sciences*, **105**, 17994–17999.

Fierer,N. *et al.* (2010) Forensic identification using skin bacterial communities. *Proceedings of the National Academy of Sciences*, **107**, 6477–6481.

Fuhrman,J.A. (1999) Marine viruses and their biogeochemical and ecological effects. *Nature*, **399**, 541–548.

KUWAHARA,T. *et al.* (2006) Sex differences in effects of physical training on sweat gland function (proceedings of the 54th meeting of japan society of physiological anthropology). *Journal of physiological anthropology*, **25**, 200.

Markle,J.G. *et al.* (2013) Sex differences in the gut microbiome drive hormone-dependent regulation of autoimmunity. *Science*, **339**, 1084–1088.

McMurdie,P.J. and Holmes,S. (2013) Phyloseq: An r package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE*, **8**, e61217.

Nichols,D. (2007) Cultivation gives context to the microbial ecologist. *FEMS microbiology ecology*, **60**, 351–357.

Wu,S. *et al.* (2012) Composition, diversity, and origin of the bacterial community in grass carp intestine. *PloS one*, **7**, e30440.
