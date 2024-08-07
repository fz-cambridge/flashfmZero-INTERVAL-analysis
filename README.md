# flashfmZero-INTERVAL-analysis
<!-- badges: start -->
[![DOI]
<!-- badges: end -->
The flashfmZero-INTERVAL factor analysis and fine mapping (scripts and instructions)

## Replication

This repository holds the scripts and instructions used to generate the results for the flashfmZero-INTERVAL factor analysis and fine mapping paper (cited below). For the repository containing the source code of flashfmZero R package, see: https://github.com/jennasimit/flashfmZero

Paper citation:

> *Improved fine-mapping resolution of high-dimensional traits through multivariate analyses of latent factors* <br />
> F Zhou, WJ Astle, AS Butterworth, JL Asimit <br />
> The preprint and its supplemental material can be found in *bioRxiv* <br />
> doi: https://...


## Background
Genome-wide association studies (GWAS) of high-dimensional traits, such as molecular phenotypes or imaging parameters, are becoming increasingly common. Often such traits are biologically related and have shared genetic association signals. Factor analysis provides a way to estimate a smaller number of latent factors underlying many original traits. GWAS of latent factors are more scalable for multi-trait approaches and can capture biological mechanisms generating variation in high-dimensional traits parsimoniously. Here, we introduce a zero-correlation multi-trait fine-mapping approach, flashfmZero, for any number of latent factors. In our application to 25 statistically uncorrelated, yet biologically related latent factors derived from 99 blood cell traits in the INTERVAL cohort of UK blood donors[^1][^2][^3] and cross-checking the results of our analysis with the UK Biobank fine-mapping results[^4], we show how GWAS of latent factors enables detection of signals that have sub-threshold associations with several blood cell traits. Fine-mapping of latent factors  reduced the size of credible sets compared to blood cell traits. These analysis techniques ease interpretation of results from many traits and highlight common underlying factors amongst them.

The *R_script_INTERVAL.R:* describes the whole process of the main analysis:
- Step_1: collect and combine the 99 blood cell traits from [^2][^3];
- Step_2: delete the rows of sample with missing values (so the sample size is reduced to 18k);
- Step_3: (important) the size-reduced 99 blood cell traits are normalised again;
- Step_4: run factor analyses (FA) based on size-reduced 99 normalised blood cell traits; 
- Step_5: get the optimal 25 FA latent factors with a matrix contains all loadings of 99 blood cell traits;
- Step_6: (optional) the 25 FA latent factors are also normalised but the effect/impact is small;
- Step_7: link all 99 blood cell traits with the 25 FA latent factors and create a network/connection visualization;
- Step_8: re-use BOLT-LMM GWAS of full-size 99 blood cell traits from the two published papers[^2][^3];
- Step_9: run BOLT-LMM to get GWAS of reduced-size 99 normalised blood cell and 25 latent factors; 
- Step_10: compare GWAS signals between different traits/factors and methods (i.e. blood cell traits vs latent factors);
- Step_11: check the (combination of) strengths and regions of GWAS signals;
- Step_12: check correlations of genotypes between GWAS significant SNPs;
- Step_13: implement conditional analyses by adding lead SNPs (from the two published papers[^2][^3]) into LMM regressions;
- Step_14: detect novel SNPs based on FA in comparison with blood cell traits (both full-size sample 43k and the reduced-size sample 18k);
- Step_15: use VEP (build 37) and connections of latent factors with blood cell traits to better understand novelty;
- Step_16: fine-mapping based on both FA latent factors and blood cell traits, compare the results.

## Reference
[^1]: https://www.donorhealth-btru.nihr.ac.uk/studies/interval-study/
[^2]: Astle, W.J., Elding, H., Jiang, T., Allen, D., Ruklisa, D., Mann, A.L., Mead, D., Bouman, H., Riveros-Mckay, F., Kostadima, M.A., et al. (2016). The Allelic Landscape of Human Blood Cell Trait Variation and Links to Common Complex Disease. Cell 167, 1415–1429.e19.
[^3]: Akbari, P., Vuckovic, D., Stefanucci, L., Jiang, T., Kundu, K., Kreuzhuber, R., Bao, E.L., Collins, J.H., Downes, K., Grassi, L., et al. (2023). A genome-wide association study of blood cell morphology identifies cellular proteins implicated in disease aetiology. Nat. Commun. 14, 5023.
[^4]: Vuckovic, D., Bao, E.L., Akbari, P., Lareau, C.A., Mousas, A., Jiang, T., Chen, M.-H., Raffield, L.M., Tardaguila, M., Huffman, J.E., et al. (2020). The Polygenic and Monogenic Basis of Blood Traits and Diseases. Cell 182, 1214–1231.e11.



