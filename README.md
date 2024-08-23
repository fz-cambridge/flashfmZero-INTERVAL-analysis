# flashfmZero-INTERVAL-analysis

The flashfmZero-INTERVAL factor analysis and fine mapping (scripts and instructions)

## Replication

This repository holds the scripts and instructions used to generate the results for the flashfmZero-INTERVAL factor analysis and fine mapping paper (cited below). For the repository containing the source code of flashfmZero R package, see: https://github.com/jennasimit/flashfmZero

Paper citation:

> *Improved genetic discovery and fine-mapping resolution of high-dimensional traits through multivariate analyses of latent factors* <br />
> F Zhou, WJ Astle, AS Butterworth, JL Asimit <br />
> The preprint and its supplemental material can be found in *bioRxiv* <br />


## Background
Genome-wide association studies (GWAS) of high-dimensional traits, such as molecular phenotypes or imaging parameters, often use univariate approaches, leaving a large carbon footprint. Biological mechanisms generating variation in high-dimensional traits can be captured parsimoniously through GWAS of a smaller number of latent factors from factor analysis. Here, we introduce a zero-correlation multi-trait fine-mapping approach, flashfmZero, for any number of latent factors. In our application to 25 latent factors derived from 99 blood cell traits in the INTERVAL cohort[^1][^2][^3], we show how GWAS of latent factors enables detection of signals that have sub-threshold associations with several blood cell traits[^4]. FlashfmZero resulted in 99% credible sets with the same size or fewer variants than those for blood cell traits in 87% of our comparisons, and all latent trait fine-mapping credible sets were subsets of those from flashfmZero. These analysis techniques give enhanced power for discovery and fine-mapping for many traits.

The *Script_INTERVAL.R:* describes the whole process of the main analysis in more detailed steps and provides the R codes for running factor analysis and fine mapping:
- Step_1: only retain individuals who have no missing trait values amongst the 99 blood cell traits[^2][^3]; the sample size is reduced to 18k;
- Step_2: use inverse normal rank transformation for each of the 99 blood cell traits in the sample of 18k;
- Step_3: run factor analysis (FA) with varimax rotation on the 99 normalised blood cell traits;
- Step_4: based on scree plot, select 25 latent factors;
- Step_5: output 25 FA latent factors estimated at each individual (factor scores) and a matrix containing all factor loadings of 99 blood cell traits;
- Step_6: use inverse normal rank transformation for each of the 25 latent factors;
- Step_7: run GWAS of each of the 25 latent traits and 99 blood cell traits in the samples of 18k using BOLT-LMM;
- Step_8: implement conditional analyses by adding lead SNPs[^2][^3] into BOLT-LMM regressions;
- Step_9: select regions for fine-mapping based on GWAS signals for both FA latent factors and blood cell traits; use JAMdynamic and flashfmZero with in-sample LD.


## Reference
[^1]: https://www.donorhealth-btru.nihr.ac.uk/studies/interval-study/
[^2]: Astle, W.J., Elding, H., Jiang, T., Allen, D., Ruklisa, D., Mann, A.L., Mead, D., Bouman, H., Riveros-Mckay, F., Kostadima, M.A., et al. (2016). The Allelic Landscape of Human Blood Cell Trait Variation and Links to Common Complex Disease. Cell 167, 1415–1429.e19.
[^3]: Akbari, P., Vuckovic, D., Stefanucci, L., Jiang, T., Kundu, K., Kreuzhuber, R., Bao, E.L., Collins, J.H., Downes, K., Grassi, L., et al. (2023). A genome-wide association study of blood cell morphology identifies cellular proteins implicated in disease aetiology. Nat. Commun. 14, 5023.
[^4]: Vuckovic, D., Bao, E.L., Akbari, P., Lareau, C.A., Mousas, A., Jiang, T., Chen, M.-H., Raffield, L.M., Tardaguila, M., Huffman, J.E., et al. (2020). The Polygenic and Monogenic Basis of Blood Traits and Diseases. Cell 182, 1214–1231.e11.



