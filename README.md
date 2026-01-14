# Code for phylogenetic analysis of passerine songs

Tree and data from: https://github.com/quentinbacquele/phylogenetic_analysis




Resources:
- _**threshml_fit_2_corr_matrix.csv**_: estimated covariance matrix between underlying liabilities from multivariate threshold model for motif presence. Extracted from ThreshML run.
- _**namefix_trait_data_100_1912.csv**_: trait data, with names aligned with tree via matching synonyms.
- _**namefix_tree_100_1912.tre**_: phylogenetic consensus tree.

Scripts:
- _**univariate_threshold_multiple_chain.R**_: script to estimate motif presence/absence for ancestral species.
- _**plot_multichain_threshold.R**_: plotting motif presence/absence reconstructions.
- _**clade_PCA_reconstruction.R**_: reconstructing acoustic features for motifs within chosen clades.
- _**multivariate_threshold_model.R**_: visualising and analysing motif correlation matrix.


