# contains bits of code adapted from https://gist.github.com/jaquol/03f41f57dc6b0eacef101e9920f24d78

library(sleuth)
library(reshape2)

# hardcode folders that contain kallisto results
base_dir <- '..'
sample_id <- dir(file.path(base_dir, 'results'))
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, 'results', id))

print(sample_id)
print(kal_dirs)

# read experimental setup
s2c <- read.table('expt_setup.tsv', header=TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path=kal_dirs)

# check that s2c is set up correctly
print(s2c)

# be careful of nomenclature! conditions are sorted in alphabetical order, and
# the first one is always the baseline

# load data, and make a regression model with 'condition' as the dependent variable
so <- sleuth_prep(s2c, ~ condition)

# 1. print normalised TPMs
tpms <- kallisto_table(so)[, c('target_id', 'sample', 'tpm')]
tpms <- dcast(tpms, target_id ~ sample, value.var='tpm')
write.table(tpms, 'normalised_abundances.tsv', sep='\t', quote=FALSE, row.names=FALSE)

# estimate parameters for the response error measurement model that is dependent on 'condition'
so <- sleuth_fit(so)

# run another model where gene expression is independent of any factors (i.e. null expectations)
so <- sleuth_fit(so, ~1, 'reduced')

# calculate likelihood ratios that gene expression is dependent on 'condition'
so <- sleuth_lrt(so, 'reduced', 'full')

# run wald test to get beta (approx equals ln FC), even though LRT is technically better than it
# colnames(so$fits[['full']]$design_matrix) has all the possible conditions
so <- sleuth_wt(so, 'condition32C')

# 2. print DEG table
res_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
res_wt <- sleuth_results(so, 'condition32C')
res <- merge(res_wt[, c('target_id', 'b', 'se_b', 'mean_obs')], res_lrt, on='target_id', sort=TRUE)
res <- na.omit(res)
write.table(res, 'sleuth_results.tsv', sep = "\t", quote = FALSE, row.names = FALSE)
