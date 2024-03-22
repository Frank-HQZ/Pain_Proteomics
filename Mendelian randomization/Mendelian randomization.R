library(TwoSampleMR)

##### Read exposure data
exposure_data <- read.table('exposure_data.txt')
exposure_data <- format_data(exposure_data,type='exposure',
                             snp_col = "rsid",beta_col = "BETA",
                             se_col = "SE",effect_allele_col = "A1",
                             other_allele_col = "A0",
                             pval_col = "P",
                             chr_col = "CHROM")

##### Read outcome data
outcome_data <- read.table('outcome_data.txt')
outcome_data <- format_data(outcome_data,type='outcome',
                            snp_col = "rsid",beta_col = "BETA",
                            se_col = "SE",effect_allele_col = "A1",
                            other_allele_col = "A0",
                            pval_col = "P",
                            chr_col = "CHROM")

##### Clumped exposure data
snp_clump <- ld_clump_local( dplyr::tibble ( rsid = exposure_data$SNP, pval = exposure_data$pval.exposure),
                             clump_kb = 1000, clump_r2 = 0.01, clump_p = 1, 
                             bfile = 'g1000_eur', plink_bin = "plink")
exposure_clump <- exposure_data[which(exposure_data$SNP %in% snp_clump$rsid),]

##### Perform MR analysis
dat <- harmonise_data(exposure_dat = exposure_clump, outcome_dat = outcome_data)
res <- mr(dat)
  
##### Test heterogeneity
heterogeneity <- mr_heterogeneity(dat)
mre <- mr(dat,method_list=c('mr_ivw_mre'))
fe <- mr(dat,method_list=c('mr_ivw_fe'))

##### Test pleiotropy
pleiotropy <- mr_pleiotropy_test(dat)
