library(clusterProfiler)
library(org.Hs.eg.db)

##### GO enrichment
enrich.go <- enrichGO(gene = gene_list,
                      universe = background_list,
                      OrgDb = 'org.Hs.eg.db',
                      ont = 'ALL',
                      pAdjustMethod = 'fdr')

##### KEGG enrichment
enrich.kegg <- enrichKEGG(gene = gene_list,
                          universe = background_list,
                          organism = 'hsa',
                          keyType = 'kegg',
                          pAdjustMethod = 'fdr')


