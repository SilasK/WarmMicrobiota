
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        scriptdir = "character"
    )
)
snakemake <- Snakemake(
    input = list('output/seqtab_with_chimeras.rds', "seqtab" = 'output/seqtab_with_chimeras.rds'),
    output = list('output/seqtab_nochimeras.rds', 'stats/Nreads_chimera_removed.txt', "rds" = 'output/seqtab_nochimeras.rds', "nreads" = 'stats/Nreads_chimera_removed.txt'),
    params = list(),
    wildcards = list(),
    threads = 1,
    log = list('logs/dada2/removeChimeras.txt'),
    resources = list(),
    config = list("java_mem" = 30, "threads" = 4, "truncLen" = c(180, 100), "maxEE" = c(2, 2), "truncQ" = 2, "learn_nbases" = '100e6', "chimera_method" = 'consensus', "idtaxa_dbs" = list("Silva" = '/Users/silas/Desktop/WarmMicrobiota/TaxonmyDBs/SILVA_SSU_r132_March2018.RData', "GTDB" = '/Users/silas/Desktop/WarmMicrobiota/TaxonmyDBs/GTDB_r86-mod_September2018.RData'), "sampletable" = '/Users/silas/Desktop/WarmMicrobiota/Sequencing_C/samples_WartmTransplanted.tsv'),
    rule = 'removeChimeras',
    scriptdir = '/Users/silas/Documents/GitHub/amplicon-seq-dada2/scripts/dada2'
)
######## Original script #########
library(dada2)

sink(snakemake@log[[1]])

seqtab.all= readRDS(snakemake@input[['seqtab']]) # seqtab.all


# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab.all, method=snakemake@config[['chimera_method']], multithread=snakemake@threads)

saveRDS(seqtab, snakemake@output[['rds']]) 


track <- rowSums(seqtab)
names(track) <- colnames(seqtab)

write.table(track,col.names = c("nonchim"), 
            snakemake@output[['nreads']],sep='\t')








