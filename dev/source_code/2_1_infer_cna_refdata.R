#setting up the genome information
creat_marts = function(genome='hg38', host='ensembl.org') {
    if (genome == 'hg19') {
        dataset = 'hsapiens_gene_ensembl'
        #host = 'grch37.ensembl.org'
    } else if (genome == 'hg38') {
        dataset = 'hsapiens_gene_ensembl'
        #host = 'useast.ensembl.org'
        #host = 'uswest.ensembl.org'
        #host = 'asia.ensembl.org'
    } else if (genome == 'mm10') {
        dataset = 'mmusculus_gene_ensembl'
    } else {
        stop('genome "', genome, '" not found')
    }
    # Choose which species to use and server to download from
    biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
             dataset = dataset,
             host = host,
             port = 80)
}

get_genomes = function(mart, add_attribute='hgnc_symbol',
                       genes.dir='/scratch/splab/zqi/htan/code/data/genes.rda') {
    filters = c('chromosome_name', 'hgnc_symbol')
    values = list(c(1:22, 'X', 'Y'), get(load(genes.dir)))
    attributes = c('chromosome_name',
                   'band',
                   'strand',
                   'start_position',
                   'end_position')
    attributes = unique(c(attributes, add_attribute))
    result = biomaRt::getBM(attributes = attributes,
                            mart = mart,
                            filters = filters,
                            values = values,
                            uniqueRows = T)
    if ('hgnc_symbol' %in% add_attribute) {
        result = dplyr::filter(result, hgnc_symbol != "")}
    if ('mgi_symbol' %in% add_attribute) {
        result = dplyr::filter(result, mgi_symbol != "")}
    return(result)
}

set_chr_arm_levels = function(dat) {
    x = unique(dat$chromosome_name)
    xnum = sort(as.numeric(x[x %in% as.character(1:100)]))
    xchar = sort(x[!x %in% xnum]) # sort(x[!x %in% as.character(1:100)])
    chrlev = c(xnum, xchar)
    dat = dplyr::mutate(dat, chromosome_name = factor(as.character(chromosome_name), levels = chrlev))
    dat = dplyr::arrange(dat, chromosome_name, start_position)
    armlev = unique(dat$arm)
    dat = dplyr::mutate(dat, arm = factor(as.character(arm), levels = armlev))
    return(dat)
}

###retrive hg38 mart and covert it to df by getBM==============================
#https://rdrr.io/github/jlaffy/infercna/src/data-raw/set_sysdata.R
#https://rdrr.io/github/jlaffy/infercna/src/data-raw/retrieve_genome_info.R
hg19 = creat_marts('hg19')
hg38 = creat_marts('hg38')
mm10 = creat_marts('mm10')
marts = list(hg19, hg38, mm10)
add_attribute = c('hgnc_symbol', 'hgnc_symbol', 'mgi_symbol')
#make sure the order matches in marts and add_attribute
dats = Map(get_genomes, mart = marts, add_attribute = add_attribute)
#dats is a list of data.frames. No list names.
#sapply(1:3, function(i) colnames(dats[[i]])[7] <<- 'symbol')
sapply(1:3, function(i) colnames(dats[[i]])[6] <<- 'symbol')
#add a new col: arm
dats[1:2] = sapply(dats[1:2], 
                   function(d) 
                   dplyr::mutate(d, arm = paste0(chromosome_name, stringr::str_extract(band, "p|q"))),
                   simplify = F)
dats[[3]] = dplyr::mutate(dats[[3]], arm = chromosome_name)
#dats[[3]] = dats[[3]] %>% dplyr::mutate(arm = chromosome_name)

dats2 = sapply(dats, function(dat) set_chr_arm_levels(dat), simplify = F)
hg19 = dats2[[1]]
hg38 = dats2[[2]]
mm10 = dats2[[3]]
# usethis::use_data(hg19, hg38, mm10, internal = T)
# Error: Path '/scratch/splab/zqi/htan/' does not appear to be inside a project or package.
# In addition: Warning message:
#     In path_file(base_path) : restarting interrupted promise evaluation
rm(dats,dats2)
rm(hg19,mm10)

###set up the environmnet=====================================================
#https://rdrr.io/github/jlaffy/infercna/src/data-raw/setGenv.R
vars = as.list(as.data.frame(hg38))
nam = c('symbol', 'chromosome_name', 'arm', 'start_position', 'end_position')
nam2 = c('gene', 'chr', 'arm', 'start', 'end')
vars = stats::setNames(vars[nam], nam2)
vars[2:5] = sapply(vars[2:5], stats::setNames, vars[['gene']], simplify = F)
vars = c(list(name = 'hg38', data = data), vars)
Genv = list2env(vars)
rm(vars,nam,nam2)


