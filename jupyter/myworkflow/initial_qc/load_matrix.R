suppressWarnings(library(SingleCellExperiment))
suppressWarnings(library(scater))

# Set the raw data and save/load it
umi <- NULL
if (!exists("generate_matrix")){
   generate_matrix = F
}

if (!exists("matrix_destination")){
   matrix_destination = "matrix.rds"
}


#genes_of_interest <- c(
#    "Eomes", "Brachyury", "Mesp1",                # meso
#    "Pou5f1", "nanog",                            # pluripotent
#    "Sox1", "Sox2", "Pou3f1", "zfp462", "slc7a3"  # neuroectoderm
#) 

if (is.null(umi) && generate_matrix){
    message("Regenerating matrix from source file ", matrix_destination ,"...")
    counts_raw <- read.table("/extra/sebastian_arnolds/analysis/trim/counts_matrix.tsv",sep="\t")
    barcodes_vector <- as.vector(
        read.table("/extra/sebastian_arnolds/analysis/celseq_barcodes_raw_all.txt")[, 1]
    )
    count_matrix <- as.matrix(counts_raw[, barcodes_vector])

#    # Change labels for genes of interest
#    name_map = list(
#        ENSMUSG00000062327 = "Brachyury",
#        ENSMUSG00000032446 = "Eomes",
#        ENSMUSG00000030544 = "Mesp1",
#        ENSMUSG00000024406 = "Pou5f1",
#        ENSMUSG00000012396 = "nanog",
#        ENSMUSG00000096014 = "Sox1",
#        ENSMUSG00000074637 = "Sox2",
#        ENSMUSG00000090125 = "Pou3f1",
#        ENSMUSG00000060206 = "zfp462",
#        ENSMUSG00000031297 = "slc7a3"
#    )
#
#    for (i in 1:length(name_map)){
#        old_gene_name <- names(name_map[i])
#        new_gene_name <- name_map[[i]]
#        #message(old_gene_name, "→", new_gene_name)
#        rownames(count_matrix)[which(rownames(count_matrix) == old_gene_name)] <- new_gene_name
#    }

    umi <- count_matrix

    # Create SCE and annotate
    sce <- SingleCellExperiment(assays = list(counts = umi))

    # Assign known/related groups
    #is.meso <- rownames(sce) %in% c("Eomes", "Brachyury", "Mesp1")
    #is.pluri <- rownames(sce) %in% c("Pou5f1", "nanog")
    #is.neuro <- rownames(sce) %in% c("Sox1", "Sox2", "Pou3f1", "zfp462", "slc7a3")
    #rowData(sce)$is_mesoderm <- is.meso
    #rowData(sce)$is_pluripotenz <- is.pluri
    #rowData(sce)$is_neuroectoderm <- is.neuro

    # Assign known/related cell batches
    colData(sce)$plate_number <- c(rep("plate1", 96), rep("plate2", 96)) 

    #head(counts(sce))

    sce <- getBMFeatureAnnos(
        sce, 
        filters = "ensembl_gene_id",
        attributes = c(
            "ensembl_gene_id",              # Gene stable ID
            "external_gene_name",           # Casual name
            "external_transcript_name",     # Transcript-specific name
            "gene_biotype",                 # Gene biotype
            "transcript_biotype",           # Trans type
            "description",                  # Gene description
            "band",                         # Karyotype band
            "refseq_mrna",
            "go_id",                        # Go Term accession (cellular domains)
            "go_linkage_type",              # Go Term evidence code
            "name_1006",                    # Go Term name
            "definition_1006",              # Go Term definition
            "namespace_1003"                # Go domain             
        ),
        feature_symbol = "mgi_symbol",
        feature_id = "ensembl_gene_id",
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "mmusculus_gene_ensembl", 
        host = "www.ensembl.org"
    )

    saveRDS(sce, matrix_destination)

} else {
    message("Loading matrix from source file ", matrix_destination, "...")
    sce <- readRDS(matrix_destination)
}


message(dim(sce)[1], " genes x ", dim(sce)[2], " cells. (", 
        length(unique(colnames(sce))), ") unique barcodes.")
