#!/usr/bin/env nextflow

params.samplesheet = null
params.fraction = 0.5

process SUBSET_CELLS {
    tag "${sample_name}"
    cpus 1
    memory 8.GB

    container 'satijalab/seurat:latest'

    input:
    tuple val(sample_name), path(matrix_dir)
    val fraction

    output:
    tuple val(sample_name), path("${sample_name}_subset"), emit: subset
    path("${sample_name}_top.txt"), emit: top_genes

    script:
    """
    #!/usr/bin/env Rscript
    library(Seurat)
    library(Matrix)

    counts <- Read10X(data.dir = "${matrix_dir}")
    obj <- CreateSeuratObject(counts = counts)

    totals <- Matrix::colSums(GetAssayData(obj, layer = "counts"))
    n_keep <- max(1, floor(length(totals) * as.numeric("${fraction}")))
    keep <- names(sort(totals, decreasing = TRUE))[seq_len(n_keep)]
    sub <- subset(obj, cells = keep)

    mat <- GetAssayData(sub, layer = "counts")

    outdir <- "${sample_name}_subset"
    dir.create(outdir, showWarnings = FALSE)
    Matrix::writeMM(mat, file.path(outdir, "matrix.mtx"))
    write.table(colnames(mat), file.path(outdir, "barcodes.tsv"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(data.frame(id = rownames(mat), name = rownames(mat), type = "Gene Expression"),
                file.path(outdir, "features.tsv"),
                sep = "\\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    gene_totals <- Matrix::rowSums(mat)
    top3 <- names(sort(gene_totals, decreasing = TRUE))[seq_len(min(3, length(gene_totals)))]
    writeLines(top3, "${sample_name}_top.txt")
    """

    stub:
    """
    mkdir -p ${sample_name}_subset
    touch ${sample_name}_subset/matrix.mtx
    touch ${sample_name}_subset/barcodes.tsv
    touch ${sample_name}_subset/features.tsv
    printf "GENE1\\nGENE2\\nGENE3\\n" > ${sample_name}_top.txt
    """
}

workflow {
    main:
    samples = channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], file(row[1], checkIfExists: true)) }

    SUBSET_CELLS(samples, params.fraction)

    SUBSET_CELLS.out.top_genes
        .splitText()
        .map { it -> it.trim() }
        .filter { it -> it }
        .unique()
        .collect()
        .view { genes -> "Union of top 3 expressed genes: ${genes.sort().join(', ')}" }

    publish:
    subset_dirs = SUBSET_CELLS.out.subset.map { _sample, dir -> dir }
    top_genes = SUBSET_CELLS.out.top_genes
}

output {
    subset_dirs {
        mode 'copy'
    }
    top_genes {
        mode 'copy'
    }
}
