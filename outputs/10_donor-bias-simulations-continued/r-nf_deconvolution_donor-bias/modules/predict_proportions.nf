#!/usr/bin/env nextflow

// defines process to perform downsampling

// set dsl version
nextflow.enable.dsl=2

process predict_proportions {
    publishDir("$params.results_folder", mode: "copy", overwrite: false)

    input:
        val sce_filepath
        val bulk_filepath
        val index_matrix_filepath
        val iterations_index
        val deconvolution_method
        val assay_name
        val celltype_variable
    output:
        path("deconvolution-results_*")

    script:
    """
    Rscript $params.predict_proportions_script -r $sce_filepath -d $deconvolution_method -a $assay_name -c $celltype_variable
    """
}