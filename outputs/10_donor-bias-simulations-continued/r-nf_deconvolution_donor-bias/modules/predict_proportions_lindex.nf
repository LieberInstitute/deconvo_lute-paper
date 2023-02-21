#!/usr/bin/env nextflow

// defines process to perform downsampling

// set dsl version
nextflow.enable.dsl=2

process predict_proportions_lindex {
    publishDir("$params.results_folder", mode: "copy", overwrite: false)

    maxForks 40

    input:
        val sce_filepath
        val bulk_filepath
        val list_index_fpath
        val iterations_index
        val deconvolution_method
        val assay_name
        val celltype_variable
    output:
        path("deconvolution-results_*")

    script:
    """
    Rscript $params.predict_proportions_script -r $sce_filepath -b $bulk_filepath -l $list_index_fpath -i $iterations_index -d $deconvolution_method -a $assay_name -c $celltype_variable
    """
}