#!/usr/bin/env nextflow

// This is the main workflow script for example 4.

// set dsl version
nextflow.enable.dsl=2

// include modules for workflow

include { predict_proportions as predict_prop } from "$launchDir/modules/predict_proportions"
include { analyze_results as analyze_res } from "$launchDir/modules/analyze_results"

// define a new workflow
workflow {
    
    // parse channels
    
    sce_filepath = channel.fromList( params.sce_filepath ) // single-cell expression data
    
    bulk_filepath = channel.fromList( params.bulk_filepath ) // pseudo-bulk expression data
    
    mi_filepath = channel.fromList( params.index_matrix_filepath ) // indices matrix filepath
    
    iterations_index = channel.fromList( params.iterations_index ) // index of the current run
    
    method = channel.fromList( params.method ) // deconvolution method name
    
    assay_name = channel.fromList( params.assay_name ) // name of the assay type to use
    
    celltype_variable = channel.fromList( params.celltype_variable ) // name of celltype variable
    
    tp_filepath = channel.fromList( params.true-proportions_filepath ) // true proportions data
    
    
    // run workflow
    predict_prop(   
                    sce_filepath, 
                    bulk_filepath, 
                    mi_filepath, 
                    iterations_index, 
                    method, 
                    assay_name, 
                    celltype_variable
                )
    predict_prop.out.view()
    analyze_res( predict_prop.out, tp_filepath )
}

workflow.onComplete {
    println "Pipeline successfully completed at: $workflow.complete"
}

workflow.onError {
    println "Pipeline stopped with following error: $workflow.errorMessage"
}

