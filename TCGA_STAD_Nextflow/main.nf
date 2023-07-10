#!/usr/bin/env nextflow

process DownloadData {
    
    tag 'r'

    output: 
       file("${params.subdirectory_name}/*.maf.gz")
    
    script:
    """
    Rscript --vanilla ${baseDir}/${params.download_maf_script} ${params.subdirectory_name}
    """
}

process RunMaftools {
    
    stageInMode 'copy'

    tag 'r'

    input:
        file("${params.subdirectory_name}/*.maf.gz")

    output: 
        file("${params.subdirectory_name}/binary_numerical_merged.csv")
    
    script:
    """
    Rscript --vanilla ${baseDir}/${params.maftools_script} ${params.subdirectory_name}
    """
}

process RunPythonScript {
    
    tag 'python'
    publishDir 'results', mode: 'copy'

    input:
        tuple val(hyperparam), val(distance_type), val(cluster)
        file("${params.subdirectory_name}/binary_numerical_merged.csv")

    output: 
        file('results/*.csv')
    
    script:
    """
    python3 ${baseDir}/${params.py_script} ${hyperparam} ${distance_type} ${cluster} "${params.subdirectory_name}/binary_numerical_merged.csv"
    """
}

// Create a channel with the combinations of parameters
combinations = Channel.from(params.hyperparams)
                      .combine(Channel.from(params.distance_types))
                      .combine(Channel.from(params.clusters))

workflow {
    files = DownloadData()
    tabular_data = RunMaftools(files)
    RunPythonScript(combinations, tabular_data)
}
