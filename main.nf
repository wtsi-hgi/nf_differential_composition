#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

VERSION = "0.0.1" // Do not edit, controlled by bumpversion.


// Modules to include.
include {
    wf__differential_composition;
} from "./modules/differential_composition.nf"

// Set default parameters.
params.output_dir           = "nf-differential_composition"
params.help                 = false
params.anndata_cell_label = [value: 'cluster']
// Default parameters for differential composition
params.differential_composition = [
    run_process: false
]

// Define the help messsage.
def help_message() {
    log.info """
    ============================================================================
     single cell differential condition ~ v${VERSION}
    ============================================================================

    Runs basic single cell preprocessing

    Usage:
    nextflow run main.nf -profile <local|lsf> -params-file params.yaml [options]

    Mandatory arguments:
        --file_anndata      Anndata file with cell type labels.

    Optional arguments:
        --output_dir        Directory name to save results to. (Defaults to
                            '${params.output_dir}')

        -params-file        YAML file containing analysis parameters. See
                            example in example_runtime_setup/params.yml.

    Profiles:
        lsf                 lsf cluster execution
    """.stripIndent()
}


// Boot message - either help message or the parameters supplied.
if (params.help){
    help_message()
    exit 0
} else {
    log.info """
    ============================================================================
     single cell differential composition ~ v${VERSION}
    ============================================================================
    file_anndata                  : ${params.file_anndata}
    output_dir (output folder)    : ${params.output_dir}
    """.stripIndent()
}


// Initalize Channels.
// anndata = Channel
//     .fromPath(params.file_anndata)


// Run the workflow.
workflow {
    main:
        // Run differential cell type composition analysis
        if (params.differential_composition.run_process) {
            wf__differential_composition(
                params.output_dir,
                params.file_anndata,
                params.anndata_cell_label.value,
                params.differential_composition.models.value
            )
        }
    
    // NOTE: One could do publishing in the workflow like so, however
    //       that will not allow one to build the directory structure
    //       depending on the input data call. Therefore, we use publishDir
    //       within a process.
    // publish:
    //     merge_samples.out.anndata to: "${params.output_dir}",
    //         mode: "copy",
    //         overwrite: "true"
}


workflow.onComplete {
    // executed after workflow finishes
    // ------------------------------------------------------------------------
    log.info """\n
    ----------------------------------------------------------------------------
     pipeline execution summary
    ----------------------------------------------------------------------------
    Completed         : ${workflow.complete}
    Duration          : ${workflow.duration}
    Success           : ${workflow.success}
    Work directory    : ${workflow.workDir}
    Exit status       : ${workflow.exitStatus}
    """.stripIndent()
}
