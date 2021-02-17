#!/usr/bin/env nextflow


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process prepare_diffcomposition_input {
    // Prepare input for R and make basic plots
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(anndata)
        val(cell_label_column)

    output:
        tuple(
            val(outdir),
            path("diff_comp-counts.tsv.gz"),
            path("diff_comp-metadata.tsv.gz"),
            emit: results
        )

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/differential_composition"
        //outdir = "${outdir_prev}/differential_composition/${condition_column}/"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        // --metadata_columns ${condition_column} \
        """
        echo "prepare_diffcomposition_input: ${process_info}"
        echo "publish_directory: ${outdir}"
        013-prepare_diffcomposition_input.py \
            --h5ad_file ${anndata} \
            --experiment_id_column 'experiment_id' \
            --cell_label_column ${cell_label_column} \
            --cell_label_analyse 'all' \
            --metadata_columns 'all' \
            --output_file diff_comp
        """
}


// process run_univariate_test {
//     // Run run_univariate_test in R
//     // ------------------------------------------------------------------------
//     scratch false        // use tmp directory
//     echo echo_mode       // echo output from script
//
//     publishDir  path: "${outdir}",
//                 saveAs: {filename -> filename.replaceAll("${runid}-", "")},
//                 mode: "${task.publish_mode}",
//                 overwrite: "true"
//
//     input:
//         tuple(
//             val(outdir_prev),
//             val(runid),
//             val(condition_column),
//             val(levels),
//             path(df_counts),
//             path(df_meta)
//         )
//
//     // NOTE: mark output as optional because no test if not exactly 2
//     // conditions
//     output:
//         val(outdir, emit: outdir)
//         tuple(
//             val(runid), // need random hex to control grouping
//             path("diff_comp-mannwhitneyu_results.tsv.gz"),
//             emit: results
//         ) optional true
//         path("plots/*.png") optional true
//         path("plots/*.pdf") optional true
//         // tuple(
//         //     val(outdir),
//         //     val(runid), //need random hex to control grouping
//         //     val(condition_column),
//         //     path("diff_comp-univariate_results.tsv.gz"),
//         //     emit: results
//         // ) optional true
//
//     script:
//         runid = random_hex(16)
//         outdir = "${outdir_prev}"
//         cmd__levels = ""
//         if (levels != "") {  // add disc cov call
//             cmd__levels = "--levels ${levels}"
//         }
//         process_info = "${runid} (runid)"
//         process_info = "${process_info}, ${task.cpus} (cpus)"
//         process_info = "${process_info}, ${task.memory} (memory)"
//         """
//         echo "run_univariate_test: ${process_info}"
//         echo "publish_directory: ${outdir}"
//         rm -fr plots
//         014-run_mannwhitneyu_test.R \
//             --counts_file ${df_counts} \
//             --metadata_file ${df_meta} \
//             ${cmd__levels} \
//             --condition "${condition_column}" \
//             --out_file diff_comp
//         if [[ -f "diff_comp-mannwhitneyu_results.tsv.gz" ]]; then
//             015-make_cell_metric_plots.R \
//                 --counts_file diff_comp-counts.tsv.gz \
//                 --metadata_file diff_comp-metadata.tsv.gz \
//                 --condition ${condition_column} \
//                 ${cmd__levels} \
//                 --results_file diff_comp-mannwhitneyu_results.tsv.gz \
//                 --model_plots_only \
//                 --out_file diff_comp-mannwhitneyu_results
//         fi
//         mkdir plots
//         mv *pdf plots/ 2>/dev/null || true
//         mv *png plots/ 2>/dev/null || true
//         """
// }


process run_dirichlet_regression {
    // Run dirichlet_regression in R
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        tuple(
            val(outdir_prev),
            path(df_counts),
            path(df_meta)
        )
        each model

    output:
        val(outdir, emit: outdir)
        path("${out_file}-dirichletreg_model.Rds.gz", emit: model)
        tuple(
            val(runid), // need random hex to control grouping
            path("${out_file}-dirichletreg_results.tsv.gz"),
            emit: results
        )
        path(
            "${out_file}-dirichletreg_forwardmodel.Rds.gz",
            emit: model_forward
        ) optional true
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true
        // tuple(
        //     val(outdir),
        //     val(runid), //need random hex to control grouping
        //     val(condition_column),
        //     path("diff_comp-dirichletreg_model.Rds.gz"),
        //     path("diff_comp-dirichletreg_results.tsv.gz")
        //     emit: results
        // )

    script:
        runid = random_hex(16)
        // Set the outdir using the manipulated formula
        formula_clean = "${model.formula}".replaceAll("_", "")
        formula_clean = "${formula_clean}".replaceAll(" ", "_")
        formula_clean = "${formula_clean}".replaceAll("~", "")
        formula_clean = "${formula_clean}".replaceAll("\\+", "_plus_")
        formula_clean = "${formula_clean}".replaceAll("_", "")
        formula_clean = "${formula_clean}".replaceAll("\\)", "")
        formula_clean = "${formula_clean}".replaceAll("I\\(", "")
        formula_clean = "${formula_clean}".replaceAll("\\^", "power")
        formula_clean = "${formula_clean}".replaceAll("\\.", "pt")
        variable_target_clean = "${model.variable_target}".replaceAll(
            "\\)", ""
        )
        variable_target_clean = "${variable_target_clean}".replaceAll(
            "I\\(", ""
        )
        variable_target_clean = "${variable_target_clean}".replaceAll(
            "\\^", "power"
        )
        variable_target_clean = "${variable_target_clean}".replaceAll(
            "\\.", "pt"
        )
        outdir = "${outdir_prev}/${variable_target_clean}"
        outdir = "${outdir}/formula__${formula_clean}"
        out_file = "diff_comp"
        // Sort out any variables that need to be cast
        cmd__varcast = ""
        if (model.variable_discrete != "") {  // add disc cov call
            cmd__varcast = "${cmd__varcast} --discrete_variables ${model.variable_discrete}"
        }
        if (model.variable_continuous != "") {  // add contin cov call
            cmd__varcast = "${cmd__varcast} --continuous_variables ${model.variable_continuous}"
        }
        // Make discrete levels command
        cmd__levels = ""
        if (model.variable_discrete_level != "") {  // add disc cov call
            cmd__levels = "--discrete_levels \"${model.variable_discrete_level}\""
        }
        // Sort out if this is a forward model
        cmd__forward = ""
        if (model.forward_selection_variables != "") {  // add disc cov call
            cmd__forward = "--forward_selection_variables \"${model.forward_selection_variables}\""
            forward_clean = "${model.forward_selection_variables}".replaceAll(
                "_", ""
            )
            forward_clean = "${forward_clean}".replaceAll("\\,", "_")
            forward_clean = "${forward_clean}".replaceAll("\\)", "")
            forward_clean = "${forward_clean}".replaceAll("I\\(", "")
            forward_clean = "${forward_clean}".replaceAll("\\^", "power")
            forward_clean = "${forward_clean}".replaceAll("\\.", "pt")
            outdir = "${outdir}/forward__${forward_clean}"
        }
        // Details on process
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "run_dirichlet_regression: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        014-run_dirichlet_regression.R \
            --counts_file ${df_counts} \
            --metadata_file ${df_meta} \
            --formula_rhs "${model.formula}" \
            --variable_target "${model.variable_target}" \
            ${cmd__varcast} \
            ${cmd__levels} \
            ${cmd__forward} \
            --optimization_var_threshold 0.05 \
            --out_file ${out_file}
        # No need to include formula info in the out file.
        015-make_cell_metric_plots.R \
            --counts_file ${df_counts} \
            --metadata_file ${df_meta} \
            --results_file ${out_file}-dirichletreg_results.tsv.gz \
            ${cmd__varcast} \
            ${cmd__levels} \
            --out_file ${out_file}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process serialize_result_files {
    // Serializes known markers for analysis
    // ------------------------------------------------------------------------
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    input:
        tuple(
            val(runid),
            file(input_file)
        )

    output:
        path("${runid}-${input_file}", emit: file)

    script:
        //runid = random_hex(16)
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "serialize_result_files: ${process_info}"
        ln --physical ${input_file} ${runid}-${input_file}
        """
}


process merge_and_plot_all_results {
    // Plots all results
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        file(result_files)

    output:
        val(outdir, emit: outdir)
        path("diff_comp-all_results.tsv.gz", emit: results)
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        result_files = result_files.join(",")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "merge_and_plot_all_results: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        # NOTE: for merge_plot_results there may be duplicate formulas due to
        # stepwise forward model selection in run_dirichlet_regression
        016-merge_plot_results.R \
            --results_files ${result_files} \
            --out_file diff_comp
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


workflow wf__differential_composition {
    take:
        outdir
        anndata
        anndata_cell_label
        model
    main:
        // Prep the data to be read into R
        prepare_diffcomposition_input(
            outdir,
            anndata,
            anndata_cell_label
        )
        // Run basic differential composition analysis
        // run_univariate_test(
        //     prepare_diffcomposition_input.out.results
        // )
        // Fit dirichlet_regression in R
        run_dirichlet_regression(
            prepare_diffcomposition_input.out.results,
            model
        )
        // Merge all of the results
        // Here we use groupTuple to keep as channel (rather than collecting
        // to a list). Then we serialize those files and re-name them based
        // in their runid so that we don't have file name clashes
        results_channel = run_dirichlet_regression.out.results.groupTuple(by: 0)
        // results_channel = results_channel.concat(
        //     run_univariate_test.out.results.groupTuple(by: 0)
        // )
        serialize_result_files(results_channel)
        merge_and_plot_all_results(
            "${outdir}/differential_composition",
            serialize_result_files.out.file.collect()
        )


    emit:
        cell_labels = run_dirichlet_regression.out.results
}
