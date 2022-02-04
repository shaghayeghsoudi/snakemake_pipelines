#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Shaghayegh Soudi
# Module Author:    Shaghayegh Soudi
# Contributors:    NA 

##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op


# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section    



# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["dpclust"]`
CFG = op.setup_module(
    name = "dpclust",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs","dpclust3p","dpclust","purity_adjusted","dpclust3p_purityadj","dpclust_purityadj","outputs"],
)


# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _dpclust_input_maf,
    _dpclust_input_battenberg,
    _dpclust_preprocess_output,
    _dpclust_master_output,
    _dpclust_master_best_purity,
    _dpclust_output,  
    _dpclust_all


##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _dpclust_input_maf:
    input:       
        maf = CFG["inputs"]["maf"]  
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched_slms-3.final.maf"  
    run:
        op.absolute_symlink(input.maf, output.maf)  


rule _dpclust_input_battenberg:
    input: 
        subclones = CFG["inputs"]["subclones"],     
        cellularity = CFG["inputs"]["cellularity"],
        rho_psi = CFG["inputs"]["rho_psi"]
    output:
        subclones = CFG["dirs"]["inputs"] + "battenberg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}.subclones.txt",
        cellularity = CFG["dirs"]["inputs"] + "battenberg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}.cellularity_ploidy.txt",
        rho_psi = CFG["dirs"]["inputs"] + "battenberg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_rho_and_psi.txt"
    run:
        op.absolute_symlink(input.subclones, output.subclones)
        op.absolute_symlink(input.cellularity, output.cellularity)
        op.absolute_symlink(input.rho_psi, output.rho_psi)


rule _dpclust_input_sex:
    input: 
        sex = CFG["inputs"]["sex"]
    output:
        sex = CFG["dirs"]["inputs"] + "battenberg/{seq_type}--{genome_build}/{normal_id}.sex"       
    run:
        op.absolute_symlink(input.sex, output.sex)
        

rule _install_dpclust:
    output:
        complete = CFG["dirs"]["inputs"] + "dpclust_dependencies_installed.success"
    conda:
        CFG["conda_envs"]["dpclust"]
    log:
        input = CFG["logs"]["inputs"] + "input.log"
    shell:
        """
        R -q -e 'BiocManager::install(c("GenomicRanges","IRanges"))' &&
        R -q -e 'devtools::install_github("shaghayeghsoudi/dpclust3p")' >> {log.input} && 
        R -q -e 'devtools::install_github("Wedge-Oxford/dpclust")' >> {log.input} &&
        touch {output.complete}"""


# make dpclust3p (makes dpclust preprocessing inputs, script modified to take new vcf format)
rule _dpclust_preprocess_output:
    input:
        maf = str(rules._dpclust_input_maf.output.maf),
        subclones = str(rules._dpclust_input_battenberg.output.subclones),
        rho_psi = str(rules._dpclust_input_battenberg.output.rho_psi),
        sex = str(rules._dpclust_input_sex.output.sex),
        installed = CFG["dirs"]["inputs"] + "dpclust_dependencies_installed.success"                          
    output:       
        dpclust_loci = CFG["dirs"]["dpclust3p"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_loci.txt",
        dpclust_allelefrequencies = CFG["dirs"]["dpclust3p"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_alleleFrequencies.txt",
        dpclust_alldirichletprocessinfo = CFG["dirs"]["dpclust3p"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_allDirichletProcessInfo.txt"   
    log:
        stderr = CFG["logs"]["dpclust3p"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}.create_preprocess_dpclus_inputs.stderr.log",
        stdout = CFG["logs"]["dpclust3p"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}.create_prepeocess_dpclus_inputs.stdout.log"
    params:                             
        out_dir_pre = CFG["dirs"]["dpclust3p"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}",
        script = CFG["inputs"]["dpclust3p_script"]    
    conda:
        CFG["conda_envs"]["dpclust"]
    shell:
        op.as_one_line("""
        sex=$(tail -n +2 {input.sex} | cut -f 4);      
        Rscript {params.script} -s {wildcards.tumour_id} -m {input.maf} -r {input.rho_psi} -c {input.subclones} -o {params.out_dir_pre} --sex $sex 
        2> {log.stderr} > {log.stdout} 
        """)
        

rule _dpclust_master_output:
    input:
        cellularity = str(rules._dpclust_input_battenberg.output.cellularity),
        rho_psi = str(rules._dpclust_input_battenberg.output.rho_psi),
        sex = str(rules._dpclust_input_sex.output.sex),
        installed = CFG["dirs"]["inputs"] + "dpclust_dependencies_installed.success",
        dpclust_alldirichletprocessinfo = CFG["dirs"]["dpclust3p"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_allDirichletProcessInfo.txt"                           
    output:
        master = CFG["dirs"]["dpclust3p"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}_master.txt"
    log:    
        stderr = CFG["logs"]["dpclust3p"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/create_master_input.stderr.log"
    params:
        script_master = CFG["inputs"]["dpclust3p_master"]
    conda:
        CFG["conda_envs"]["dpclust"]    
    shell:
        op.as_one_line("""
        sex=$(tail -n +2 {input.sex} | cut -f 4); 
        Rscript {params.script_master} -s {wildcards.tumour_id} -d {wildcards.tumour_id} -p {input.cellularity} -r {input.rho_psi} --output {output.master} --sex $sex
        2> {log.stderr}
        """)


# This rule runs the entire dpclust pipeline.
dpclust_extension1 = "_DPoutput_" + str(CFG["options"]["dpclust_run"]["iters"]) + "iters_" + str(CFG["options"]["dpclust_run"]["burnin"]) + "burnin_seed" + str(CFG["options"]["dpclust_run"]["seed"])
dpclust_extension2 =  str(CFG["options"]["dpclust_run"]["iters"]) + "iters_" + str(CFG["options"]["dpclust_run"]["burnin"]) + "burnin"
rule _dpclust_run: 
    input:
        master = CFG["dirs"]["dpclust3p"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}_master.txt",
        dpclust_alldirichletprocessinfo = CFG["dirs"]["dpclust3p"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_allDirichletProcessInfo.txt" 
    output:    
        bestcluster = CFG["dirs"]["dpclust"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}" + dpclust_extension1 + "/{tumour_id}_" + dpclust_extension2 + "_bestClusterInfo.txt",
        bestconsensus_bed = CFG["dirs"]["dpclust"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}" + dpclust_extension1 + "/{tumour_id}_" + dpclust_extension2 + "_bestConsensusAssignments.bed",
        bestconsensus_results = CFG["dirs"]["dpclust"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}" + dpclust_extension1 + "/{tumour_id}_" + dpclust_extension2 + "_bestConsensusResults.RData",
        mutation_assignments = CFG["dirs"]["dpclust"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}" + dpclust_extension1 + "/{tumour_id}_" + dpclust_extension2 + "_mutation_assignments.png",
        mutclus_likelihoods = CFG["dirs"]["dpclust"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}" + dpclust_extension1 + "/{tumour_id}_" + dpclust_extension2 + "_mutationClusterLikelihoods.bed",
        cluster_locations = CFG["dirs"]["dpclust"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}" + dpclust_extension1 + "/{tumour_id}_DirichletProcessplot_with_cluster_locations_2.png"    
    log:
        stderr = CFG["logs"]["dpclust"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}" + dpclust_extension2 + "/dpclust.stderr.log", 
        stdout = CFG["logs"]["dpclust"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}" + dpclust_extension2 + "/dpclust.stdout.log"   
    params:      
        script_main= CFG["inputs"]["dpclust_script"],
        input_dir = CFG["dirs"]["dpclust3p"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}",
        output_dir = CFG["dirs"]["dpclust"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}"
    conda: 
         CFG["conda_envs"]["dpclust"]
    shell:
        op.as_one_line("""
        Rscript {params.script_main} -r 1 -d {params.input_dir} -i {input.master} -o {params.output_dir}
        2> {log.stderr} > {log.stdout}
        """)

### this rule estimates the best purity from DPclust output
rule _dpclust_adjust_purity: 
    input:
        clonal_location = CFG["dirs"]["dpclust"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}" + dpclust_extension1 + "/{tumour_id}_" + dpclust_extension2 + "_bestClusterInfo.txt",
        rho_psi = str(rules._dpclust_input_battenberg.output.rho_psi)       
    output: 
        estimated_purity = CFG["dirs"]["purity_adjusted"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_dpclust_estimated_purity.txt"
    log:
        stderr = CFG["logs"]["purity_adjusted"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}.dpclust.stderr.log"  
    params:      
        script = CFG["scripts"]["dpclust_best_purity"],
        best_purity_out_dir = CFG["dirs"]["purity_adjusted"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}"
    conda:
        CFG["conda_envs"]["dpclust"]    
    shell:
        op.as_one_line("""
        Rscript {params.script} -n {wildcards.tumour_id} -r {input.rho_psi} -c {input.clonal_location} -o {params.best_purity_out_dir}
        """)


# make dpclust3p (makes dpclust preprocessing inputs, from new estimated purity)
rule _dpclust3p_adjusted_purity:
    input:
        maf = str(rules._dpclust_input_maf.output.maf),
        subclones = str(rules._dpclust_input_battenberg.output.subclones),
        rho_psi_best = str(rules._dpclust_adjust_purity.output.estimated_purity),     
        sex = str(rules._dpclust_input_sex.output.sex)
    output:       
        bp_loci = CFG["dirs"]["dpclust3p_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_loci.txt",
        bp_allelefrequencies = CFG["dirs"]["dpclust3p_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_alleleFrequencies.txt",
        bp_alldirichletprocessinfo = CFG["dirs"]["dpclust3p_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_allDirichletProcessInfo.txt"  
    log:
        stderr = CFG["logs"]["dpclust3p_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_create_preprocess_best_purity_dpclus_inputs.stderr.log",
        stdout = CFG["logs"]["dpclust3p_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_create_prepeocess_best_purity_dpclus.stdout.log"
    params:
        out_dir_pre_bestpurity = CFG["dirs"]["dpclust3p_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}",
        script = CFG["inputs"]["dpclust3p_script"]    
    conda:
        CFG["conda_envs"]["dpclust"]
    shell:
        op.as_one_line("""
        sex=$(tail -n +2 {input.sex} | cut -f 4);
        Rscript {params.script} -s {wildcards.tumour_id} -m {input.maf} -r {input.rho_psi_best} -c {input.subclones} -o {params.out_dir_pre_bestpurity} --sex $sex
        2> {log.stderr} > {log.stdout}
        """)


rule _dpclust_master_best_purity:
    input:
        cellularity = str(rules._dpclust_input_battenberg.output.cellularity),
        rho_psi = str(rules._dpclust_adjust_purity.output.estimated_purity),   
        sex = str(rules._dpclust_input_sex.output.sex)
    output:
        master_bp = CFG["dirs"]["dpclust3p_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}_best_purity_master.txt"
    log:    
        stderr = CFG["logs"]["dpclust3p_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/create_master_input.stderr.log"
    params:
        script_master = CFG["inputs"]["dpclust3p_master"]
    conda:
        CFG["conda_envs"]["dpclust"]    
    shell:
        op.as_one_line("""
        sex=$(tail -n +2 {input.sex} | cut -f 4); 
        Rscript {params.script_master} -s {wildcards.tumour_id} -d {wildcards.tumour_id} -p {input.cellularity} -r {input.rho_psi} --output {output.master_bp} --sex $sex
        2> {log.stderr}
        """)


# This rule runs the entire dpclust pipeline (bp == "best_purity").
dpclust_extension1 = "_DPoutput_" + str(CFG["options"]["dpclust_run"]["iters"]) + "iters_" + str(CFG["options"]["dpclust_run"]["burnin"]) + "burnin_seed" + str(CFG["options"]["dpclust_run"]["seed"])
dpclust_extension2 = str(CFG["options"]["dpclust_run"]["iters"]) + "iters_" + str(CFG["options"]["dpclust_run"]["burnin"]) + "burnin"

rule _dpclust_bp_run: 
    input:
        master = CFG["dirs"]["dpclust3p_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}_best_purity_master.txt",  
        dpclust_alldirichletprocessinfo_bp = CFG["dirs"]["dpclust3p_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_allDirichletProcessInfo.txt"      
    output:      
        bp_bestcluster = CFG["dirs"]["dpclust_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}" + dpclust_extension1 + "/{tumour_id}_" + dpclust_extension2 + "_bestClusterInfo.txt",
        bp_bestconsensus_bed = CFG["dirs"]["dpclust_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}" + dpclust_extension1 + "/{tumour_id}_" + dpclust_extension2 + "_bestConsensusAssignments.bed",
        bp_bestconsensus_results = CFG["dirs"]["dpclust_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}" + dpclust_extension1 + "/{tumour_id}_" + dpclust_extension2 + "_bestConsensusResults.RData",
        bp_mutation_assignments = CFG["dirs"]["dpclust_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}" + dpclust_extension1 + "/{tumour_id}_" + dpclust_extension2 + "_mutation_assignments.png",
        bp_mutclus_likelihoods = CFG["dirs"]["dpclust_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}" + dpclust_extension1 + "/{tumour_id}_" + dpclust_extension2 + "_mutationClusterLikelihoods.bed",
        bp_cluster_locations = CFG["dirs"]["dpclust_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}" + dpclust_extension1 + "/{tumour_id}_DirichletProcessplot_with_cluster_locations_2.png"    
    log:
        stderr = CFG["logs"]["dpclust_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}" + dpclust_extension1 + "/dpclust.purity_adjusted.stderr.log", 
        stdout = CFG["logs"]["dpclust_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}" + dpclust_extension1 + "/dpclust.purity_adjusted.stdout.log"  
    params:      
        script_main= CFG["inputs"]["dpclust_script"],
        input_dir = CFG["dirs"]["dpclust3p_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}",
        out_dir = CFG["dirs"]["dpclust_purityadj"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}"
    conda: 
         CFG["conda_envs"]["dpclust"]   
    shell:
        op.as_one_line("""
        Rscript {params.script_main} -r 1 -d {params.input_dir} -i {input.master} -o {params.out_dir}
        2> {log.stdout} > {log.stderr}
        """)
    

# Symlinks the final output files into the module results directory (und er '99-outputs/')
# All plots generated by Battenberg are symlinked using a glob for convenience
dpclust_extension2 =  str(CFG["options"]["dpclust_run"]["iters"]) + "iters_" + str(CFG["options"]["dpclust_run"]["burnin"]) + "burnin"

rule _dpclust_output:
    input:
        bestcluster = rules._dpclust_bp_run.output.bp_bestcluster,
        bestconsensus_bed = rules._dpclust_bp_run.output.bp_bestconsensus_bed,
        mutation_assignments = rules._dpclust_bp_run.output.bp_mutation_assignments,
        likelihoods = rules._dpclust_bp_run.output.bp_mutclus_likelihoods,
        cluster_locations = rules._dpclust_bp_run.output.bp_cluster_locations
    output:
        bestcluster = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}_" + dpclust_extension2 + "_bestClusterInfo.txt",
        bestconsensus_bed = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}_" + dpclust_extension2 + "_bestConsensusAssignments.bed", 
        mutation_assignments = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}_" + dpclust_extension2 + "_mutation_assignments.png",
        mutclus_likelihoods = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}_" + dpclust_extension2 + "_mutationClusterLikelihoods.bed",
        cluster_locations = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}_DirichletProcessplot_with_cluster_locations_2.png"  
    run:
        op.relative_symlink(input.bestcluster, output.bestcluster,in_module=True)
        op.relative_symlink(input.bestconsensus_bed, output.bestconsensus_bed,in_module=True)
        op.relative_symlink(input.mutation_assignments, output.mutation_assignments,in_module=True)
        op.relative_symlink(input.likelihoods, output.mutclus_likelihoods,in_module=True)
        op.relative_symlink(input.cluster_locations, output.cluster_locations,in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _dpclust_all:
    input:
        expand(
            [
                rules._dpclust_output.output.bestcluster,
                rules._dpclust_output.output.bestconsensus_bed,
                rules._dpclust_output.output.mutation_assignments,
                rules._dpclust_output.output.mutclus_likelihoods,
                rules._dpclust_output.output.cluster_locations
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])



##### CLEANUP #####

# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable

