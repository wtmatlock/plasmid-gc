#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// -------------------------
// PROCESS: COMPUTE CONTIG GC AND LENGTH
// -------------------------
process COMPUTE_CONTIG_STATISTICS {
    container "quay.io/biocontainers/biopython:1.79"
    cpus params.threads
    memory params.memory
    
    publishDir "${params.outdir}", mode: "copy"

    output:
    path "contig_statistics.tsv", emit: contig_statistics

    script:
    """
    ls -d "${params.project_dir}/${params.contigs_dir}"/*.fasta > fasta_paths.txt
    python "${params.project_dir}/scripts/contig_statistics.py" fasta_paths.txt contig_statistics.tsv
    """
}

// -------------------------
// METADATA CHANNEL
// -------------------------
meta_ch = Channel.fromPath(params.metadata)
    .splitCsv(header: true, sep: ',')
    .map { row ->
        row['contig_id'] = row.get('contig_id') ?: row.get('contig') ?: ''
        row['type']      = (row.get('type') ?: '').toLowerCase()
        row['fasta']     = file("${params.contigs_dir}/${row.contig_id}.fasta")
        return row
    }

// -------------------------
// COLLECT PLASMID FASTAS
// -------------------------
plasmid_fastas_ch = meta_ch
    .filter { it.type == 'plasmid' }
    .map { it.fasta.toString() }
    .collect()
    .map { it.sort() }
    
// -------------------------
// PROCESS: RUN PLASMIDFINDER VIA ABRICATE
// -------------------------
process RUN_PLASMIDFINDER {
    container "quay.io/biocontainers/abricate:1.2.0--h05cac1d_0"
    cpus params.threads
    memory params.memory

    publishDir "${params.outdir}", mode: "copy"

    input:
    val fasta_list   // list of plasmid FASTA paths

    output:
    path "plasmidfinder.tsv"

    script:
    """
    # Concatenate all FASTAs into one file
    for f in ${fasta_list.join(' ')}; do
        cat "\$f"
        echo
    done > all_plasmids.fasta
    
    # Run abricate
    abricate --threads $task.cpus --db plasmidfinder all_plasmids.fasta > plasmidfinder.tsv
    """
}

// -------------------------
// PROCESS: RUN MOB-TYPER
// -------------------------
process RUN_MOBTYPER {
    container "quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1"
    cpus params.threads
    memory params.memory

    publishDir "${params.outdir}", mode: "copy"

    input:
    val fasta_list   // list of plasmid FASTA paths

    output:
    path "mobtyper.tsv"

    script:
    """
    # Concatenate all FASTAs into one file
    for f in ${fasta_list.join(' ')}; do
        cat "\$f"
        echo
    done > all_plasmids.fasta
    
    # Run MOB-typer
    mob_typer --num_threads $task.cpus --multi --infile all_plasmids.fasta --out_file mobtyper.tsv
    """
}

// -------------------------
// PROCESS: RUN PLING
// -------------------------
process RUN_PLING {
    container "staphb/pling:2.0.1"
    cpus params.threads
    memory params.memory

    publishDir "${params.outdir}", mode: "copy"

    input:
    val fasta_list

    output:
    path "pling_results/dcj_thresh_4_graph/objects/typing.tsv", emit: plasmid_clusters
    path "pling_results/**"

    script:
    """
    mkdir -p pling_results
    
    ls -d ${fasta_list.join(' ')} > input.txt

    pling input.txt pling_results align \
        --containment_distance 0.3 \
        --cores $task.cpus \
        --batch_size 2000
    """
}

// -------------------------
// COLLECT CHROMOSOME FASTAS
// -------------------------
chromosome_fastas_ch = meta_ch
    .filter { it.type == 'chromosome' }
    .map { it.fasta.toString() }
    .collect()
    .map { it.sort() }
    
// -------------------------
// PROCESS: RUN MASHTREE
// -------------------------
process RUN_MASHTREE {
    container "quay.io/biocontainers/mashtree:1.4.6--pl5321h7b50bb2_3"
    cpus params.threads
    memory params.memory

    publishDir "${params.outdir}", mode: "copy"

    input:
    val fasta_list   // list of chromosome FASTA paths

    output:
    path "chromosome_mashtree.dnd"

    script:
    """  
    # Run mashtree
    mashtree --mindepth 0 --numcpus $task.cpus ${fasta_list.join(' ')} > chromosome_mashtree.dnd
    """
}
    
// -------------------------
// PROCESS: RUN MLST
// -------------------------
process RUN_MLST {
    container "quay.io/biocontainers/mlst:2.28.1--hdfd78af_0"
    cpus params.threads
    memory params.memory

    publishDir "${params.outdir}", mode: "copy"

    input:
    val fasta_list   // list of chromosome FASTA paths

    output:
    path "mlst.tsv", emit: mlst_groups

    script:
    """  
    # Run mlst
    mlst ${fasta_list.join(' ')} > mlst.tsv
    """
}


// -------------------------
// WORKFLOW
// -------------------------
workflow {

    // --- Compute contig statistics
    contig_stats = COMPUTE_CONTIG_STATISTICS()
    
    // --- Run PlasmidFinder on all plasmids
    plasmid_fastas_ch
        | RUN_PLASMIDFINDER
        
    // --- Run MOB-typer on all plasmids
    plasmid_fastas_ch
        | RUN_MOBTYPER

    // --- Plasmid clustering
    plasmid_clustering_ch = plasmid_fastas_ch 
        | RUN_PLING
                
    // --- Run mashtree on all chromosomes
    mashtree_results = chromosome_fastas_ch
        | RUN_MASHTREE
        
    // --- Run MLST on all chromosomes
    mlst_results = chromosome_fastas_ch
        | RUN_MLST

    emit:
    chromosome_gene_igr_gc_ch
}
