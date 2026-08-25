// ── CGMLST_DISTANCES ──────────────────────────────────────────────────────────
//
// Step 10b: Pairwise cgMLST allele distance matrix from ChewBBACA output.
//
// Calls bin/cgmlst_distances.py (Tier 3, unit-tested separately).
// Runs after CHEWBBACA; receives cgmlst_allele_matrix.tsv.
//
// Distance: count of loci where BOTH samples have valid allele calls and the
// calls differ. Loci where either sample has a non-called ChewBBACA designator
// (LNF, PLOT3, PLOT5, LOTSC, NIPH, NIPHEM, ALM, ASM, PAMA, -) are excluded.
// INF-N (inferred new alleles) are treated as valid calls.
// When ALL loci are non-called for a pair, distance is emitted as NA.
//
// If the input matrix is the not_run placeholder (emitted by CHEWBBACA when no
// schema is available), the script re-emits a not_run placeholder and exits 0.
// ─────────────────────────────────────────────────────────────────────────────

process CGMLST_DISTANCES {
    tag "cgmlst_distances"
    label 'process_single'

    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path allele_matrix   // cgmlst_allele_matrix.tsv from CHEWBBACA

    output:
    path 'cgmlst_distance_matrix.tsv', emit: distance_matrix

    script:
    """
    export PATH="${moduleDir.parent.parent}/bin:\${PATH}"

    cgmlst_distances.py \\
        --allele-matrix ${allele_matrix} \\
        --output        cgmlst_distance_matrix.tsv
    """
}
