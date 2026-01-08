nextflow.enable.dsl=2

// ---------------------- PARAMS ----------------------
params.exp                 = ''
params.ptid                = ''
params.sample_id           = ''
params.workflow            = ''
params.species             = ''
params.outdir              = ''
params.cellranger          = ''

// References
params.human_transcriptome = '/datacommons/dhvi/10X_Genomics_data/references/GEX/refdata-gex-GRCh38-2020-A'
params.rhesus_transcriptome= '/datacommons/dhvi/10X_Genomics_data/references/GEX/MMul10_GEX_release115'
params.mouse_transcriptome = '/datacommons/dhvi/10X_Genomics_data/references/GEX/refdata-gex-mm10-2020-A'
params.human_vdj_ref       = '/datacommons/dhvi/10X_Genomics_data/references/VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0'
params.rhesus_vdj_ref      = '/datacommons/dhvi/10X_Genomics_data/references/VDJ/Rhesus_Compiled'
params.mouse_vdj_ref       = '/datacommons/dhvi/10X_Genomics_data/references/VDJ/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0'
params.inner_enrichment    = "/datacommons/dhvi/mb488/scripts/10x/cellranger/inner_enrichment_primers.txt"
params.feature_ref         = ''

// Analysis toggles
params.control             = ''
params.clonal              = 'c'   // c=cloanalyst, p=partis, b=both
params.germline            = ''

// mkfastq inputs
params.gex_raw             = ''
params.vdj_raw             = ''
params.beam_raw            = ''
params.gex_csv             = ''
params.vdj_csv             = ''
params.beam_csv            = ''

// Mapping-file inputs (CSV/TSV; allows multiple rows per sample)
params.gex_fastq_map       = ''    // columns: sample,fastq_dir
params.vdj_fastq_map       = ''
params.beam_fastq_map      = ''


// ---------------------- HELPERS ----------------------
/** Write a sample_ids.txt (quoted, one per line) */
process write_sample_ids {
  input:
    val sample_ids
  output:
    path "sample_ids.txt", emit: sample_file
  script:
  """
  printf "%s\\n" ${ sample_ids.collect{"\\\"${it}\\\""}.join(' ') } | xargs -n1 > sample_ids.txt
  """
}

/** Extract sample IDs from a mkfastq CSV (2nd column, skip header) */
process extract_sample_ids_from_csv {
  input:
    path csv
  output:
    path "sample_ids.txt", emit: sample_file
  shell:
  '''
  awk -F, 'NR>1 && $2!="" {print $2}' !{csv} > sample_ids.txt
  '''
}

/** Parse a CSV/TSV map file into:
 *  - tuplesCh:  (sample, "dir1,dir2,...")
 *  - namesCh:   sample
 * We avoid channel-level groupBy by collecting then grouping in Groovy.
 */
def parse_map_to_tuples(mapPath) {
  // choose separator by extension (override if you like)
  def p   = mapPath.toString().toLowerCase()
  def sep = p.endsWith('.tsv') ? '\t' : ','

  // tuples: (sample, "dir1,dir2,...")
  def tuples = nextflow.Channel
      .fromPath(mapPath)
      .splitCsv(header:true, sep: sep)
      .filter { rec -> rec.sample && rec.fastq_dir }               // keep valid rows
      .map    { rec -> tuple(rec.sample.trim(), rec.fastq_dir.trim()) }
      .groupTuple()                                                // group by sample
      .map    { sample, dirs -> tuple(sample, dirs.join(',')) }

  // names: just the sample names (distinct by construction)
  def names = tuples.map { it[0] }

  return [tuples, names]
}


// ---------------------- GEX ----------------------
process makefastq_gex {
  tag "${sample_id}"
  publishDir "${params.outdir}/GEX/cellranger", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    path gex_raw
    val sample_id
    path gex_csv
    val cellranger

  output:
    path "${sample_id}_mkfastq_outs", emit: fastqs
    path "sample_ids.txt",           emit: sample_file

  shell:
  '''
  /datacommons/dhvi/10X_Genomics_data/cellranger/cellranger-!{cellranger}/cellranger mkfastq \
    --run=!{gex_raw} --csv=!{gex_csv} --output-dir=!{sample_id}_mkfastq_outs --delete-undetermined
  awk -F, 'NR > 1 {print $2}' !{gex_csv} > sample_ids.txt
  '''
}

process count_gex {
  tag "${sample}"
  publishDir "${params.outdir}/GEX/cellranger", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    tuple val(sample), val(fastqs)
    val species
    val cellranger

  output:
    path "${sample}", emit: current_sample

  script:
    def count_script = cellranger.startsWith("9") ? "count_new.sh" : "count.sh"
    if (species == "human")
      """
      echo "Using human transcriptome"
      ${count_script} ${fastqs} ${sample} ${params.human_transcriptome} ${cellranger}
      """
    else if (species == "rhesus")
      """
      echo "Using rhesus transcriptome"
      ${count_script} ${fastqs} ${sample} ${params.rhesus_transcriptome} ${cellranger}
      """
    else if (species == "mouse")
      """
      echo "Using mouse transcriptome"
      ${count_script} ${fastqs} ${sample} ${params.mouse_transcriptome} ${cellranger}
      """
}

process agg_gex {
  publishDir "${params.outdir}/GEX/cellranger", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    val counts_paths
    path sample_file
    val sample_id
    val cellranger

  output:
    path "${sample_id}_Agg"
    path "${sample_id}_Agg.csv"

  script:
  """
  echo ${counts_paths} > counts_paths.txt
  sed -i 's/[][]//g; s/,\\s*/,/g' counts_paths.txt
  make_agg_csv2.sh counts_paths.txt ${sample_id}_Agg.csv
  /datacommons/dhvi/10X_Genomics_data/cellranger/cellranger-${cellranger}/cellranger aggr \
     --id=${sample_id}_Agg --csv=${sample_id}_Agg.csv --normalize=none --nosecondary
  """
}


// ---------------------- VDJ ----------------------
process makefastq_vdj {
  tag "${sample_id}"
  publishDir "${params.outdir}/VDJ/cellranger", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    path vdj_raw
    val sample_id
    path vdj_csv
    val cellranger

  output:
    path "${sample_id}_mkfastq_outs", emit: fastqs
    path "sample_ids.txt",           emit: sample_file

  shell:
  '''
  /datacommons/dhvi/10X_Genomics_data/cellranger/cellranger-!{cellranger}/cellranger mkfastq \
    --run=!{vdj_raw} --csv=!{vdj_csv} --output-dir=!{sample_id}_mkfastq_outs --delete-undetermined
  awk -F, 'NR > 1 {print $2}' !{vdj_csv} > sample_ids.txt
  '''
}

process vdj {
  tag "${sample}"
  publishDir "${params.outdir}/VDJ/cellranger", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    tuple val(sample), val(fastqs)
    val species
    val cellranger

  output:
    path "${sample}", emit: current_sample

  script:
    if (species == "human")
      """
      echo "Using human VDJ reference"
      vdj.sh ${fastqs} ${sample} ${params.human_vdj_ref} ${cellranger}
      """
    else if (species == "rhesus")
      """
      echo "Using rhesus VDJ reference and inner enrichment primers"
      vdj_rhesus.sh ${fastqs} ${sample} ${params.rhesus_vdj_ref} ${cellranger} ${params.inner_enrichment}
      """
    else if (species == "mouse")
      """
      echo "Using mouse VDJ reference"
      vdj.sh ${fastqs} ${sample} ${params.mouse_vdj_ref} ${cellranger}
      """
}


// ---------------------- Downstream VDJ ----------------------
process fxnl_1to1_clonal_analysis_human {
  publishDir "${params.outdir}/VDJ", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    val vdj_paths
    val sample_count

  output:
    path 'analysis', emit: analysis

  script:
    if ( sample_count == 1 )
      """
      echo "one lane cloanalyst human"
      echo ${vdj_paths} > vdj_paths.txt
      sed -i 's/[][]//g; s/,\\s*/,/g' vdj_paths.txt
      VDJ_FOLDER=\$(cat vdj_paths.txt)
      get_fxnl_1to1_clonal_analysis_human_CR7.1_single_lane.sh \$VDJ_FOLDER
      """
    else
      """
      echo "multiple lanes cloanalyst human"
      echo ${vdj_paths} > vdj_paths.txt
      sed -i 's/[][]//g; s/,\\s*/,/g' vdj_paths.txt
      tr ',' '\\n' < vdj_paths.txt > separated_paths.txt
      get_fxnl_1to1_clonal_analysis_human_CR7.1.sh \$PWD/separated_paths.txt
      """
}

process fxnl_1to1_clonal_analysis_rhesus {
  publishDir "${params.outdir}/VDJ", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    val vdj_paths
    val sample_count

  output:
    path 'analysis', emit: analysis

  script:
    if ( sample_count == 1 )
      """
      echo "one lane cloanalyst rhesus"
      echo ${vdj_paths} > vdj_paths.txt
      sed -i 's/[][]//g; s/,\\s*/,/g' vdj_paths.txt
      VDJ_FOLDER=\$(cat vdj_paths.txt)
      get_fxnl_1to1_clonal_analysis_rhesus_CR7.1_single_lane.sh \$VDJ_FOLDER
      """
    else
      """
      echo "multiple lanes cloanalyst rhesus"
      echo ${vdj_paths} > vdj_paths.txt
      sed -i 's/[][]//g; s/,\\s*/,/g' vdj_paths.txt
      tr ',' '\\n' < vdj_paths.txt > separated_paths.txt
      get_fxnl_1to1_clonal_analysis_rhesus_CR7.1.sh \$PWD/separated_paths.txt
      """
}

process fxnl_1to1_clonal_analysis_mouse {
  publishDir "${params.outdir}/VDJ", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    val vdj_paths
    val sample_count

  output:
    path 'analysis', emit: analysis

  script:
    if ( sample_count == 1 )
      """
      echo "one lane cloanalyst mouse"
      echo ${vdj_paths} > vdj_paths.txt
      sed -i 's/[][]//g; s/,\\s*/,/g' vdj_paths.txt
      VDJ_FOLDER=\$(cat vdj_paths.txt)
      get_fxnl_1to1_clonal_analysis_mouse_CR6_single_lane.sh \$VDJ_FOLDER
      """
    else
      """
      echo "multiple lanes cloanalyst mouse"
      echo ${vdj_paths} > vdj_paths.txt
      sed -i 's/[][]//g; s/,\\s*/,/g' vdj_paths.txt
      tr ',' '\\n' < vdj_paths.txt > separated_paths.txt
      get_fxnl_1to1_clonal_analysis_mouse_CR6.sh \$PWD/separated_paths.txt
      """
}

process partis_fxnl_1to1_clonal_analysis_human {
  conda '/hpc/group/dhvi/WieheLab/scripts/nextflow_pipelines/conda_envs/partis'
  publishDir "${params.outdir}/VDJ", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    val vdj_paths
    val sample_count
    val germline

  output:
    path 'partis_analysis', emit: partis_analysis

  script:
    if ( sample_count == 1 )
      """
      echo "one lane partis human"
      echo ${vdj_paths} > vdj_paths.txt
      sed -i 's/[][]//g; s/,\\s*/,/g' vdj_paths.txt
      VDJ_FOLDER=\$(cat vdj_paths.txt)
      partis_fxnl_1to1_clonal_analysis_human_single.sh \$VDJ_FOLDER ${germline}
      """
    else
      """
      echo "multiple lanes partis human"
      echo ${vdj_paths} > vdj_paths.txt
      sed -i 's/[][]//g; s/,\\s*/,/g' vdj_paths.txt
      tr ',' '\\n' < vdj_paths.txt > separated_paths.txt
      partis_fxnl_1to1_clonal_analysis_human.sh \$PWD/separated_paths.txt ${germline}
      """
}

process partis_fxnl_1to1_clonal_analysis_rhesus {
  conda '/hpc/group/dhvi/WieheLab/scripts/nextflow_pipelines/conda_envs/partis'
  publishDir "${params.outdir}/VDJ", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    val vdj_paths
    val sample_count
    val germline

  output:
    path 'partis_analysis', emit: partis_analysis

  when:
    sample_count > 0

  script:
    def single_script = (germline == 'kimdb') ? 'partis_kimdb_single.sh'
                                              : 'partis_fxnl_1to1_clonal_analysis_rhesus_light_single.sh'
    def multi_script  = (germline == 'kimdb') ? 'partis_kimdb.sh'
                                              : 'partis_fxnl_1to1_clonal_analysis_rhesus_light.sh'

    if ( sample_count == 1 ) {
      """
      set -euo pipefail
      echo "one lane partis (${germline})"
      echo ${vdj_paths} > vdj_paths.txt
      sed -i 's/[][]//g; s/,\\s*/,/g' vdj_paths.txt
      VDJ_FOLDER=\$(cat vdj_paths.txt)
      ${single_script} "\$VDJ_FOLDER"
      """
    } else {
      """
      set -euo pipefail
      echo "multiple lanes partis (${germline})"
      echo ${vdj_paths} > vdj_paths.txt
      sed -i 's/[][]//g; s/,\\s*/,/g' vdj_paths.txt
      tr ',' '\\n' < vdj_paths.txt > separated_paths.txt
      ${multi_script} "\$PWD/separated_paths.txt"
      """
    }
}

process partis_fxnl_1to1_clonal_analysis_mouse {
  conda '/hpc/group/dhvi/WieheLab/scripts/nextflow_pipelines/conda_envs/partis'
  publishDir "${params.outdir}/VDJ", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    val vdj_paths
    val sample_count
    val germline

  output:
    path 'partis_analysis', emit: partis_analysis

  script:
    if ( sample_count == 1 )
      """
      echo "one lane"
      echo ${vdj_paths} > vdj_paths.txt
      sed -i 's/[][]//g; s/,\\s*/,/g' vdj_paths.txt
      VDJ_FOLDER=\$(cat vdj_paths.txt)
      get_fxnl_1to1_clonal_analysis_mouse_CR6_single_lane.sh \$VDJ_FOLDER ${germline}
      """
    else
      """
      echo "multiple lanes"
      echo ${vdj_paths} > vdj_paths.txt
      sed -i 's/[][]//g; s/,\\s*/,/g' vdj_paths.txt
      tr ',' '\\n' < vdj_paths.txt > separated_paths.txt
      get_fxnl_1to1_clonal_analysis_mouse_CR6.sh \$PWD/separated_paths.txt ${germline}
      """
}

process merged_cloanalyst_partis {
  beforeScript 'module load R/4.1.1-rhel8'
  publishDir "${params.outdir}/VDJ", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    val vdj_path
    val partis_path

  output:
    path 'merged_clonal_analysis_with_partis.csv', emit: merged_csv

  script:
  """
  Rscript /hpc/group/dhvi/WieheLab/scripts/nextflow_pipelines/10x-nextflow-pipeline/bin/add_Partis_to_clonal_analysis.R \
     ${vdj_path[0]} ${partis_path[0]}
  """
}


// ---------------------- BEAM ----------------------
process makefastq_beam {
  tag "${sample_id}"
  publishDir "${params.outdir}/BEAM/cellranger", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    path beam_raw
    val sample_id
    path beam_csv
    val cellranger
    path analysis   // only to enforce dependency

  output:
    path "${sample_id}_mkfastq_outs", emit: fastqs
    path "sample_ids.txt",           emit: sample_file

  shell:
  '''
  /datacommons/dhvi/10X_Genomics_data/cellranger/cellranger-!{cellranger}/cellranger mkfastq \
    --run=!{beam_raw} --csv=!{beam_csv} --output-dir=!{sample_id}_mkfastq_outs --delete-undetermined
  awk -F, 'NR > 1 {print $2}' !{beam_csv} > sample_ids.txt
  '''
}

process count_beam {
  tag "${sample}"
  publishDir "${params.outdir}/BEAM/cellranger", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    tuple val(sample), val(fastqs)
    val species
    val cellranger
    path feature_ref

  output:
    path "${sample}", emit: current_sample

  script:
    def count_script = cellranger.startsWith("9") ? "count_beam_new.sh" : "count_beam.sh"
    if (species == "human")
      """
      echo "Using human transcriptome"
      ${count_script} ${fastqs} ${sample} ${params.human_transcriptome} ${cellranger} ${feature_ref}
      """
    else if (species == "rhesus")
      """
      echo "Using rhesus transcriptome"
      ${count_script} ${fastqs} ${sample} ${params.rhesus_transcriptome} ${cellranger} ${feature_ref}
      """
    else if (species == "mouse")
      """
      echo "Using mouse transcriptome"
      ${count_script} ${fastqs} ${sample} ${params.mouse_transcriptome} ${cellranger} ${feature_ref}
      """
}

process mat2csv {
  tag "${sample}"
  publishDir "${params.outdir}/BEAM/cellranger", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    path sample
    val cellranger

  output:
    path "${sample}.mat2csv.out.csv", emit: current_csv

  script:
  """
  echo ${sample} > sample.txt
  /datacommons/dhvi/10X_Genomics_data/cellranger/cellranger-${cellranger}/cellranger mat2csv \
    ${sample}/outs/raw_feature_bc_matrix ${sample}.mat2csv.out.csv
  """
}

process calc_beam_scores {
  beforeScript 'module load R/4.1.1-rhel8'
  publishDir "${params.outdir}/BEAM", mode: 'copy', overwrite: true
  cache 'lenient'

  input:
    val mat_paths
    path feature_ref

  output:
    path "*"

  script:
  """
  echo ${mat_paths} > mat_paths.txt
  sed -i 's/[][]//g; s/,\\s*/,/g' mat_paths.txt
  tr ',' '\\n' < mat_paths.txt > separated_paths.txt
  Rscript /hpc/group/dhvi/WieheLab/scripts/nextflow_pipelines/10x-nextflow-pipeline/bin/calculate_BEAM_scores_1.R \
     ${params.outdir} separated_paths.txt ${feature_ref} ${params.outdir}/VDJ ${params.exp} ${params.ptid} ${params.control}
  """
}


// ---------------------- WORKFLOW ----------------------
workflow {

  // Local copies
  def exp        = params.exp
  def ptid       = params.ptid
  def sample_id  = params.sample_id
  def wf         = params.workflow
  def species    = params.species
  def cellranger = params.cellranger
  def control    = params.control
  def clonal     = params.clonal

  // Optional paths -> channels
  if (params.feature_ref?.trim()) feature_ref = nextflow.Channel.fromPath(params.feature_ref)
  if (params.gex_raw?.trim())     gex_raw     = nextflow.Channel.fromPath(params.gex_raw)
  if (params.vdj_raw?.trim())     vdj_raw     = nextflow.Channel.fromPath(params.vdj_raw)
  if (params.beam_raw?.trim())    beam_raw    = nextflow.Channel.fromPath(params.beam_raw)
  if (params.gex_csv?.trim())     gex_csv     = nextflow.Channel.fromPath(params.gex_csv)
  if (params.vdj_csv?.trim())     vdj_csv     = nextflow.Channel.fromPath(params.vdj_csv)
  if (params.beam_csv?.trim())    beam_csv    = nextflow.Channel.fromPath(params.beam_csv)

  // ---------- GEX ----------
  if (wf == "gex") {

    def gexTuples
    def gexNames
    def sampleFileCh

    if (params.gex_fastq_map?.trim()) {
      (gexTuples, gexNames) = parse_map_to_tuples(params.gex_fastq_map)
      // sample_ids.txt for aggr
      write_sample_ids( gexNames.collect() )
      sampleFileCh = write_sample_ids.out.sample_file
    }
    else if (params.gex_raw?.trim() && params.gex_csv?.trim()) {
      makefastq_gex(gex_raw, sample_id, gex_csv, cellranger)
      def sampleIds = makefastq_gex.out.sample_file.splitText().map{ it.trim() }
      def fastqsDir = makefastq_gex.out.fastqs
      def tuples    = sampleIds.combine(fastqsDir).map { s, d -> tuple(s, d) }

      gexTuples = tuples
      gexNames  = tuples.map{ it[0] }
      sampleFileCh = makefastq_gex.out.sample_file
    }
    else {
      exit 1, "For workflow=gex supply either --gex_fastq_map OR (--gex_raw and --gex_csv)."
    }

    count_gex(gexTuples, species, cellranger)
    agg_gex(count_gex.out.current_sample.collect(), sampleFileCh, sample_id, cellranger)
  }

  // ---------- VDJ ----------
  if (wf == "vdj") {

    // Germline setup
    if (params.germline == '' && species == "human") {
      germline = "/datacommons/dhvi/partis_germlines/human/"
    }
    if (params.germline != '' && species == "human") {
      germline = nextflow.Channel.fromPath(params.germline)
    }
    if (species == "rhesus") {
      germline = (params.germline?.trim()) ? params.germline.trim() : "bw"
    }

    def vdjTuples
    def vdjNames

    if (params.vdj_fastq_map?.trim()) {
      (vdjTuples, vdjNames) = parse_map_to_tuples(params.vdj_fastq_map)
    }
    else if (params.vdj_raw?.trim() && params.vdj_csv?.trim()) {
      makefastq_vdj(vdj_raw, sample_id, vdj_csv, cellranger)
      def sampleIds = makefastq_vdj.out.sample_file.splitText().map{ it.trim() }
      def fastqsDir = makefastq_vdj.out.fastqs
      def tuples    = sampleIds.combine(fastqsDir).map { s, d -> tuple(s, d) }
      vdjTuples = tuples
      vdjNames  = tuples.map{ it[0] }
    }
    else {
      exit 1, "For workflow=vdj supply either --vdj_fastq_map OR (--vdj_raw and --vdj_csv)."
    }

    vdj(vdjTuples, species, cellranger)

    vdjNames.collect().map{ it.size() }.set{ sampleCount }

    if (species == "human" && clonal == "c")
      fxnl_1to1_clonal_analysis_human(vdj.out.current_sample.collect(), sampleCount)
    if (species == "rhesus" && clonal == "c")
      fxnl_1to1_clonal_analysis_rhesus(vdj.out.current_sample.collect(), sampleCount)
    if (species == "mouse" && clonal == "c")
      fxnl_1to1_clonal_analysis_mouse(vdj.out.current_sample.collect(), sampleCount)

    if (species == "human" && clonal == "p")
      partis_fxnl_1to1_clonal_analysis_human(vdj.out.current_sample.collect(), sampleCount, germline)
    if (species == "rhesus" && clonal == "p")
      partis_fxnl_1to1_clonal_analysis_rhesus(vdj.out.current_sample.collect(), sampleCount, germline)
    if (species == "mouse" && clonal == "p")
      partis_fxnl_1to1_clonal_analysis_mouse(vdj.out.current_sample.collect(), sampleCount, germline)

    if (species == "human" && clonal == "b") {
      fxnl_1to1_clonal_analysis_human(vdj.out.current_sample.collect(), sampleCount)
      partis_fxnl_1to1_clonal_analysis_human(vdj.out.current_sample.collect(), sampleCount, germline)
      merged_cloanalyst_partis(
        fxnl_1to1_clonal_analysis_human.out.analysis.collect(),
        partis_fxnl_1to1_clonal_analysis_human.out.partis_analysis.collect()
      )
    }
    if (species == "rhesus" && clonal == "b") {
      fxnl_1to1_clonal_analysis_rhesus(vdj.out.current_sample.collect(), sampleCount)
      partis_fxnl_1to1_clonal_analysis_rhesus(vdj.out.current_sample.collect(), sampleCount, germline)
      merged_cloanalyst_partis(
        fxnl_1to1_clonal_analysis_rhesus.out.analysis.collect(),
        partis_fxnl_1to1_clonal_analysis_rhesus.out.partis_analysis.collect()
      )
    }
    if (species == "mouse" && clonal == "b") {
      fxnl_1to1_clonal_analysis_mouse(vdj.out.current_sample.collect(), sampleCount)
      partis_fxnl_1to1_clonal_analysis_mouse(vdj.out.current_sample.collect(), sampleCount, germline)
      merged_cloanalyst_partis(
        fxnl_1to1_clonal_analysis_mouse.out.analysis.collect(),
        partis_fxnl_1to1_clonal_analysis_mouse.out.partis_analysis.collect()
      )
    }
  }

  // ---------- BEAM ----------
  if (wf == "beam") {

    if (!params.feature_ref?.trim())
      exit 1, "workflow=beam requires --feature_ref"

    // Germline setup (same as VDJ)
    if (params.germline == '' && species == "human") {
      germline = "/datacommons/dhvi/partis_germlines/human/"
    }
    if (params.germline != '' && species == "human") {
      germline = nextflow.Channel.fromPath(params.germline)
    }
    if (species == "rhesus") {
      germline = (params.germline?.trim()) ? params.germline.trim() : "bw"
    }

    // VDJ upstream (map or mkfastq)
    def vdjTuplesB
    def vdjNamesB

    if (params.vdj_fastq_map?.trim()) {
      (vdjTuplesB, vdjNamesB) = parse_map_to_tuples(params.vdj_fastq_map)
    }
    else if (params.vdj_raw?.trim() && params.vdj_csv?.trim()) {
      makefastq_vdj(vdj_raw, sample_id, vdj_csv, cellranger)
      def sampleIds = makefastq_vdj.out.sample_file.splitText().map{ it.trim() }
      def fastqsDir = makefastq_vdj.out.fastqs
      def tuples    = sampleIds.combine(fastqsDir).map { s, d -> tuple(s, d) }
      vdjTuplesB = tuples
      vdjNamesB  = tuples.map{ it[0] }
    }
    else {
      exit 1, "For workflow=beam, provide VDJ via --vdj_fastq_map OR (--vdj_raw and --vdj_csv)."
    }

    vdj(vdjTuplesB, species, cellranger)
    vdjNamesB.collect().map{ it.size() }.set{ sampleCount }

    // Branch choose (to enforce BEAM mkfastq dependency if needed)
    def analysisWait = null
    if (species == "human" && clonal == "c") {
      fxnl_1to1_clonal_analysis_human(vdj.out.current_sample.collect(), sampleCount)
      analysisWait = fxnl_1to1_clonal_analysis_human.out.analysis
    }
    if (species == "rhesus" && clonal == "c") {
      fxnl_1to1_clonal_analysis_rhesus(vdj.out.current_sample.collect(), sampleCount)
      analysisWait = fxnl_1to1_clonal_analysis_rhesus.out.analysis
    }
    if (species == "mouse" && clonal == "c") {
      fxnl_1to1_clonal_analysis_mouse(vdj.out.current_sample.collect(), sampleCount)
      analysisWait = fxnl_1to1_clonal_analysis_mouse.out.analysis
    }
    if (species == "human" && clonal == "p") {
      partis_fxnl_1to1_clonal_analysis_human(vdj.out.current_sample.collect(), sampleCount, germline)
      analysisWait = partis_fxnl_1to1_clonal_analysis_human.out.partis_analysis
    }
    if (species == "rhesus" && clonal == "p") {
      partis_fxnl_1to1_clonal_analysis_rhesus(vdj.out.current_sample.collect(), sampleCount, germline)
      analysisWait = partis_fxnl_1to1_clonal_analysis_rhesus.out.partis_analysis
    }
    if (species == "mouse" && clonal == "p") {
      partis_fxnl_1to1_clonal_analysis_mouse(vdj.out.current_sample.collect(), sampleCount, germline)
      analysisWait = partis_fxnl_1to1_clonal_analysis_mouse.out.partis_analysis
    }
    if (species == "human" && clonal == "b") {
      fxnl_1to1_clonal_analysis_human(vdj.out.current_sample.collect(), sampleCount)
      partis_fxnl_1to1_clonal_analysis_human(vdj.out.current_sample.collect(), sampleCount, germline)
      merged_cloanalyst_partis(
        fxnl_1to1_clonal_analysis_human.out.analysis.collect(),
        partis_fxnl_1to1_clonal_analysis_human.out.partis_analysis.collect()
      )
      analysisWait = merged_cloanalyst_partis.out.merged_csv
    }
    if (species == "rhesus" && clonal == "b") {
      fxnl_1to1_clonal_analysis_rhesus(vdj.out.current_sample.collect(), sampleCount)
      partis_fxnl_1to1_clonal_analysis_rhesus(vdj.out.current_sample.collect(), sampleCount, germline)
      merged_cloanalyst_partis(
        fxnl_1to1_clonal_analysis_rhesus.out.analysis.collect(),
        partis_fxnl_1to1_clonal_analysis_rhesus.out.partis_analysis.collect()
      )
      analysisWait = merged_cloanalyst_partis.out.merged_csv
    }
    if (species == "mouse" && clonal == "b") {
      fxnl_1to1_clonal_analysis_mouse(vdj.out.current_sample.collect(), sampleCount)
      partis_fxnl_1to1_clonal_analysis_mouse(vdj.out.current_sample.collect(), sampleCount, germline)
      merged_cloanalyst_partis(
        fxnl_1to1_clonal_analysis_mouse.out.analysis.collect(),
        partis_fxnl_1to1_clonal_analysis_mouse.out.partis_analysis.collect()
      )
      analysisWait = merged_cloanalyst_partis.out.merged_csv
    }

    // BEAM FASTQs (map or mkfastq)
    def beamTuples
    def beamNames

    if (params.beam_fastq_map?.trim()) {
      (beamTuples, beamNames) = parse_map_to_tuples(params.beam_fastq_map)
    }
    else if (params.beam_raw?.trim() && params.beam_csv?.trim()) {
      if (analysisWait == null)
        exit 1, "BEAM mkfastq route requires an upstream analysis (choose c/p/b)."
      makefastq_beam(beam_raw, sample_id, beam_csv, cellranger, analysisWait)
      def beamIds = makefastq_beam.out.sample_file.splitText().map{ it.trim() }
      def beamDir = makefastq_beam.out.fastqs
      def tuples  = beamIds.combine(beamDir).map { s, d -> tuple(s, d) }
      beamTuples = tuples
      beamNames  = tuples.map{ it[0] }
    }
    else {
      exit 1, "For workflow=beam supply --beam_fastq_map OR (--beam_raw and --beam_csv)."
    }

    // Run BEAM
    count_beam(beamTuples, species, cellranger, feature_ref)
    mat2csv(count_beam.out.current_sample, cellranger)
    calc_beam_scores(mat2csv.out.current_csv.collect(), feature_ref)
  }
}


