#!/usr/bin/env bash

project_dir="$($(dirname $0)/abs_path.py $(dirname $0)/../)"
export __PIPELINE_PROJECT_DIR__=${project_dir}
export PATH="${project_dir}/bin:$PATH"
. "${project_dir}/bin/load_tools_paths.sh"

. $(dirname $0)/get_common_paths.sh

logs_dir="${project_dir}/logs"
data_dir="${project_dir}/data"
__dir="${project_dir}/.pipeline"
checkpoints_dir="${__dir}/checkpoints"
dbs_automatic_dir=${dbs_dir}/automatic
cfg_file=${1-${dbs_dir}/analysis_configuration}

# Runs commands in parallel
par_tsv () {
  parallel --colsep '\t' '--header' ':' "${@}"
}

# Returns 1 if there is the checkpoint for the step with the name passed as
#   parameter.
is_step_to_run () {
  local cur_step=${1}
  local res=0
  if [[ -e "${checkpoints_dir}/${cur_step}" ]]
  then
    res=1
  fi
  return $res
}

# Checks the exit status of the current step and creates the checkpoint if
#   the step was succesful. It will exit with an error otherwise.
check_status_and_create_checkpoint () {
  local cur_step=${1}
  local exit_status=${2}
  local checkpoint_file="${checkpoints_dir}/${cur_step}"
  if [[ ! -e $checkpoint_file ]]
  then
    if (( exit_status == 0 ))
    then
      touch $checkpoint_file
    else
      echo "Error at step ${step} exiting"
      exit 1
    fi
  fi
}

# Saves the data abount the step that was run. The saved data are:
#   - name of the step
#   - the entire command that was executed
#   - the date and time at which the tool was started
#   - the execution duration of the tool in seconds
#   - the exit code returned by the tool
save_tool_metadata () {
  local step_name=${1}
  local command_line="${2}"
  local start_time="${3}"
  local duration="${4}"
  local exit_status=${5}
  local tool_metadata_file="${__dir}/tool_metadata.tsv"
  if [[ ! -e $tool_metadata_file ]]
  then
    echo -e "step\tcommand_line\tstart_time\tduration\texit_status" > $tool_metadata_file
  fi
  echo -e "${step_name}\t${command_line}\t${start_time}\t${duration}\t${exit_status}" >> $tool_metadata_file
}


# Runs each step of the pipeline. This function create the folders for the step
# if the create_folder parameter is set to true.
# Then the command that is passed is run, the metadata is saved and the pipeline
# proceeds the execution if the step is completed succesfully.
_run_step () {
  local step=${1}; shift
  local create_folder=${1}; shift

  if is_step_to_run $step
  then
    if (( create_folder ))
    then
      mk_step.sh $step
    fi
    # echo "internal cmd:" ${@}
    # echo "cmd:" $(eval echo ${@})
    start_seconds=$(date +%s)
    start_time=$(date +%Y%m%d%H%M%S)
    cmd="$(eval echo "${@}")"
    $(eval echo "${@}")
    local exit_status=${?}
    # echo ${exit_status}
    duration=$(( $(date +%s) - ${start_seconds} ))
    save_tool_metadata "${step}" "${cmd}" "${start_time}" "${duration}" "${exit_status}"
    check_status_and_create_checkpoint $step $exit_status
  fi
}

# Execute a step of the pipeline.
run_step () {
  local step=${1}; shift
  _run_step $step 1 ${@}
}

# Execute a step of the pipeline without creating the folder.
run_step_no_dir () {
  local step=${1}; shift
  _run_step $step 0 ${@}
}

if [[ -e "$cfg_file" ]]
then
  . "$cfg_file"
else
  echo "Required configuration file: \"${cfg_file}\" doesn't exist"
  exit 1
fi

all_vars_defined=1
for v in project_name tools_dir resources_dir
do
  if eval "[[ ! \$$v ]]"
  then
    if (( $all_vars_defined ))
    then
      echo "ERROR:"
    fi
    echo -e "\tConfiguration variable \"$v\" must have a value"
    all_vars_defined=0
  fi
done

(( $all_vars_defined )) || exit 1

unset all_vars_defined

sif=${sif-${dbs_dir}/sample_pairs.tsv}
sample_paths=${sample_paths-${dbs_dir}/sample_paths.tsv}

segmentation=${segmentation-cnvkit}
use_clonet_auto=${use_clonet_auto-0}
run_data_aggregation=${run_data_aggregation-0}
create_panel_of_normal=${create_panel_of_normal-0}

sif_auto=${dbs_automatic_dir}/samples.tsv
sif_auto_single=${dbs_automatic_dir}/samples_single.tsv
sif_logs_auto_dir="${project_dir}/logs/auto"

samples_dir="${data_dir}/samples"
normal_id="{analysis_id}-n"
normal_bam="${samples_dir}/${normal_id}.bam"
tumor_id="{analysis_id}"
tumor_bam="${samples_dir}/${tumor_id}.bam"
single_bam="${samples_dir}/{analysis_id}.bam"
analysis_id="$tumor_id"
step_data_dir="${data_dir}/\${step}"
step_log_dir="${logs_dir}/\${step}"
data_logs_dirs="\"${step_data_dir}\" \"${step_log_dir}\""

threads_many=${threads_many-40}
threads_mid=${threads_mid-10}
threads_few=${threads_few-5}
threads_per_sample_many=${threads_per_sample_many-5}
threads_per_sample_few=${threads_per_sample_few-3}

for d in "${dbs_automatic_dir}" "${sif_logs_auto_dir}" "${samples_dir}" "${checkpoints_dir}"
do
  [[ -e "${d}" ]] || mkdir -p "${d}"
done

for d in "tools" "resources"
do
  cur_dir="${project_dir}/${d}"
  [[ -e "${cur_dir}" ]] || ln -s "$(eval echo \$${d}'_dir')" "${cur_dir}"
done

. ${bin_dir}/load_tools_paths.sh

[[ -e "$sif_auto"        ]] || generate_sample_ids.sh "$sif" "$sif_auto" "${logs_dir}/auto/generate_ids.log"

[[ -e "$sif_auto_single" ]] || generate_sample_info_file_single.sh "$sif_auto" "$sif_auto_single"

export _JAVA_OPTIONS='-Xmx8g -XX:ParallelGCThreads=1'
export SINGULARITYENV__JAVA_OPTIONS='-Xmx8g -XX:ParallelGCThreads=1'



# This step creates a folder containing the links to the input BAM files.
run_step 'samples'      \
  par_tsv -q            \
  link_sample.sh        \
  "{file}"              \
  "$single_bam"         \
  :::: $sif_auto_single


# This step creates the links to the BAI files if the files are found along to
#   the BAM files. Otherwise the indices are created using samtools.
run_step 'samples_indices'     \
  par_tsv -j$threads_few       \
  run_samtools_index_single.sh \
    "$data_logs_dirs"          \
    "$analysis_id"             \
    "$single_bam"              \
  :::: $sif_auto_single


# Computes the pileup of all the SNPs positions available in the DBSnp version
#   passed as input.
run_step 'snps_pileup'                   \
  par_tsv -j$threads_few                 \
  ionice -c3                             \
  run_pacbam_single.sh                   \
  --threads "${threads_per_sample_many}" \
  --kit_name "{kit_name}"                \
  --dbsnp_version "{dbsnp_version}"      \
  "$data_logs_dirs"                      \
  "$analysis_id"                         \
  "$single_bam"                          \
  "{reference_genome}"                   \
  :::: $sif_auto_single


# Runs the picard CollectHsMetrics on the BAM files to collect the sequencing
#   statistics for the samples.
run_step 'picard_hsmetrics'    \
  par_tsv -j$threads_mid       \
  ionice -c3                   \
  run_picard_metrics_single.sh \
  "$data_logs_dirs"            \
  "$analysis_id"               \
  "$single_bam"                \
  "{reference_genome}"         \
  "{kit_name}"                 \
  :::: $sif_auto_single


if (( ${create_panel_of_normal} ))
then

  echo "Creation of the panel of normal."

  # If enabled with the specific option runs the mutect PON on the normal BAMs
  #   to collect the data necessary to construct a Panel of Normals.
  run_step 'mutect2_pon'                                          \
    par_tsv -j$threads_many                                       \
    run_mutect2_pon_single.sh                                     \
    "$data_logs_dirs"                                             \
    "${normal_id}"                                                \
    "$normal_bam"                                                 \
    "{reference_genome}"                                          \
    "{kit_name}"                                                  \
    "{dbsnp_version}"                                             \
    "{cosmic_version}"                                            \
    :::: $sif_auto

fi


# Prepares the input files to run the SPIA tool starting from the pileup
#   computed in a preceding step.
run_step 'spia_genotypes'                                                                         \
  par_tsv -j$threads_many                                                                         \
  generate_genotype_vcf_single.sh                                                                 \
  --snps_list_file "${res_dir}/spia/{reference_genome}/snps_spia_default-{reference_genome}.vcf" \
  "$data_logs_dirs"                                                                               \
  "$analysis_id"                                                                                  \
  "${data_dir}/snps_pileup/${analysis_id}/${analysis_id}.snps"                                    \
  :::: $sif_auto_single


# Runs SPIA: a tool that checks the genotype distance (fraction of SNP sites that
#   are different between the paired samples) to assess the proper pairing of
#   matched samples.
run_step 'spia'                                    \
  par_tsv -j$threads_many                          \
  run_spia_single.sh                               \
  "$data_logs_dirs"                                \
  "$analysis_id"                                   \
  "${data_dir}/spia_genotypes/{analysis_id}-n.vcf" \
  "${data_dir}/spia_genotypes/{analysis_id}.vcf"   \
  :::: $sif_auto


# Prepares the input files to run the EthSEQ tool starting from the pileup
#   computed in a preceding step.
run_step 'ethseq_genotypes'                                                                                     \
  par_tsv -j$threads_many                                                                                       \
  generate_genotype_vcf_single.sh                                                                               \
  --snps_list_file "${res_dir}/ethseq/{reference_genome}/ethseq-universal_exonic_model-{reference_genome}.vcf" \
  "$data_logs_dirs"                                                                                             \
  "$analysis_id"                                                                                                \
  "${data_dir}/snps_pileup/${analysis_id}/${analysis_id}.snps"                                                  \
  :::: $sif_auto_single


# Runs EthSEQ: a tool that is able to infer the ethnicity of a sample.
run_step 'ethseq'                                                                             \
  par_tsv -j$threads_few                                                                      \
  ethseq ${project_dir}/bin/run_ethseq_single.R                                               \
  --log_file "${step_log_dir}/${analysis_id}.log"                                             \
  ${ethseq_snps_model}                                                                        \
  "${res_dir}/ethseq/{reference_genome}/ethseq-universal_exonic_model-{reference_genome}.gds" \
  "${data_dir}/ethseq_genotypes/${analysis_id}.vcf"                                           \
  "${data_dir}/\${step}/${analysis_id}/"                                                      \
  :::: $sif_auto_single


# Runs MuTect to call SNVs and indels.
run_step 'mutect2'                                              \
  par_tsv -j$threads_many                                       \
  run_mutect2_single.sh                                         \
  "$data_logs_dirs"                                             \
  "$analysis_id"                                                \
  "$normal_bam"                                                 \
  "$tumor_bam"                                                  \
  "{reference_genome}"                                          \
  "{kit_name}"                                                  \
  "{dbsnp_version}"                                             \
  "{cosmic_version}"                                            \
  :::: $sif_auto


# Filters SNVs calls that were annotated as PASS by mutect2.
run_step 'mutect2_filter'                \
  par_tsv -j$threads_many                \
  run_mutect2_filter_single.sh           \
  "$data_logs_dirs"                      \
  "$analysis_id"                         \
  "$data_dir/mutect2/${analysis_id}.vcf" \
  :::: $sif_auto


# Creates a file containing only SNVs calls (removes indel calls)
run_step 'mutect2_filter_snvs'                  \
  par_tsv -j$threads_many                       \
  run_mutect2_filter_snvs.sh                    \
  "$data_logs_dirs"                             \
  "$analysis_id"                                \
  "$data_dir/mutect2_filter/${analysis_id}.vcf" \
  :::: $sif_auto


# Creates a bed file with the positions where the SNVs were found.
run_step 'snv_bed'                                   \
  par_tsv -j$threads_mid                             \
  generate_bed_from_vcf_single.sh                    \
  "$data_logs_dirs"                                  \
  "$analysis_id"                                     \
  "$data_dir/mutect2_filter_snvs/${analysis_id}.vcf" \
  :::: $sif_auto


# Annotates the variants found by mutect.
run_step 'vep'                                    \
  par_tsv -j$threads_mid                          \
  run_vep_single.sh                               \
  "$data_logs_dirs"                               \
  "$analysis_id"                                  \
  "${data_dir}/mutect2_filter/${analysis_id}.vcf" \
  "{reference_genome}"                            \
  :::: $sif_auto


# Filter the pileup of the tumor and the normal sample in parallel
#   to find the SNPs that are heterozygous (in the normal sample) and have
#   sufficient coverage in both samples.
run_step 'informative_snps'                                \
  par_tsv -j$threads_mid                                   \
  ionice -c3                                               \
  generate_informative_snps.sh                             \
  "$data_logs_dirs"                                        \
  "$analysis_id"                                           \
  "${data_dir}/snps_pileup/${normal_id}/${normal_id}.snps" \
  "${data_dir}/snps_pileup/${tumor_id}/${tumor_id}.snps"   \
  :::: $sif_auto


# Prepares the input files to run FACETS
run_step 'facets_input'                                    \
  par_tsv -j$threads_mid                                   \
  ionice -c3                                               \
  generate_facets_input.sh                                 \
  "$data_logs_dirs"                                        \
  "$analysis_id"                                           \
  "${data_dir}/snps_pileup/${normal_id}/${normal_id}.snps" \
  "${data_dir}/snps_pileup/${tumor_id}/${tumor_id}.snps"   \
  :::: $sif_auto


# Runs FACETS and extracts the copy number segmentation.
run_step 'segfiles_facets'                          \
  par_tsv -j$threads_many                           \
  facets ${project_dir}/bin/segment_facets_single.R \
  --reference_genome {reference_genome}             \
  --log_file "${step_log_dir}/${analysis_id}.log"   \
  --min_depth 10                                    \
  "${data_dir}/facets_input/${analysis_id}.txt"     \
  "${data_dir}/\${step}/${analysis_id}.seg"         \
  :::: $sif_auto


# Runs CNVKit in paired mode to obtain the copy number segmentation for the
#   samples.
run_step 'cnvkit'         \
  par_tsv -j$threads_many \
  run_cnvkit_single.sh    \
  --sex "{sex}"           \
  "$data_logs_dirs"       \
  "$analysis_id"          \
  "$normal_id"            \
  "$normal_bam"           \
  "$tumor_id"             \
  "$tumor_bam"            \
  "{reference_genome}"    \
  "{kit_name}"            \
  :::: $sif_auto


# Extract the file containing the segments generated from the CNVKit run data.
run_step 'segfiles_cnvkit'                               \
  par_tsv -j$threads_many                                \
  run_cnvkit_export_segfile_single.sh                    \
  "$data_logs_dirs"                                      \
  "$analysis_id"                                         \
  "${data_dir}/cnvkit/${analysis_id}/${analysis_id}.cns" \
  :::: $sif_auto


# Runs a HMM based segmentation (see SLMSuite).
run_step 'cnvkit_slm'                                    \
  par_tsv -j$threads_many                                \
  run_cnvkit_segment_single.sh                           \
  --method slm                                           \
  "$data_logs_dirs"                                      \
  "$analysis_id"                                         \
  "${data_dir}/cnvkit/${analysis_id}/${analysis_id}.cnr" \
  :::: $sif_auto


# Extract the file containing the segments generated from the CNVKit SLM run
#   data.
run_step 'segfiles_slm'                       \
  par_tsv -j$threads_many                     \
  run_cnvkit_export_segfile_single.sh         \
  "$data_logs_dirs"                           \
  "$analysis_id"                              \
  "${data_dir}/cnvkit_slm/${analysis_id}.cns" \
  :::: $sif_auto


case "$segmentation" in
  "cnvkit")
    selected_segfiles="segfiles_cnvkit"
    ;;
  "slm")
    selected_segfiles="segfiles_slm"
    ;;
  "facets")
    selected_segfiles="segfiles_facets"
    ;;
  *)
    echo "Segmentation '$segmentation' is not supported."
    exit 1
    ;;
esac
segfiles_dir="${data_dir}/segfiles"
if [[ ! ( -h "$segfiles_dir" || -e "$segfiles_dir" ) ]]
then
  ln -s "${selected_segfiles}" "$segfiles_dir"
else
  echo "WARNING: link to the segfiles already exists. If you want the pipeline to update the link you should delete it."
fi


# Computes the pileup at the SNVs positions identified by MuTect.
run_step 'snv_read_counts'                            \
  par_tsv -j$threads_few                              \
  ionice -c3                                          \
  run_pacbam_single.sh                                \
  --threads "${threads_per_sample_many}"              \
  --bed "$data_dir/snv_bed/{pair_id}.bed"             \
  --vcf "$data_dir/mutect2_filter_snvs/{pair_id}.vcf" \
  "$data_logs_dirs"                                   \
  "{analysis_id}"                                     \
  "$single_bam"                                       \
  "{reference_genome}"                                \
  :::: $sif_auto_single


# Integrates the pileup informations with the variant annotated by VEP.
run_step 'merge_vep_pacbam'                                    \
  par_tsv -j$threads_many                                      \
  generate_merged_vep_pacbam.sh                                \
  "$data_logs_dirs"                                            \
  "$analysis_id"                                               \
  "${data_dir}/vep/${analysis_id}.tsv"                         \
  "${data_dir}/snv_read_counts/${normal_id}/${normal_id}.snps" \
  "${data_dir}/snv_read_counts/${tumor_id}/${tumor_id}.snps"   \
  :::: $sif_auto


# Prepares the input files to run clonet with the admixture inference method
#   based on the SNVs.
run_step 'snv_read_counts_clonet'                            \
  par_tsv -j$threads_many                                    \
  generate_read_count_single.sh                              \
  "$data_logs_dirs"                                          \
  "$analysis_id"                                             \
  "${data_dir}/snv_read_counts/${tumor_id}/${tumor_id}.snps" \
  :::: $sif_auto


# Runs CLONET: a tool that is able to mitigate the effects that ploidy and
#   admixture have on copy number.
run_step 'clonet_auto'                        \
  par_tsv -j$threads_mid                      \
  run_clonet_single.sh                        \
  --threads $threads_per_sample_many          \
  "$data_logs_dirs"                           \
  "$analysis_id"                              \
  "$normal_id"                                \
  "$normal_bam"                               \
  "$tumor_id"                                 \
  "$tumor_bam"                                \
  "${data_dir}/informative_snps/$analysis_id" \
  "${data_dir}/segfiles/${analysis_id}.seg"   \
  :::: $sif_auto


# Runs CLONET-SNV: a tool that is able to mitigate the effects that ploidy and
#   admixture have on copy number. (This version of the tool estimate the
#   admixture using the SNVs)
run_step 'clonet_snv_auto'                                \
  par_tsv -j$threads_mid                                  \
  run_clonet_snv_single.sh                                \
  --threads $threads_per_sample_many                      \
  "$data_logs_dirs"                                       \
  "$analysis_id"                                          \
  "$normal_id"                                            \
  "$normal_bam"                                           \
  "$tumor_id"                                             \
  "$tumor_bam"                                            \
  "${data_dir}/informative_snps/$analysis_id"             \
  "${data_dir}/segfiles/${analysis_id}.seg"               \
  "${data_dir}/snv_read_counts_clonet/${analysis_id}.txt" \
  :::: $sif_auto


clonet_postfix_dir='manual'
if (( use_clonet_auto ))
then
  clonet_postfix_dir='auto'
fi

for n in 'clonet' 'clonet_snv'
do
  eval "have_${n}=0"
  selected_clonet="${n}_${clonet_postfix_dir}"
  clonet_dir="${data_dir}/${n}"
  if [[ ! -h "${clonet_dir}" ]]
  then
    ln -s "$selected_clonet" "${clonet_dir}"
  fi
  if [[ -d "${clonet_dir}" && ( ${use_clonet_auto} -eq 0 || -e "${checkpoints_dir}/${selected_clonet}" ) ]]
  then
    eval "have_${n}=1"
    echo "WARNING: link to the ${n} folder already exists. If you want the pipeline to update the link you should delete it."
  fi
done
unset clonet_dir clonet_selected_dir

if (( have_clonet ))
then

  # Computes corrected AF and clonality for the SNVs based on the ploidy and
  #   admixture estimates obtained from CLONET.
  run_step 'snv_clonality_clonet'                                         \
    par_tsv -j$threads_many                                               \
    run_clonet_tool.sh ${project_dir}/bin/run_snv_clonality.R             \
    --log_file "${step_log_dir}/${analysis_id}.log"                       \
    "/usr/share/clonet/src"                                               \
    "$normal_id"                                                          \
    "$tumor_id"                                                           \
    "${data_dir}/clonet/${analysis_id}/Results/allelicImbalanceTable.txt" \
    "${data_dir}/merge_vep_pacbam/${analysis_id}.txt"                     \
    "$step_data_dir"                                                      \
    :::: $sif_auto

fi


echo -e "Pipeline successfully completed."

