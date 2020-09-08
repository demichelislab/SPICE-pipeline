project_dir="$( abs_path.py $(dirname $0)/../ )"

function check_file_exists_or_exit {
  if [[ ! -e $1 ]]
  then
    echo -e "$2"
    exit 1
  fi
}

bin_dir="${project_dir}/bin"
dbs_dir="${project_dir}/dbs"
cfg_dir="${project_dir}/cfg"
res_dir="${project_dir}/resources"

if [[ ! -d ${bin_dir} ]]
then
  echo "Trying to load get_common_paths.sh from the wrong place."
fi

is_reference_genome_specified=0
is_kit_specified=0

if [[ ${reference_genome_name} ]]
then
  reference_genome_fasta="${res_dir}/reference_genome/${reference_genome_name}/reference_genome-${reference_genome_name}.fasta"
  reference_genome_chromosome_size="${res_dir}/reference_genome/${reference_genome_name}/chromosome_size-${reference_genome_name}.tsv"
  check_file_exists_or_exit "${reference_genome_fasta}" \
                            "Reference genome ${reference_genome_name} does not exists.\nCannot find file ${reference_genome_fasta}"
  check_file_exists_or_exit "${reference_genome_chromosome_size}" \
                            "Chromosome sizes file for ${reference_genome_name} does not exists.\nCannot find file ${reference_genome_chromosome_size}"
  is_reference_genome_specified=1
fi

if [[ ${is_reference_genome_specified} && ${kit_name} ]]
then
  is_kit_specified=1
  kit_target_bed="${res_dir}/kit/${reference_genome_name}/kit-${kit_name}-target-${reference_genome_name}.bed"
  kit_bait_bed="${res_dir}/kit/${reference_genome_name}/kit-${kit_name}-bait-${reference_genome_name}.bed"
  kit_target_interval_list="${res_dir}/kit/${reference_genome_name}/kit-${kit_name}-target-${reference_genome_name}.interval_list"
  kit_bait_interval_list="${res_dir}/kit/${reference_genome_name}/kit-${kit_name}-bait-${reference_genome_name}.interval_list"
  check_file_exists_or_exit "${kit_target_bed}" "Kit ${kit_name} target does not exists.\nCannot find file ${kit_target_bed}"
  check_file_exists_or_exit "${kit_bait_bed}" "Kit ${kit_name} bait does not exists.\nCannot find file ${kit_bait_bed}"
  check_file_exists_or_exit "${kit_target_interval_list}" "Interval file for kit ${kit_name} target does not exists.\nCannot find file ${kit_target_interval_list}"
  check_file_exists_or_exit "${kit_bait_interval_list}" "Interval file for kit ${kit_name} bait does not exists.\nCannot find file ${kit_bait_interval_list}"
fi

is_dbsnp_version_specified=0
if [[ ${is_reference_genome_specified} && ${dbsnp_version} ]]
then
  is_dbsnp_version_specified=1
  dbsnp_base_path="${res_dir}/dbsnp/${reference_genome_name}/${dbsnp_version}"
  dbsnp_full_vcf="${dbsnp_base_path}/dbsnp_full-${dbsnp_version}-${reference_genome_name}.vcf.gz"
  dbsnp_full_vcf_decompressed="${dbsnp_base_path}/dbsnp_full-${dbsnp_version}-${reference_genome_name}.vcf"
  check_file_exists_or_exit "${dbsnp_full_vcf}" \
                            "Dbsnp version ${dbsnp_version} for reference genome ${reference_genome_name} does not exists.\nCannot find file ${dbsnp_full_vcf}"
  if [[ ${is_kit_specified} ]]
  then
    dbsnp_restricted_kit_target_vcf="${dbsnp_base_path}/restricted_kit/${kit_name}/dbsnp_restricted_kit-${kit_name}-target-${dbsnp_version}-${reference_genome_name}.vcf"
    dbsnp_restricted_kit_bait_vcf="${dbsnp_base_path}/restricted_kit/${kit_name}/dbsnp_restricted_kit-${kit_name}-bait-${dbsnp_version}-${reference_genome_name}.vcf"
    dbsnp_restricted_kit_target_only_snps_vcf="${dbsnp_base_path}/restricted_kit/${kit_name}/dbsnp_restricted_kit_onlysnp-${kit_name}-target-${dbsnp_version}-${reference_genome_name}.vcf"
    dbsnp_restricted_kit_bait_only_snps_vcf="${dbsnp_base_path}/restricted_kit/${kit_name}/dbsnp_restricted_kit_onlysnp-${kit_name}-bait-${dbsnp_version}-${reference_genome_name}.vcf"
    dbsnp_restricted_kit_target_only_snps_onealt_vcf="${dbsnp_base_path}/restricted_kit/${kit_name}/dbsnp_restricted_kit_onlysnp_onealt-${kit_name}-target-${dbsnp_version}-${reference_genome_name}.vcf"
    dbsnp_restricted_kit_bait_only_snps_onealt_vcf="${dbsnp_base_path}/restricted_kit/${kit_name}/dbsnp_restricted_kit_onlysnp_onealt-${kit_name}-bait-${dbsnp_version}-${reference_genome_name}.vcf"
    check_file_exists_or_exit "${dbsnp_restricted_kit_target_vcf}" \
                              "Dbsnp version ${dbsnp_version} for reference genome ${reference_genome_name} restricted to kit ${kit_name} target regions does not exists.\nCannot find file ${dbsnp_restricted_kit_target_vcf}"
    check_file_exists_or_exit "${dbsnp_restricted_kit_bait_vcf}" \
                              "Dbsnp version ${dbsnp_version} for reference genome ${reference_genome_name} restricte to kit ${kit_name} bait regions does not exists.\nCannot find file ${dbsnp_restricted_kit_bait_vcf}"
    check_file_exists_or_exit "${dbsnp_restricted_kit_target_only_snps_vcf}" \
                              "Dbsnp version ${dbsnp_version} for reference genome ${reference_genome_name} restricted to kit ${kit_name} target regions does not exists.\nCannot find file ${dbsnp_restricted_kit_target_only_snps_vcf}"
    check_file_exists_or_exit "${dbsnp_restricted_kit_bait_only_snps_vcf}" \
                              "Dbsnp version ${dbsnp_version} for reference genome ${reference_genome_name} restricte to kit ${kit_name} bait regions does not exists.\nCannot find file ${dbsnp_restricted_kit_bait_only_snps_vcf}"
    check_file_exists_or_exit "${dbsnp_restricted_kit_target_only_snps_onealt_vcf}" \
                              "Dbsnp version ${dbsnp_version} for reference genome ${reference_genome_name} restricted to kit ${kit_name} target regions does not exists.\nCannot find file ${dbsnp_restricted_kit_target_only_snps_onealt_vcf}"
    check_file_exists_or_exit "${dbsnp_restricted_kit_bait_only_snps_onealt_vcf}" \
                              "Dbsnp version ${dbsnp_version} for reference genome ${reference_genome_name} restricte to kit ${kit_name} bait regions does not exists.\nCannot find file ${dbsnp_restricted_kit_bait_only_snps_onealt_vcf}"
  fi
fi 

# is_cosmic_version_specified=0
# if [[ ${is_reference_genome_specified} && ${cosmic_version} ]]
# then
#   is_cosmic_version_specified=0
#   cosmic_vcf="${res_dir}/cosmic/${reference_genome_name}/${cosmic_version}/cosmic-${cosmic_version}-${reference_genome_name}.vcf"
#   check_file_exists_or_exit "${cosmic_vcf}" \
#                             "Cosmic version ${cosmic_version} for reference genome ${reference_genome_name} does not exists.\nCannot find file ${cosmic_vcf}"
# fi

if [[ ${analysis_name} ]]
then
  tool_res_dir="${res_dir}/${analysis_name}"
fi

