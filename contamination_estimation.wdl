version development

#import "util.wdl" as util
import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/util.wdl" as util
#import "PileupSummaries.wdl" as ps
import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/PileupSummaries.wdl" as ps


workflow ContaminationEstimation {
    input {
        File? interval_list
        Array[File]? scattered_interval_list
        File tumor_bam
        File tumor_bai
        String tumor_sample_name
        File? normal_bam
        File? normal_bai
        String? normal_sample_name
        File ref_dict

        File? variants
        File? variants_idx
        File? vcf
        File? vcf_idx
        String? getpileupsummaries_extra_args

        Float minimum_population_allele_frequency = 0.01
        Float maximum_population_allele_frequency = 0.2

        Runtime? runtime_params

        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        Int emergency_extra_diskGB = 0

        Int mem_get_pileup_summaries = 4096  # needs at least 2G
        Int mem_gather_pileup_summaries = 512  # 64
        Int mem_calculate_contamination = 8192  # depends on the variants_for_contamination resource
        Int time_startup = 10
        Int time_get_pileup_summaries = 4500  # 3 d / scatter_count
        Int time_gather_pileup_summaries = 5
        Int time_calculate_contamination = 10
    }

    Int scatter_count = if defined(scattered_interval_list) then length(scattered_interval_list) else 1
    Int gatk_override_size = if defined(gatk_override) then ceil(size(gatk_override, "GB")) else 0
    Int disk_padGB = 1 + gatk_override_size + emergency_extra_diskGB

    if (!defined(runtime_params)) {
        Runtime standard_runtime = {
            "docker": gatk_docker,
            "jar_override": gatk_override,
            "max_retries": max_retries,
            "preemptible": preemptible,
            "cpu": 1,
            "machine_mem": 2024,
            "command_mem": 2024,
            "runtime_minutes": 60,
            "disk": 1 + disk_padGB,
            "boot_disk_size": 12  # needs to be > 10
        }
    }

    call ps.PileupSummaries as TumorPileupSummaries {
        input:
            interval_list = interval_list,
            scattered_interval_list = scattered_interval_list,
            bam = tumor_bam,
            bai = tumor_bai,
            ref_dict = ref_dict,
            variants = variants,
            variants_idx = variants_idx,
            vcf = vcf,
            vcf_idx = vcf_idx,
            minimum_population_allele_frequency = minimum_population_allele_frequency,
            maximum_population_allele_frequency = maximum_population_allele_frequency,
            getpileupsummaries_extra_args = getpileupsummaries_extra_args,
            output_base_name = tumor_sample_name,
            runtime_params = select_first([runtime_params, standard_runtime]),
            mem_get_pileup_summaries = mem_get_pileup_summaries,
            mem_gather_pileup_summaries = mem_gather_pileup_summaries,
            time_startup = time_startup,
            time_get_pileup_summaries = time_startup + ceil(time_get_pileup_summaries / scatter_count),
            time_gather_pileup_summaries = time_startup + time_gather_pileup_summaries
    }

    if (defined(normal_bam)) {
        call ps.PileupSummaries as NormalPileupSummaries {
            input:
                interval_list = interval_list,
                scattered_interval_list = scattered_interval_list,
                bam = tumor_bam,
                bai = tumor_bai,
                ref_dict = ref_dict,
                variants = variants,
                variants_idx = variants_idx,
                vcf = vcf,
                vcf_idx = vcf_idx,
                minimum_population_allele_frequency = minimum_population_allele_frequency,
                maximum_population_allele_frequency = maximum_population_allele_frequency,
                getpileupsummaries_extra_args = getpileupsummaries_extra_args,
                output_base_name = select_first([normal_sample_name, "normal"]),
                runtime_params = select_first([runtime_params, standard_runtime]),
                mem_get_pileup_summaries = mem_get_pileup_summaries,
                mem_gather_pileup_summaries = mem_gather_pileup_summaries,
                time_startup = time_startup,
                time_get_pileup_summaries = time_startup + ceil(time_get_pileup_summaries / scatter_count),
                time_gather_pileup_summaries = time_startup + time_gather_pileup_summaries
        }
    }

    call CalculateContamination {
        input:
            tumor_pileups = TumorPileupSummaries.pileup_summaries,
            normal_pileups = NormalPileupSummaries.pileup_summaries,
            runtime_params = select_first([runtime_params, standard_runtime]),
            memoryMB = mem_calculate_contamination,
            runtime_minutes = time_startup + time_calculate_contamination
    }

    output {
        File tumor_pileups = TumorPileupSummaries.pileup_summaries
        File? normal_pileups = NormalPileupSummaries.pileup_summaries
        File contamination_table = CalculateContamination.contamination_table
        File segmentation = CalculateContamination.segmentation
    }
}

task CalculateContamination {
    # This tool is set, by default, to operate over the intersection of your germline
    # biallelic resource with your specified intervals. The tool works by pulling down
    # alt and ref allele counts at locuses defined in the biallelic_germline_resource
    # and comparing the number of observed alt counts to the number of expected
    # alt counts at sites with low population allele frequencies. Using this comparison
    # of observed vs expected, the tool estimates contamination. It gains further
    # accuracy if a matched normal is specified.

    input {
        File tumor_pileups
        File? normal_pileups

        Runtime runtime_params
        Int? memoryMB
        Int? runtime_minutes
    }

    # Optional localization leads to cromwell error.
    # parameter_meta {
    #     tumor_pileups: {localization_optional: true}
    #     normal_pileups: {localization_optional: true}
    # }

    String tumor_sample_id = basename(tumor_pileups, ".pileup")
    String normal_sample_id = if defined(normal_pileups) then "." + basename(select_first([normal_pileups]), ".pileup") else ""
    String output_contamination = tumor_sample_id + normal_sample_id + ".contamination"
    String output_segments = tumor_sample_id + normal_sample_id + ".segments"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            CalculateContamination \
            --input '~{tumor_pileups}' \
            ~{"--matched-normal '" + normal_pileups + "'"} \
            --output '~{output_contamination}' \
            --tumor-segmentation '~{output_segments}'
    >>>

    output {
        File contamination_table = output_contamination
        File segmentation = output_segments
    }

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: select_first([memoryMB, runtime_params.machine_mem]) + " MB"
        runtime_minutes: select_first([runtime_minutes, runtime_params.runtime_minutes])
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}