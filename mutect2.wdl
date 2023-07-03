version development

## For documentation, see mutect2_multi_sample.wdl

#import "mutect2_multi_sample.wdl" as msm2
import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/mutect2_multi_sample.wdl" as msm2


workflow Mutect2 {
    input {
        File? interval_list
        Array[File]? interval_lists
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        String? individual_id
        File tumor_bam
        File tumor_bai
        File? normal_bam
        File? normal_bai

        # workflow options
        Boolean run_contamination_model = true
        Boolean run_orientation_bias_mixture_model = true
        Boolean run_variant_filter = true
        Boolean run_realignment_filter = true
        Boolean run_realignment_filter_only_on_high_confidence_variants = false
        Boolean run_cnn_scoring_model = false  # likely leads to failure if true; better performance with make_bamout = true
        Boolean run_funcotator = true

        # runtime
        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 2
    }

    Runtime standard_runtime = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "max_retries": max_retries,
        "preemptible": preemptible,
        "cpu": 1,
        "machine_mem": 2024,
        "command_mem": 2024,
        "runtime_minutes": 60,
        "disk": 1,
        "boot_disk_size": 12  # needs to be > 10
    }

    call msm2.GetSampleName {
        input:
            bam = tumor_bam,
            runtime_params = standard_runtime
    }

    call msm2.MultiSampleMutect2 {
        input:
            interval_list = interval_list,
            interval_lists = interval_lists,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,

            individual_id = select_first([individual_id, GetSampleName.sample_name]),
            tumor_bams = [tumor_bam],
            tumor_bais = [tumor_bai],
            normal_bams = if defined(normal_bam) then select_all([normal_bam]) else None,
            normal_bais = if defined(normal_bai) then select_all([normal_bai]) else None,

            run_contamination_model = run_contamination_model,
            run_orientation_bias_mixture_model = run_orientation_bias_mixture_model,
            run_variant_filter = run_variant_filter,
            run_realignment_filter = run_realignment_filter,
            run_realignment_filter_only_on_high_confidence_variants = run_realignment_filter_only_on_high_confidence_variants,
            run_cnn_scoring_model = run_cnn_scoring_model,
            run_funcotator = run_funcotator,

            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries
    }

    output {
        File unfiltered_vcf = MultiSampleMutect2.unfiltered_vcf
        File unfiltered_vcf_idx = MultiSampleMutect2.unfiltered_vcf_idx
        File? filtered_vcf = MultiSampleMutect2.filtered_vcf
        File? filtered_vcf_idx = MultiSampleMutect2.filtered_vcf_idx
        File merged_vcf = MultiSampleMutect2.merged_vcf
        File merged_vcf_idx = MultiSampleMutect2.merged_vcf_idx
        File mutect_stats = MultiSampleMutect2.mutect_stats
        File? bamout = MultiSampleMutect2.bamout
        File? bamout_index = MultiSampleMutect2.bamout_index
        File? filtering_stats = MultiSampleMutect2.filtering_stats
        File? read_orientation_model_params = MultiSampleMutect2.read_orientation_model_params
        Array[File]? contamination_table = MultiSampleMutect2.contamination_table
        Array[File]? tumor_segmentation = MultiSampleMutect2.tumor_segmentation
        Array[File?]? scored_file = MultiSampleMutect2.scored_file
        Array[File?]? scored_file_idx = MultiSampleMutect2.scored_file_idx
        Array[File?]? funcotated_file = MultiSampleMutect2.funcotated_file
        Array[File?]? funcotated_file_index = MultiSampleMutect2.funcotated_file_index
    }
}