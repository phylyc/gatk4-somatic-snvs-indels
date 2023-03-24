version development

## For documentation, see mutect2_multi_sample.wdl

#import "mutect2_multi_sample.wdl" as msm2
import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/mutect2_multi_sample.wdl" as msm2

workflow Mutect2 {
    input {
        String? individual_id
        File tumor_bam
        File tumor_bai
        File? normal_bam
        File? normal_bai

        # runtime
        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 2
    }

    GATKRuntime standard_runtime = {
        "gatk_docker": gatk_docker,
        "gatk_override": gatk_override,
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
            individual_id = select_first([individual_id, GetSampleName.sample_name]),
            tumor_bams = [tumor_bam],
            tumor_bais = [tumor_bai],
            normal_bams = if defined(normal_bam) then select_all([normal_bam]) else None,
            normal_bais = if defined(normal_bai) then select_all([normal_bai]) else None,

            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries
    }

    output {
        Array[File]? covered_intervals = MultiSampleMutect2.covered_intervals
        File evaluation_intervals = MultiSampleMutect2.evaluation_intervals
        File unfiltered_vcf = MultiSampleMutect2.unfiltered_vcf
        File unfiltered_vcf_idx = MultiSampleMutect2.unfiltered_vcf_idx
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