version development

## For documentation, see mutect2_multi_sample.wdl

#import "mutect2_multi_sample.wdl" as msm2
import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/mutect2_multi_sample.wdl" as msm2

workflow Mutect2 {
    input {
        File? interval_list
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
        Boolean run_cnn_scoring_model = false  # probably leads to failure if true
        Boolean run_funcotator = true

        Boolean compress_output = true
        Boolean make_bamout = false
        Boolean keep_germline = false  # not currently supported
        Boolean genotype_germline_sites = false  # use with care!
        Boolean genotype_pon_sites = false  # use with care!
        Boolean native_pair_hmm_use_double_precision = true
        Boolean use_linked_de_bruijn_graph = true
        Boolean recover_all_dangling_branches = true
        Boolean funcotator_use_gnomad = true

        # expose extra arguments for import of this workflow
        String? split_intervals_extra_args
        String? m2_extra_args
        String? m2_filter_extra_args
        String? select_variants_extra_args
        String? select_low_conficence_variants_jexl_arg = "'(vc.getAttribute(\"GERMQ\") < 30) || (vc.getAttribute(\"DP\") < 4) || (vc.getAttribute(\"MBQ\").0 == 0) || (vc.getAttribute(\"MFRL\").0 == 0)'"
        String? realignment_extra_args
        String? funcotate_extra_args

        # resources
        File? panel_of_normals
        File? panel_of_normals_idx
        File? germline_resource
        File? germline_resource_tbi
        File? variants_for_contamination
        File? variants_for_contamination_idx
        File? bwa_mem_index_image
        File? funcotator_transcript_list
        File? funcotator_data_sources_tar_gz

        # runtime
        Int scatter_count = 42
        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 2
        Int emergency_extra_diskGB = 0

        # memory assignments in MB
        Int additional_per_sample_mem = 256  # this actually can depend on sample size (WES vs WGS)
        Int split_intervals_mem = 512  # 64
        Int get_sample_name_mem = 512  # 256
        Int variant_call_base_mem = 4096
        Int learn_read_orientation_model_base_mem = 6144
        Int get_pileup_summaries_mem = 2048  # needs at least 2G
        Int gather_pileup_summaries_mem = 512  # 64
        Int calculate_contamination_mem = 512
        Int filter_mutect_calls_mem = 4096
        Int select_variants_mem = 2048
        Int filter_alignment_artifacts_mem = 4096
        Int merge_vcfs_mem = 512
        Int merge_mutect_stats_mem = 512 # 64
        Int merge_bams_mem = 8192  # wants at least 6G
        Int cnn_scoring_mem = 4096
        Int funcotate_mem = 4096

        # Increasing cpus likely increases costs by the same factor.
        Int variant_call_cpu = 1  # good for PairHMM: 2
        Int filter_alignment_artifacts_cpu = 1  # good for PairHMM: 4
        Int cnn_scoring_cpu = 1
    }

    Runtime standard_runtime = {
        "gatk_docker": gatk_docker,
        "gatk_override": gatk_override,
        "max_retries": max_retries,
        "preemptible": preemptible,
        "cpu": 1,
        "machine_mem": 2024,
        "command_mem": 2024,
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

            compress_output = compress_output,
            make_bamout = make_bamout,

            keep_germline = keep_germline,
            genotype_germline_sites = genotype_germline_sites,
            genotype_pon_sites = genotype_pon_sites,
            native_pair_hmm_use_double_precision = native_pair_hmm_use_double_precision,
            use_linked_de_bruijn_graph = use_linked_de_bruijn_graph,
            recover_all_dangling_branches = recover_all_dangling_branches,
            funcotator_use_gnomad = funcotator_use_gnomad,

            split_intervals_extra_args = split_intervals_extra_args,
            m2_extra_args = m2_extra_args,
            m2_filter_extra_args = m2_filter_extra_args,
            select_variants_extra_args = select_variants_extra_args,
            select_low_conficence_variants_jexl_arg = select_low_conficence_variants_jexl_arg,
            realignment_extra_args = realignment_extra_args,
            funcotate_extra_args = funcotate_extra_args,

            panel_of_normals = panel_of_normals,
            panel_of_normals_idx = panel_of_normals_idx,
            germline_resource = germline_resource,
            germline_resource_tbi = germline_resource_tbi,
            variants_for_contamination = variants_for_contamination,
            variants_for_contamination_idx = variants_for_contamination_idx,
            bwa_mem_index_image = bwa_mem_index_image,
            funcotator_transcript_list = funcotator_transcript_list,
            funcotator_data_sources_tar_gz = funcotator_data_sources_tar_gz,

            scatter_count = scatter_count,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries,
            emergency_extra_diskGB = emergency_extra_diskGB,

            additional_per_sample_mem = additional_per_sample_mem,
            split_intervals_mem = split_intervals_mem,
            get_sample_name_mem = get_sample_name_mem,
            variant_call_base_mem = variant_call_base_mem,
            learn_read_orientation_model_base_mem = learn_read_orientation_model_base_mem,
            get_pileup_summaries_mem = get_pileup_summaries_mem,
            gather_pileup_summaries_mem = gather_pileup_summaries_mem,
            calculate_contamination_mem = calculate_contamination_mem,
            filter_mutect_calls_mem = filter_mutect_calls_mem,
            select_variants_mem = select_variants_mem,
            filter_alignment_artifacts_mem = filter_alignment_artifacts_mem,
            merge_vcfs_mem = merge_vcfs_mem,
            merge_mutect_stats_mem = merge_mutect_stats_mem,
            merge_bams_mem = merge_bams_mem,
            cnn_scoring_mem = cnn_scoring_mem,
            funcotate_mem = funcotate_mem,

            variant_call_cpu = variant_call_cpu,
            filter_alignment_artifacts_cpu = filter_alignment_artifacts_cpu,
            cnn_scoring_cpu = cnn_scoring_cpu
    }

    output {
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