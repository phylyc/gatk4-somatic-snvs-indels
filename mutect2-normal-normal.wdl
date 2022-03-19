version development

## Copyright Broad Institute, 2017
## Philipp Hähnel, 2022
##
## This WDL workflow calls pairs of replicates as tumor-normal pairs, 
## counts the number of variants (i.e. false positives) and reports false 
## positive rates.
##
## Main requirements/expectations :
## - One analysis-ready BAM file (and its index) for each replicate
##
## Outputs :
## - False Positive VCF files and its index with summary
##
## Cromwell version support
## - Successfully tested on v71
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## pages at https://hub.docker.com/r/broadinstitute/* for detailed licensing information
## pertaining to the included programs.

# import "mutect2.wdl" as m2
import "https://github.com/phylyc/gatk4-somatic-snvs-indels/mutect2.wdl" as m2

workflow Mutect2NormalNormal {
	input {
		File? interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Array[File]+ bams
        Array[File]+ bais

        # workflow options
        Boolean run_contaminanation_model = true
        Boolean run_orientation_bias_mixture_model_filter = true
        Boolean run_variant_filter = true
        Boolean run_realignment_filter = true
        Boolean keep_germline = false  # not currently supported
        Boolean compress_output = true
        Boolean make_bamout = false

        # expose extra arguments for import of this workflow
        String? split_intervals_extra_args
        String? m2_extra_args
        String? m2_filter_extra_args
        String? select_variants_extra_args
        String? realignment_extra_args

        # resources
        File? panel_of_normals
        File? panel_of_normals_idx
        File? germline_resource
        File? germline_resource_tbi
		File? variants_for_contamination
		File? variants_for_contamination_idx
        File? bwa_mem_index_image

        # runtime
        Int scatter_count = 42
        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 2
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
        Int filter_mutect_calls_mem = 512
        Int select_variants_mem = 512
        Int filter_alignment_artifacts_mem = 4096
        Int merge_vcfs_mem = 512
        Int merge_mutect_stats_mem = 512 # 64
        Int merge_bams_mem = 8192  # wants at least 6G

        # Increasing cpus likely increases costs by the same factor.
        Int variant_call_cpu = 1  # good for PairHMM: 2
        Int filter_alignment_artifacts_cpu = 1  # good for PairHMM: 4
	}

	Array[Pair[File, File]] bam_pairs = cross(bams, bams)
	Array[Pair[File, File]] bai_pairs = cross(bais, bais)

	scatter(n in range(length(bam_pairs))) {
	    File tumor_bam = bam_pairs[n].left
	    File normal_bam = bam_pairs[n].right
	    File tumor_bai = bai_pairs[n].left
		File normal_bai = bai_pairs[n].right

        if (tumor_bam != normal_bam) {
            call m2.Mutect2 {
                input:
                    interval_list = interval_list,
					ref_fasta = ref_fasta,
					ref_fasta_index = ref_fasta_index,
					ref_dict = ref_dict,

					tumor_bam = tumor_bam,
					tumor_bai = tumor_bai,
					normal_bam = normal_bam,
					normal_bai = normal_bai,

					run_contaminanation_model = run_contaminanation_model,
					run_orientation_bias_mixture_model_filter = run_orientation_bias_mixture_model_filter,
					run_variant_filter = run_variant_filter,
					run_realignment_filter = run_realignment_filter,
					run_funcotator = false,
					keep_germline = keep_germline,
					compress_output = compress_output,
					make_bamout = make_bamout,

					split_intervals_extra_args = split_intervals_extra_args,
					m2_extra_args = m2_extra_args,
					m2_filter_extra_args = m2_filter_extra_args,
					select_variants_extra_args = select_variants_extra_args,
					realignment_extra_args = realignment_extra_args,

					panel_of_normals = panel_of_normals,
					panel_of_normals_idx = panel_of_normals_idx,
					germline_resource = germline_resource,
					germline_resource_tbi = germline_resource_tbi,
					variants_for_contamination = variants_for_contamination,
					variants_for_contamination_idx = variants_for_contamination_idx,
					bwa_mem_index_image = bwa_mem_index_image,

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

					variant_call_cpu = variant_call_cpu,
					filter_alignment_artifacts_cpu = filter_alignment_artifacts_cpu
            }

            call CountFalsePositives {
                input:
                    interval_list = interval_list,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_dict = ref_dict,
                    vcf = Mutect2.merged_vcf,
                    vcf_idx = Mutect2.merged_vcf_idx,
                    gatk_override = gatk_override,
                    gatk_docker = gatk_docker
            }
		}
	}

	call GatherTables {
		input:
			tables = select_all(CountFalsePositives.false_positive_counts)
	}

	output {
		File summary = GatherTables.summary
		Array[File] false_positives_vcfs = select_all(Mutect2.merged_vcf)
		Array[File] false_positives_vcf_indices = select_all(Mutect2.merged_vcf_idx)
	}
}

task GatherTables {
	input {
		# we assume that each table consists of two lines: one header line and one record
		Array[File] tables
	}

	String dollar = "$"

	command <<<
	    # extract the header from one of the files
		head -n 1 ~{tables[0]} > summary.txt

		# then append the record from each table
		for table in ~{sep=" " tables}; do
			tail -n +2 ~{dollar}table >> summary.txt
		done
	>>>

	runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
        memory: "1 GB"
        disks: "local-disk " + 1 + " HDD"
    }

	output {
		File summary = "summary.txt"
	}
}

task CountFalsePositives {
	input {
		File? interval_list
		File ref_fasta
		File ref_fasta_index
		File ref_dict
		File vcf
		File vcf_idx

		File? gatk_override
		String gatk_docker
	}

	command <<<
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

	    gatk --java-options "-Xmx4g" CountFalsePositives \
		    -V ${vcf} \
		    -R ${ref_fasta} \
		    ~{"-L " + interval_list} \
		    -O false-positives.txt \
	>>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: "5 GB"
        disks: "local-disk " + 5 + " HDD"
    }

	output {
		File false_positive_counts = "false-positives.txt"
	}
}