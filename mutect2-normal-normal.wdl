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
import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/mutect2.wdl" as m2


workflow Mutect2NormalNormal {
	input {
		File? interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Array[File]+ bams
        Array[File]+ bais

        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
	}

	Array[Pair[File, File]] bam_pairs = cross(bams, bams)
	Array[Pair[File, File]] bai_pairs = cross(bais, bais)

	scatter (n in range(length(bam_pairs))) {
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

					run_cnn_scoring_model = false,
					run_funcotator = false
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

		gatk --java-options "-Xmx4g" \
			CountFalsePositives \
			-V ~{vcf} \
			-R ~{ref_fasta} \
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