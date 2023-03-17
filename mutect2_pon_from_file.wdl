version development

## Create a Mutect2 panel of normals
##
## Description of inputs
## intervals: genomic intervals
## ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
## normal_bams_file: file with list of normal bam paths
##      This is necessary if the array of bam paths consists of more than 16384 characters
##      because Terra is stupid.
## scatter_count: number of parallel jobs when scattering over intervals
## pon_name: the resulting panel of normals is {pon_name}.vcf
## m2_extra_args: additional command line parameters for Mutect2. This should not
## involve --max-mnp-distance, which the wdl hard-codes to 0 because GenpmicsDBImport
## can't handle MNPs

# import "mutect2_multi_sample.wdl" as msm2
import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/mutect2_pon.wdl" as m2pon

workflow Mutect2_Panel_from_File {
    input {
        File? interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File normal_bams_file
        File normal_bais_file

        File germline_resource
        File germline_resource_idx

        Boolean compress_output = true
        String mutect2_extra_args = ""
        String pon_name

        Int min_contig_size = 1000000
        Int num_contigs = 24

        # runtime
        Int scatter_count = 10
        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 2
        Int max_retries = 2

        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int emergency_extra_diskGB = 10
    }

    call m2pon.Mutect2_Panel {
        input:
            interval_list = interval_list,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,

            normal_bams = read_lines(normal_bams_file),
            normal_bais = read_lines(normal_bais_file),

            germline_resource = germline_resource,
            germline_resource_idx = germline_resource_idx,

            compress_output = compress_output,
            mutect2_extra_args = mutect2_extra_args,
            pon_name = pon_name,

            min_contig_size = min_contig_size,
            num_contigs = num_contigs,

            scatter_count = scatter_count,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries,

            emergency_extra_diskGB = emergency_extra_diskGB
    }

    output {
        File pon = Mutect2_Panel.pon
        File pon_idx = Mutect2_Panel.pon_idx
        Array[File] normal_calls = Mutect2_Panel.normal_calls
        Array[File] normal_calls_idx = Mutect2_Panel.normal_calls_idx
    }
}