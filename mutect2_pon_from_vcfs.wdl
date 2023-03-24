version development

## Create a Mutect2 panel of normals
##
## Description of inputs
## intervals: genomic intervals
## ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
## normal_vcfs_file: file with list of Mutect2 vcf outputs from normal bams
##      This is necessary if the array of bam paths consists of more than 16384 characters
##      because Terra is stupid.
## scatter_count: number of parallel jobs when scattering over intervals
## pon_name: the resulting panel of normals is {pon_name}.vcf
## m2_extra_args: additional command line parameters for Mutect2. This should not
## involve --max-mnp-distance, which the wdl hard-codes to 0 because GenpmicsDBImport
## can't handle MNPs

# import "mutect2_multi_sample.wdl" as msm2
import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/mutect2_multi_sample.wdl" as msm2
import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/mutect2_pon.wdl" as m2pon


workflow Mutect2_Panel_from_VCFs {
    input {
        File? interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File normal_vcfs_file
        File normal_vcf_indices_file

        Boolean compress_output = true
        String pon_name

        Int min_contig_size = 1000000

        # runtime
        Int scatter_count = 24
        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 2
        Int max_retries = 2

        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int emergency_extra_diskGB = 10
    }

    Int gatk_override_size = if defined(gatk_override) then ceil(size(gatk_override, "GB")) else 0
    Int disk_padGB = 1 + gatk_override_size + emergency_extra_diskGB

    Runtime standard_runtime = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "max_retries": max_retries,
        "preemptible": preemptible,
        "cpu": 1,
        "machine_mem": 4096,
        "command_mem": 2048,
        "runtime_minutes": 60,
        "disk": 1 + disk_padGB,
        "boot_disk_size": 12  # needs to be > 10
    }

    String split_intervals_extra_args = (
        "--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION"
        + " --min-contig-size " + min_contig_size
    )
    call msm2.SplitIntervals {
        input:
            interval_list = interval_list,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            scatter_count = scatter_count,
            split_intervals_extra_args = split_intervals_extra_args,
            runtime_params = standard_runtime
    }

    scatter (scattered_intervals in SplitIntervals.interval_files) {
        call m2pon.CreatePanel {
            input:
                input_vcfs = read_lines(normal_vcfs_file),
                input_vcf_indices = read_lines(normal_vcf_indices_file),
                interval_list = scattered_intervals,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                compress_output = compress_output,
                output_vcf_name = pon_name,
                runtime_params = standard_runtime
        }
    }

    call msm2.MergeVCFs {
        input:
            input_vcfs = CreatePanel.output_vcf,
            input_vcf_indices = CreatePanel.output_vcf_index,
            output_name = pon_name,
            compress_output = compress_output,
            runtime_params = standard_runtime
    }

    output {
        File pon = MergeVCFs.merged_vcf
        File pon_idx = MergeVCFs.merged_vcf_idx
    }
}