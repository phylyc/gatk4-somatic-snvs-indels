version development

## Create a Mutect2 panel of normals
##
## Description of inputs
## intervals: genomic intervals
## ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
## normal_bams: arrays of normal bams
## scatter_count: number of parallel jobs when scattering over intervals
## pon_name: the resulting panel of normals is {pon_name}.vcf
## m2_extra_args: additional command line parameters for Mutect2. This should not
## involve --max-mnp-distance, which the wdl hard-codes to 0 because GenpmicsDBImport
## can't handle MNPs

# import "mutect2_multi_sample.wdl" as msm2
import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/mutect2_multi_sample.wdl" as msm2

workflow Mutect2_Panel {
    input {
        File? interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Array[File]+ normal_bams
        Array[File]+ normal_bais

        File germline_resource
        File germline_resource_idx

        Boolean compress_output = true
        String m2_extra_args = ""
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

    Int gatk_override_size = if defined(gatk_override) then ceil(size(gatk_override, "GB")) else 0
    Int disk_padGB = 1 + gatk_override_size + emergency_extra_diskGB

    Runtime standard_runtime = {
        "gatk_docker": gatk_docker,
        "gatk_override": gatk_override,
        "max_retries": max_retries,
        "preemptible": preemptible,
        "cpu": 1,
        "machine_mem": 2024,
        "command_mem": 2024,
        "disk": 1 + disk_padGB,
        "boot_disk_size": 12  # needs to be > 10
    }

    scatter (normal_bam in zip(normal_bams, normal_bais)) {
        call msm2.GetSampleName {
            input:
                bam = normal_bam.left,
                runtime_params = standard_runtime,
        }

        call msm2.MultiSampleMutect2 {
            input:
                interval_list = interval_list,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                individual_id = GetSampleName.sample_name,
                tumor_bams = [normal_bam.left],
                tumor_bais = [normal_bam.right],
                run_contamination_model = false,
                run_orientation_bias_mixture_model = false,
                run_variant_filter = false,
                run_realignment_filter = false,
                run_funcotator = false,
                compress_output = compress_output,
                scatter_count = scatter_count,
                m2_extra_args = m2_extra_args + " --max-mnp-distance 0",
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                preemptible = preemptible,
                max_retries = max_retries
        }
    }

    String split_intervals_extra_args = (
        "--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION"
        + " --min-contig-size " + min_contig_size
    )
    call msm2.SplitIntervals {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            scatter_count = num_contigs,
            split_intervals_extra_args = split_intervals_extra_args,
            runtime_params = standard_runtime
    }

    scatter (scattered_intervals in SplitIntervals.interval_files) {
        call CreatePanel {
            input:
                input_vcfs = MultiSampleMutect2.merged_vcf,
                input_vcf_indices = MultiSampleMutect2.merged_vcf_idx,
                interval_list = scattered_intervals,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                compress_output = compress_output,
                gnomad = germline_resource,
                gnomad_idx = germline_resource_idx,
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
        Array[File] normal_calls = MultiSampleMutect2.merged_vcf
        Array[File] normal_calls_idx = MultiSampleMutect2.merged_vcf
    }
}

task CreatePanel {
    input {
        File interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Array[File]+ input_vcfs
        Array[File]+ input_vcf_indices
        String output_vcf_name

        Boolean compress_output = false
        File gnomad
        File gnomad_idx

        # UMCCR found that 5 is optimizing the F2 score, but not by much compared to 2:
        # https://umccr.org/blog/panel-of-normals/
        Int min_sample_count = 2
        String? create_pon_extra_args

        # runtime
        Runtime runtime_params
        Int memoryMB = 8186
    }

    # GenomicsDB requires that the reference be a local file.
    parameter_meta{
        gnomad: {localization_optional: true}
        gnomad_idx: {localization_optional: true}
    }

    Int vcf_size = 2 * ceil(size(input_vcfs, "GB"))
    Int disk_size = runtime_params.disk + vcf_size

    String output_file = output_vcf_name + ".vcf.gz"
    String output_file_idx = output_file + ".tbi"

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            GenomicsDBImport \
            --genomicsdb-workspace-path pon_db \
            -R ~{ref_fasta} \
            -V ~{sep=' -V ' input_vcfs} \
            -L ~{interval_list}

        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            CreateSomaticPanelOfNormals \
            -R ~{ref_fasta} \
            --germline-resource ~{gnomad} \
            -V gendb://pon_db \
            -O ~{output_file} \
            --min-sample-count ~{min_sample_count} \
            ~{create_pon_extra_args}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: memoryMB + " MB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File output_vcf = output_file
        File output_vcf_index = output_file_idx
    }
}