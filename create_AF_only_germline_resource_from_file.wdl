version development

import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/create_AF_only_germline_resource.wdl" as singleAFonly

workflow CreateAFonlyVcf_from_File {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File vcfs_file
        File vcfs_idx_file

        Boolean compress_output = true

        # runtime
        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        Int emergency_extra_diskGB = 0
    }

    Array[File]+ vcfs = read_lines(vcfs_file)
    Array[File]+ vcfs_idx = read_lines(vcfs_idx_file)

    scatter (vcf_pair in zip(vcfs, vcfs_idx)) {
        call singleAFonly.CreateAFonlyVcf {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                vcf = vcf_pair.left,
                vcf_idx = vcf_pair.right,
                compress_output = compress_output,
                gatk_docker = gatk_docker,
                gatk_override = gatk_override,
                preemptible = preemptible,
                max_retries = max_retries,
                emergency_extra_diskGB = emergency_extra_diskGB
        }
    }

    call MergeVCFs as MergeAFonlyVCFs {
        input:
            input_vcfs = CreateAFonlyVcf.af_only_vcf,
            input_vcf_indices = CreateAFonlyVcf.af_only_vcf_idx,
            output_name = "af_only.filtered",
            compress_output = compress_output,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries,
    }

    call MergeVCFs as MergeBiallelicAFonlyVCFs {
        input:
            input_vcfs = CreateAFonlyVcf.biallelic_af_only_vcf,
            input_vcf_indices = CreateAFonlyVcf.biallelic_af_only_vcf_idx,
            output_name = "af_only.filtered.biallelic",
            compress_output = compress_output,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries,
    }

    call MergeVCFs as MergeMultiallelicAFonlyVCFs {
        input:
            input_vcfs = CreateAFonlyVcf.multiallelic_af_only_vcf,
            input_vcf_indices = CreateAFonlyVcf.multiallelic_af_only_vcf_idx,
            output_name = "af_only.filtered.multiallelic",
            compress_output = compress_output,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries,
    }

    output {
        File af_only_vcf = MergeAFonlyVCFs.merged_vcf
        File af_only_vcf_idx = MergeAFonlyVCFs.merged_vcf_idx
  		File biallelic_af_only_vcf = MergeBiallelicAFonlyVCFs.merged_vcf
  		File biallelic_af_only_vcf_idx = MergeBiallelicAFonlyVCFs.merged_vcf
        File multiallelic_af_only_vcf = MergeMultiallelicAFonlyVCFs.merged_vcf
        File multiallelic_af_only_vcf_idx = MergeMultiallelicAFonlyVCFs.merged_vcf
    }
}

task MergeVCFs {
    # Consider replacing MergeVcfs with GatherVcfsCloud once the latter is out of beta.

	input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        String output_name
        Boolean compress_output = false

        String gatk_docker
        File? gatk_override
        Int max_retries = 1
        Int preemptible = 1
        Int cpu = 1
        Int? disk
        Int boot_disk_size = 12
        Int machine_memoryMB = 512
        Int command_memoryMB = 512
    }

    # Optional localization leads to cromwell error.
    # parameter_meta {
    #     input_vcfs: {localization_optional: true}
    #     input_vcf_indices: {localization_optional: true}
    # }

    Int diskGB = ceil(1.2 * size(input_vcfs, "GB"))

    String output_vcf = output_name + ".vcf" + if compress_output then ".gz" else ""
    String output_vcf_idx = output_vcf + if compress_output then ".tbi" else ".idx"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
        gatk --java-options "-Xmx~{command_memoryMB}m" \
            MergeVcfs \
            ~{sep="' " prefix("-I '", input_vcfs)}' \
            -O '~{output_vcf}'
    >>>

    output {
    	File merged_vcf = output_vcf
        File merged_vcf_idx = output_vcf_idx
    }

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: boot_disk_size
        memory: machine_memoryMB + " MB"
        disks: "local-disk " + select_first([disk, diskGB]) + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
