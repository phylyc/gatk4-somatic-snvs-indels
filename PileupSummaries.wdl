version development

#import "util.wdl" as util
import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/util.wdl" as util


workflow PileupSummaries {
	input {
        File? interval_list
        Array[File]? scattered_interval_list
        File bam
        File bai
        File ref_dict

        File? variants
        File? variants_idx
        File? vcf
        File? vcf_idx
        String? getpileupsummaries_extra_args

        String output_name

        Runtime? runtime_params

        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        Int emergency_extra_diskGB = 0

        Int mem_get_pileup_summaries = 4096  # needs at least 2G
        Int mem_gather_pileup_summaries = 512  # 64
        Int time_startup = 10
        Int time_get_pileup_summaries = 90  # 1.5 h
        Int time_gather_pileup_summaries = 5
	}

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

    if (defined(vcf) && !defined(variants)) {
        call ToPileupVCF {
            input:
                vcf = select_first([vcf]),
                vcf_idx = select_first([vcf_idx]),
                runtime_params = select_first([runtime_params, standard_runtime]),
        }
    }

    if (defined(scattered_interval_list)) {
        scatter (scattered_intervals in select_all(select_first([scattered_interval_list]))) {
            call GetPileupSummaries as ScatteredGetPileupSummaries {
                input:
                    input_bam = bam,
                    input_bai = bai,
                    interval_list = scattered_intervals,
                    variants = select_first([variants, ToPileupVCF.variants]),
                    variants_idx = select_first([variants_idx, ToPileupVCF.variants_idx]),
                    runtime_params = select_first([runtime_params, standard_runtime]),
                    memoryMB = mem_get_pileup_summaries,
                    runtime_minutes = time_startup + time_get_pileup_summaries
            }
        }

        call GatherPileupSummaries {
            input:
                input_tables = flatten(ScatteredGetPileupSummaries.pileup_summaries),
                ref_dict = ref_dict,
                output_name = output_name,
                runtime_params = select_first([runtime_params, standard_runtime]),
                memoryMB = mem_gather_pileup_summaries,
                runtime_minutes = time_startup + time_gather_pileup_summaries
        }
    }
    # else
    if (!defined(scattered_interval_list)) {
        call GetPileupSummaries {
            input:
                input_bam = bam,
                input_bai = bai,
                interval_list = interval_list,
                variants = select_first([variants, ToPileupVCF.variants]),
                variants_idx = select_first([variants_idx, ToPileupVCF.variants_idx]),
                runtime_params = select_first([runtime_params, standard_runtime]),
                memoryMB = mem_get_pileup_summaries,
                runtime_minutes = time_startup + time_get_pileup_summaries
        }
    }

    output {
        File pileup_summaries = select_first([GatherPileupSummaries.merged_pileup_summaries, select_first(GetPileupSummaries.pileup_summaries)])
    }
}

task ToPileupVCF {
    # Input: a (multi-sample) VCF, e.g. from Mutect2
    # create a VCF with AF=1.0 for all variants in the input VCF, dropping all
    # samples, resulting in a gnomad-style VCF with only AF=1.0 in the INFO field.

    input {
        File vcf
        File vcf_idx

        Float AF = 1.0

        Runtime runtime_params
        String? docker = "dceoy/bcftools"
    }

    Int diskGB = ceil(2 * size(vcf, "GB"))

    String output_base_name = basename(basename(basename(vcf, ".gz"), ".bgz"), ".vcf")
    String tmp_vcf = output_base_name + ".tmp.vcf"
    String uncompressed_vcf = output_base_name + ".af_only.vcf"
    String af_only_vcf = uncompressed_vcf + ".gz"
    String af_only_vcf_idx = af_only_vcf + ".tbi"

    String dollar = "$"

    command <<<
        set -e
        bcftools view -G '~{vcf}' \
            | bcftools annotate -x FORMAT,^INFO/AF \
            > '~{tmp_vcf}'
        grep "^#" '~{tmp_vcf}' > '~{uncompressed_vcf}'
        grep -v "^#" '~{tmp_vcf}' \
            | awk 'BEGIN{OFS="\t"}{~{dollar}8 = "AF=~{AF}"; print}' \
            >> '~{uncompressed_vcf}'
        bgzip -c '~{uncompressed_vcf}' > '~{af_only_vcf}'
        bcftools index -t -o ~{af_only_vcf_idx} ~{af_only_vcf}
        rm -f '~{tmp_vcf}' '~{uncompressed_vcf}'
    >>>

    output {
        File variants = af_only_vcf
        File variants_idx = af_only_vcf_idx
    }

    runtime {
        docker: select_first([docker, runtime_params.docker])
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        runtime_minutes: runtime_params.runtime_minutes
        disks: "local-disk " + diskGB + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task GetPileupSummaries {
    # If the variants for contamination and the intervals for this scatter don't
    # intersect, GetPileupSummaries throws an error. However, there is nothing wrong
    # with an empty intersection for our purposes; it simply doesn't contribute to the
    # merged pileup summaries that we create downstream. We implement this with an array
    # outputs. If the tool errors, no table is created and the glob yields an empty array.

	input {
        File? interval_list
        File input_bam
        File input_bai
        File? variants
        File? variants_idx
        String? getpileupsummaries_extra_args

        Runtime runtime_params
        Int? memoryMB
        Int? runtime_minutes
        Int? disk_spaceGB
	}

    parameter_meta {
        interval_list: {localization_optional: true}
        input_bam: {localization_optional: true}
        input_bai: {localization_optional: true}
        variants: {localization_optional: true}
        variants_idx: {localization_optional: true}
    }

    String sample_id = basename(input_bam, ".bam")
    String output_name = sample_id + ".pileup"

    command <<<
        set +e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}

        if ~{!defined(variants)} ; then
            echo "ERROR: variants must be supplied."
            false
        fi

        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            GetPileupSummaries \
            --input '~{input_bam}' \
            ~{"--intervals '" +  interval_list + "'"} \
            --intervals '~{variants}' \
            --interval-set-rule INTERSECTION \
            --variant '~{variants}' \
            --output '~{output_name}' \
            ~{getpileupsummaries_extra_args}

        # It only fails due to empty intersection between variants and intervals, which is ok.
        exit 0
    >>>

    # need to use glob in case GetPileupSummaries fails and does not produce output files
    output {
        Array[File] pileup_summaries = glob(output_name)
    }

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: select_first([memoryMB, runtime_params.machine_mem]) + " MB"
        runtime_minutes: select_first([runtime_minutes, runtime_params.runtime_minutes])
        disks: "local-disk " + select_first([disk_spaceGB, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task GatherPileupSummaries {
    input {
        Array[File] input_tables
        File ref_dict
        String output_name

        Runtime runtime_params
        Int? memoryMB
        Int? runtime_minutes
    }

    # Optional localization leads to cromwell error.
    # parameter_meta {
    #     input_tables: {localization_optional: true}
    #     ref_dict: {localization_optional: true}
    # }

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            GatherPileupSummaries \
            --sequence-dictionary '~{ref_dict}' \
            ~{sep="' " prefix("-I '", input_tables)}' \
            -O '~{output_name}'
    >>>

    output {
        File merged_pileup_summaries = output_name
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
