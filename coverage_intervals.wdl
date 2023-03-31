version development


import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/util.wdl" as util


workflow CoverageIntervals {
    input {
        File? interval_list
        Array[File]? interval_lists
        File? intervals_for_coverage
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        String collection_name
        Array[String]+ sample_names
        Array[File]+ input_bams
        Array[File]+ input_bais

        # arguments
        Boolean paired_end = false
        Int min_read_depth_threshold = 1
        Int intervals_bin_length = 0
        Int intervals_padding = 0

        # runtime
        String bedtools_docker = "staphb/bedtools"
        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        Int emergency_extra_diskGB = 0

        # memory assignments in MB
        Int mem_get_genome_coverage = 8192
        Int mem_get_evaluation_intervals = 8192

        # runtime assignments in minutes (for HPC cluster)
        Int time_startup = 10
        Int time_apply_read_filters = 120
        Int time_get_genome_coverage = 60
        Int time_get_evaluation_intervals = 60
    }

    Runtime standard_runtime = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "max_retries": max_retries,
        "preemptible": preemptible,
        "cpu": 1,
        "machine_mem": 4096,
        "command_mem": 4096,
        "runtime_minutes": 60,
        "disk": 3 + emergency_extra_diskGB,
        "boot_disk_size": 12  # needs to be > 10
    }

    scatter (sample in zip(sample_names, zip(input_bams, input_bais))) {
        String sample_name = sample.left
        File bam = sample.right.left
        File bai = sample.right.right

        call ApplyReadFilters {
            input:
                input_bam = bam,
                input_bai = bai,
                interval_list = intervals_for_coverage,
                runtime_params = standard_runtime,
                runtime_minutes = time_startup + time_apply_read_filters
        }

        call GetGenomeCoverage {
            input:
#                ref_fasta = ref_fasta,
#                ref_fasta_index = ref_fasta_index,
#                ref_dict = ref_dict,
                sample_name = sample_name,
                input_bam = ApplyReadFilters.bam,
                input_bai = ApplyReadFilters.bai,
                min_read_depth_threshold = min_read_depth_threshold,
                paired_end = paired_end,
                runtime_params = standard_runtime,
                bedtools_docker = bedtools_docker,
                memoryMB = mem_get_genome_coverage,
                runtime_minutes = time_startup + time_get_genome_coverage
        }

        call BedToIntervalList as GetGenomeCoverageIntervals {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                bed_file = GetGenomeCoverage.covered_regions,
                runtime_params = standard_runtime,
        }
    }

    call GetEvaluationIntervals {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            collection_name = collection_name,
            interval_list = interval_list,
            interval_lists = interval_lists,
            covered_intervals = GetGenomeCoverageIntervals.interval_list,
            bin_length = intervals_bin_length,
            padding = intervals_padding,
            runtime_params = standard_runtime,
            memoryMB = mem_get_evaluation_intervals,
            runtime_minutes = time_startup + time_get_evaluation_intervals
    }

    output {
        Array[File] covered_intervals = GetGenomeCoverageIntervals.interval_list
        File evaluation_intervals = GetEvaluationIntervals.evaluation_interval_list
    }
}

task ApplyReadFilters {
    input {
        File? interval_list
        File input_bam
        File input_bai

        Array[String]? read_filters
        String? extra_args

        Runtime runtime_params
        Int? memoryMB
        Int? runtime_minutes
    }

    Int diskGB = ceil(1.3 * size(input_bam, "GB")) + runtime_params.disk
    String filtered_bam = basename(input_bam, ".bam") + ".filtered.bam"
    String filtered_bai = basename(input_bai, ".bai") + ".filtered.bai"

    parameter_meta {
        input_bam: {localization_optional: true}
        input_bai: {localization_optional: true}
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}

        echo ""
        echo "    Apply Read Filters that are automatically applied to the data by the Engine before processing by Mutect2. "
        echo ""

        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            PrintReads \
            -I '~{input_bam}' \
            -O '~{filtered_bam}' \
            ~{"-L '" + interval_list + "'"} \
            --read-filter NonChimericOriginalAlignmentReadFilter \
            --read-filter NotSecondaryAlignmentReadFilter \
            --read-filter GoodCigarReadFilter \
            --read-filter NonZeroReferenceLengthAlignmentReadFilter \
            --read-filter PassesVendorQualityCheckReadFilter \
            --read-filter MappedReadFilter \
            --read-filter MappingQualityAvailableReadFilter \
            --read-filter NotDuplicateReadFilter \
            --read-filter MappingQualityReadFilter \
            --read-filter MappingQualityNotZeroReadFilter \
            --read-filter WellformedReadFilter \
            ~{true="--read-filter " false="" defined(read_filters)}~{default="" sep=" --read-filter " read_filters} \
            ~{extra_args}
    >>>

    output {
        File bam = filtered_bam
        File bai = filtered_bai
    }

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: select_first([memoryMB, runtime_params.machine_mem]) + " MB"
        runtime_minutes: select_first([runtime_minutes, runtime_params.runtime_minutes])
        disks: "local-disk " + diskGB + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task GetGenomeCoverage {
    input {
        String? sample_name
        File input_bam
        File input_bai

        Int min_read_depth_threshold = 1
        Boolean paired_end = false

        Runtime runtime_params
        String? bedtools_docker = "staphb/bedtools"
        Int? memoryMB
        Int? runtime_minutes
    }

    Int diskGB = ceil(1.3 * size(input_bam, "GB")) + runtime_params.disk
    Int max = if paired_end then ceil(min_read_depth_threshold / 2) else min_read_depth_threshold

    String filtered_bam = basename(input_bam, ".bam") + ".filtered.bam"
    String covered_regions_bed = sample_name + ".coveraged_regions.bed"
    String interval_list_file = basename(covered_regions_bed, ".bed") + ".interval_list"

    parameter_meta {
    }

    String dollar = "$"

    command <<<
        set -e

        echo "..."
        echo "    Collecting BEDGRAPH summaries of feature coverage"
        echo "    and removing regions below ~{min_read_depth_threshold} read depth:"
        echo "..."

        bedtools genomecov \
            -ibam '~{input_bam}' \
            -max ~{max} \
            -bg \
            ~{true="-pc " false="" paired_end} \
        | awk '~{dollar}4>=~{max}' \
        | bedtools merge -c 4 -o min -d 1 -i stdin \
        > '~{covered_regions_bed}'
    >>>

    output {
        File covered_regions = covered_regions_bed
    }

    runtime {
        docker: select_first([bedtools_docker, runtime_params.docker])
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: select_first([memoryMB, runtime_params.machine_mem]) + " MB"
        runtime_minutes: select_first([runtime_minutes, runtime_params.runtime_minutes])
        disks: "local-disk " + diskGB + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task BedToIntervalList {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File bed_file

        Runtime runtime_params
        Int? memoryMB
        Int? runtime_minutes
    }

    String interval_list_file = basename(bed_file, ".bed") + ".interval_list"

    parameter_meta {
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}

        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            BedToIntervalList \
            -I '~{bed_file}' \
            -O '~{interval_list_file}' \
            -SD '~{ref_dict}'
    >>>

    output {
        File interval_list = interval_list_file
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

#task GetGenomeCoverage {
#    input {
#        File ref_fasta
#        File ref_fasta_index
#        File ref_dict
#
#        String? sample_name
#        File input_bam
#        File input_bai
#
#        Int min_read_depth_threshold = 1
#
#        Runtime runtime_params
#        Int? memoryMB
#        Int? runtime_minutes
#    }
#
#    Int diskGB = ceil(1.2 * size(input_bam, "GB")) + runtime_params.disk
#
#    String filtered_bam = basename(input_bam, ".bam") + ".filtered.bam"
#    String covered_regions_bed = sample_name + ".coveraged_regions.bed"
#    String interval_list_file = basename(covered_regions_bed, ".bed") + ".interval_list"
#
#    parameter_meta {
##        ref_fasta: {localization_optional: true}
##        ref_fasta_index: {localization_optional: true}
##        ref_dict: {localization_optional: true}
#        input_bam: {localization_optional: true}
#        input_bai: {localization_optional: true}
#    }
#
#    String dollar = "$"
#
#    command <<<
#        set -e
#        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
#
#        echo "..."
#        echo "    Apply Read Filters that are automatically applied to the data by the Engine before processing by Mutect2:"
#        echo "..."
#
#        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
#            PrintReads \
#            -I '~{input_bam}' \
#            -O '~{filtered_bam}' \
#            --read-filter NonChimericOriginalAlignmentReadFilter \
#            --read-filter NotSecondaryAlignmentReadFilter \
#            --read-filter GoodCigarReadFilter \
#            --read-filter NonZeroReferenceLengthAlignmentReadFilter \
#            --read-filter PassesVendorQualityCheckReadFilter \
#            --read-filter MappedReadFilter \
#            --read-filter MappingQualityAvailableReadFilter \
#            --read-filter NotDuplicateReadFilter \
#            --read-filter MappingQualityReadFilter \
#            --read-filter MappingQualityNotZeroReadFilter \
#            --read-filter WellformedReadFilter
#
#        echo "..."
#        echo "    Collecting BEDGRAPH summaries of feature coverage"
#        echo "    and removing regions below ~{min_read_depth_threshold} read depth:"
#        echo "..."
#
#        bedtools genomecov \
#            -ibam '~{filtered_bam}' \
#            -max ~{min_read_depth_threshold} \
#            -bg \
#            -pc \
#        | awk '~{dollar}4>=~{min_read_depth_threshold}' \
#        | bedtools merge -c 4 -o min -d 1 -i stdin \
#        > '~{covered_regions_bed}'
#
#        echo "..."
#        echo "    Create interval_list from bed file:"
#        echo "..."
#
#        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
#            BedToIntervalList \
#            -I '~{covered_regions_bed}' \
#            -O '~{interval_list_file}' \
#            -SD '~{ref_dict}'
#    >>>
#
#    output {
##        File covered_regions = covered_regions_bed
#        File interval_list = interval_list_file
#    }
#
#    runtime {
#        docker: runtime_params.gatk_docker
#        bootDiskSizeGb: runtime_params.boot_disk_size
#        memory: select_first([memoryMB, runtime_params.machine_mem]) + " MB"
#        runtime_minutes: select_first([runtime_minutes, runtime_params.runtime_minutes])
#        disks: "local-disk " + diskGB + " HDD"
#        preemptible: runtime_params.preemptible
#        maxRetries: runtime_params.max_retries
#        cpu: runtime_params.cpu
#    }
#}

task GetEvaluationIntervals {
    input {
        File? interval_list
        Array[File]? interval_lists
        Array[File]+ covered_intervals
        File? interval_blacklist
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        String collection_name

        Int bin_length = 0
        Int padding = 0
        String? split_intervals_extra_args

        Runtime runtime_params
        Int? memoryMB
        Int? runtime_minutes
    }

    Int diskGB = (
        if defined(interval_list) then ceil(size(interval_list, "GB")) else 0
        + if defined(interval_lists) then ceil(2 * size(interval_lists, "GB")) else 0
        + ceil(2 * size(covered_intervals, "GB"))
        + runtime_params.disk
    )

    String merged_input_intervals = collection_name + ".merged_input.interval_list"
    String merged_covered_regions = collection_name + ".merged_covered_regions.interval_list"
    String intersection_intervals = collection_name + ".merged_input.CAP.merged_covered_regions.interval_list"
    String preprocessed_evaluation_intervals = collection_name + ".preprocessed_evaluation_regions.interval_list"

    parameter_meta {
#        interval_list: {localization_optional: true}
#        interval_lists: {localization_optional: true}
#        covered_intervals: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}

        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            IntervalListTools \
            --ACTION UNION \
            ~{"-I '" + interval_list + "'"} \
            ~{true="-I '" false="" defined(interval_lists)}~{default="" sep="' -I '" interval_lists}~{true="'" false="" defined(interval_lists)} \
            -O '~{merged_input_intervals}'

        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            IntervalListTools \
            --ACTION UNION \
            ~{sep="' " prefix("-I '", covered_intervals)}' \
            -O '~{merged_covered_regions}'

        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            IntervalListTools \
            --ACTION INTERSECT \
            -I '~{merged_input_intervals}' \
            -I '~{merged_covered_regions}' \
            -O '~{intersection_intervals}'

        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            PreprocessIntervals \
            -R '~{ref_fasta}' \
            -L '~{intersection_intervals}' \
            ~{"-XL '" + interval_blacklist + "'"} \
            --bin-length ~{bin_length} \
            --padding ~{padding} \
            --interval-merging-rule OVERLAPPING_ONLY \
            -O '~{preprocessed_evaluation_intervals}'
    >>>

    output {
        File evaluation_interval_list = preprocessed_evaluation_intervals
    }

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: select_first([memoryMB, runtime_params.machine_mem]) + " MB"
        runtime_minutes: select_first([runtime_minutes, runtime_params.runtime_minutes])
        disks: "local-disk " + diskGB + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}