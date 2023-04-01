version development


import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/util.wdl" as util


workflow CoverageIntervalsToBam {
    input {
        File ref_fasta_index
        File interval_lists_file
        String collection_name

        # runtime
        String bedtools_gatk_docker
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        Int emergency_extra_diskGB = 0

        # memory assignments in MB
        Int mem_intervals_to_bam = 2048
        Int mem_merge_bams = 8192

        # runtime assignments in minutes (for HPC cluster)
        Int time_startup = 10
        Int time_intervals_to_bam = 60
        Int time_merge_bams = 240
    }

    Runtime standard_runtime = {
        "docker": bedtools_gatk_docker,
        "jar_override": gatk_override,
        "max_retries": max_retries,
        "preemptible": preemptible,
        "cpu": 1,
        "machine_mem": 2048,
        "command_mem": 2048 - 256,
        "runtime_minutes": 60,
        "disk": 10 + emergency_extra_diskGB,
        "boot_disk_size": 12  # needs to be > 10
    }

    scatter (intervals in read_lines(interval_lists_file)) {
        call IntervalsToBam {
            input:
                interval_list = intervals,
                ref_fasta_index = ref_fasta_index,
                runtime_params = standard_runtime,
                memoryMB = mem_intervals_to_bam,
                runtime_minutes = time_intervals_to_bam + time_startup
        }
    }

    call MergeBams {
        input:
            bam_files = select_all(IntervalsToBam.bam),
            collection_name = collection_name,
            runtime_params = standard_runtime,
            memoryMB = mem_merge_bams,
            runtime_minutes = time_merge_bams + time_startup
    }

    output {
        File bam = MergeBams.bam
        File bai = MergeBams.bai
    }
}

task IntervalsToBam {
    input {
        File ref_fasta_index
        File interval_list

        Runtime runtime_params
        Int? memoryMB
        Int? runtime_minutes
        Int diskGB = 0
        Boolean use_ssd = false
    }

    parameter_meta {
#        ref_fasta_index: {localization_optional: true}
        interval_list: {localization_optional: true}
    }

    String name = basename(interval_list, ".interval_list")
    String bed = name + ".bed"
    String tmp_bed = "tmp." + bed
    String bam_out = name + ".bam"
    String tmp_bam = "tmp." + bam_out
    String bai_out = name + ".bai"

    Int disk_spaceGB = 5 * ceil(size(interval_list, "GB")) + diskGB + runtime_params.disk

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}

        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            IntervalListToBed \
            --INPUT '~{interval_list}' \
            --OUTPUT '~{tmp_bed}' \
            --SCORE 0

        # Annotate reads with sample name
        sed "s/4\t0\t+/~{name}/g" '~{tmp_bed}' > '~{bed}'
        bedtools bedtobam -i '~{bed}' -g '~{ref_fasta_index}' > '~{tmp_bam}'

        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            SortSam \
            --INPUT '~{tmp_bam}' \
            --OUTPUT '~{bam_out}' \
            --SORT_ORDER coordinate

        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            BuildBamIndex \
            --INPUT '~{bam_out}' \
            --OUTPUT '~{bai_out}'

        rm ~{bed} ~{tmp_bed} ~{tmp_bam}
    >>>

    output {
        File bam = bam_out
        File bai = bai_out
    }

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: select_first([memoryMB, runtime_params.machine_mem]) + " MB"
        runtime_minutes: select_first([runtime_minutes, runtime_params.runtime_minutes])
        disks: "local-disk " + disk_spaceGB + if use_ssd then " SSD" else " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task MergeBams {
    input {
        Array[File]+ bam_files
        String collection_name

        Runtime runtime_params
        Int? memoryMB
        Int? runtime_minutes
        Int diskGB = 0
        Boolean use_ssd = false
    }

    parameter_meta {
#        bam_files: {localization_optional: true}  # samtools requires localization
    }

    String bam_out = collection_name + ".bam"
    String bai_out = collection_name + ".bai"

    Int disk_spaceGB = 3 * ceil(size(bam_files, "GB")) + diskGB + runtime_params.disk

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            MergeSamFiles \
            ~{sep="' " prefix("-I '", bam_files)}' \
            -O `~{bam_out}' \
            --SORT_ORDER coordinate

        gatk --java-options "-Xmx~{select_first([memoryMB, runtime_params.command_mem])}m" \
            BuildBamIndex \
            -I '~{bam_out}' \
            -O '~{bai_out}' \
            --VALIDATION_STRINGENCY LENIENT
    >>>

    output {
        File bam = bam_out
        File bai = bai_out
    }

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: select_first([memoryMB, runtime_params.machine_mem]) + " MB"
        runtime_minutes: select_first([runtime_minutes, runtime_params.runtime_minutes])
        disks: "local-disk " + disk_spaceGB + if use_ssd then " SSD" else " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}