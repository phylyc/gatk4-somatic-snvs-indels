version development


struct Runtime {
    String docker
    File? jar_override
    Int max_retries
    Int preemptible
    Int cpu
    Int machine_mem
    Int command_mem
    Int runtime_minutes
    Int disk
    Int boot_disk_size
}


struct Sample {
    File bam
    File bai
    String bam_sample_name
    String? assigned_sample_name
}
