version development

## Create a Mutect2 panel of normals
##
## Description of inputs
## normal_bams_file: file with list of normal bam paths
##      This is necessary if the array of bam paths consists of more than 16384 characters
##      because Terra is stupid.

# import "mutect2_multi_sample.wdl" as msm2
import "https://github.com/phylyc/gatk4-somatic-snvs-indels/raw/master/mutect2_pon.wdl" as m2pon

workflow Mutect2_Panel_from_File {
    input {
        File normal_bams_file
        File normal_bais_file
    }

    call m2pon.Mutect2_Panel {
        input:
            normal_bams = read_lines(normal_bams_file),
            normal_bais = read_lines(normal_bais_file)
    }

    output {
        File pon = Mutect2_Panel.pon
        File pon_idx = Mutect2_Panel.pon_idx
        Array[File] normal_calls = Mutect2_Panel.normal_calls
        Array[File] normal_calls_idx = Mutect2_Panel.normal_calls_idx
    }
}