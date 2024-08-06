DEPRECATED. Please check for updated workflow [here](https://github.com/phylyc/somatic_workflow).


# somatic-snvs-indels

### Purpose : 
Workflows for somatic short variant analysis with GATK4. 

## mutect2 :
Implements somatic short variant discovery using the GATK Best Practices. If no patient-matched normal is specified, it should be expected that the variant calls of tumor samples contain ~30k rare germline variants.

#### Requirements/expectations (optional)
- tumor bam and index
- (normal bam and index)
- (Panel of Normals to filter sequencing platform artifacts)
- (germline resource (gnomAD) to filter common germline variants)
- (biallelic variants for contamination to estimate contamination between samples)
- (bwa_mem index image to filter alignment artifacts)

#### Outputs (optional)
- unfiltered vcf and index
- Mutect2 stats
- (bam file as aligned by Mutect2)
- (filtered vcf and index)
- (filtering stats)
- (read orientation bias model parameters)
- (contamination tables)
- (segmentation tables)
- (funcotated VCF or MAF)

## mutect2_multi_sample :
Implements somatic short variant discovery using the GATK Best Practices. Using multiple 
samples of the same patient increases the sensitivity for variant detection in each sample.

#### Requirements/expectations (optional)
- array of tumor bams and indices
- (array of normal bams and indices)

The other optional arguments are the same as above.

#### Outputs 
Same as above.

## mutect2_pon :
Creates a Panel of Normals to be implemented in somatic short variant discovery. 

#### Requirements/expectations
- array of normal bams and indices
- germline resource (gnomAD) to filter common germline variants

#### Outputs 
- PoN vcf and index
- normal variant call vcfs and indices

## mutect2_pon_from_file :
Creates a Panel of Normals to be implemented in somatic short variant discovery. Sometimes
it is easier to supply a file containing the path to a normal bam on each line.

#### Requirements/expectations
- file with paths to normal bams and indices
- germline resource (gnomAD) to filter common germline variants

#### Outputs 
Same as above.

## mutect2_pon_from_vcfs :
Creates a Panel of Normals to be implemented in somatic short variant discovery. This workflow can be used if variant calls from normal samples are already available.

#### Requirements/expectations
- file with paths to vcfs from Mutect2 calls on normal bams and indices
- germline resource (gnomAD) to filter common germline variants

#### Outputs 
Same as above.

## mutect2-normal-normal :
Used to validate the mutect2 workflow. Calls the mutect2 workflow on all combination of
pairs as tumor-normal pairs of the bam file and its replicates.

#### Requirements/expectations
- array of one analysis-read bam file and its replicates (and their indices)

#### Outputs
- false positive VCF files and their indices with summary 
- false positive counts

### Software version requirements :
- GATK >=4.2.6.1

Cromwell version support 
- Successfully tested on v71


## Parameter descriptions :
#### mutect2 (single pair/sample)
Primary inputs:
- ``Mutect2.intervals`` -- A file listing genomic intervals to search for somatic mutations. This should be in the standard GATK4 format.
- ``Mutect2.ref_fasta`` -- reference fasta. In google bucket:  ``gs://gatk-best-practices/somatic-b37/Homo_sapiens_assembly19.fasta``
- ``Mutect2.ref_fasta_index`` -- In google bucket: ``gs://gatk-best-practices/somatic-b37/Homo_sapiens_assembly19.fasta.fai``
- ``Mutect2.ref_dict`` -- In google bucket: ``gs://gatk-best-practices/somatic-b37/Homo_sapiens_assembly19.dict``
- ``Mutect2.tumor_bam`` -- File path or storage location (depending on backend) of the tumor bam file.
- ``Mutect2.tumor_bam_index`` -- File path or storage location (depending on backend) of the tumor bam file index.
- ``Mutect2.normal_bam`` -- (optional) File path or storage location (depending on backend) of the normal bam file.
- ``Mutect2.normal_bam_index`` -- (optional, but required if ``Mutect2.normal_bam`` is specified) File path or storage location (depending on backend) of the normal bam file index.

Workflow options:
- ``Mutect2.run_contamination_model`` -- ``true``/``false`` whether the cross-contamination model should be run.
- ``Mutect2.run_orientation_bias_mixture_model`` -- ``true``/``false`` whether the orientation bias model should be run.
- ``Mutect2.run_variant_filter`` -- ``true``/``false`` whether the Mutect2 filter model should be run. This task takes as optional input the contamination model and the orientation bias model output.
- ``Mutect2.run_realignment_filter`` -- ``true``/``false`` whether realignment filter should be run. This tasks realigns the reads at the called variant positions to hg38 and compares alignment scores.
- ``Mutect2.run_realignment_filter_only_on_high_confidence_variants`` -- ``true``/``false`` whether realignment filter should be run only on variants for which we can say with high confidence that they are somatic. The selection is determined by the ``Mutect2.select_low_conficence_variants_jexl_arg`` argument below.
- ``Mutect2.run_cnn_scoring_model`` -- ``true``/``false`` whether CNNScoreVariants should be run. This removes paired normal sample annotations from the funcotated MAF as the CNN model was trained to run on single sample VCFs. Running the CNN model currently likely leads to a failure of the workflow. 
- ``Mutect2.run_funcotator`` -- ``true``/``false`` whether Funcotator should be run. 
- ``Mutect2.compress_output`` -- ``true``/``false`` whether the vcfs should be compressed. There is no real reason not to. 
- ``Mutect2.make_bamout`` -- ``true``/``false`` whether Mutect2 should return a bam file to investigate support of called variants. (Useful for debugging.)
- ``Mutect2.keep_germline`` -- ``true``/``false`` whether germline variants should not be filtered. This is currently not supported. 
- ``Mutect2.genotype_germline_variants`` -- ``true``/``false`` whether germline variants should appear in the Mutect2 output. [Issue](https://github.com/broadinstitute/gatk/issues/7391)
- ``Mutect2.genotype_pon_sites`` -- ``true``/``false`` whether panel of normal sites should appear in the Mutect2 output. They carry a PON flag and will be filtered. 
- ``Mutect2.use_linked_de_bruijn_graph`` -- ``true``/``false`` whether to use a linked de Bruijn graph for the local assembly. This should lead to better haplotypes (gain of specificity), but is prone to obvious false negatives (loss of sensitivity) if used without recovery of all dangling branches.
- ``Mutect2.recover_all_dangling_branches`` -- ``true``/``false`` whether to recover dangling branches that fork after splitting from the reference. This is recommended to use with a linked de Bruijn graph as otherwise some highly covered known hotspot mutations (and others) may be missed.
- ``Mutect2.funcotator_use_gnomad`` -- ``true``/``false`` whether Funcotator should use gnomAD as annotation resource. 

Resources:
- ``Mutect2.panel_of_normals`` -- (optional) Panel of normals VCF to use for false positive reduction by filtering common sequencing artifacts.
- ``Mutect2.panel_of_normals_index`` -- (optional, but required if ``Mutect2.panel_of_normals`` is specified) VCF index for the panel of normals. Please see GATK4 tool ``IndexFeatureFile`` for creation of an index.
- ``Mutect2.germline_resource`` -- (optional) gnomAD vcf containing population allele frequencies (AF) of common and rare alleles. Download an exome or genome sites vcf [here](http://gnomad.broadinstitute.org/downloads). Essential for determining possible germline variants in tumor.
- ``Mutect2.germline_resource_index`` -- (optional, but required if ``Mutect2.germline_resource`` is specified) VCF index for gnomAD. Please see GATK4 tool ``IndexFeatureFile`` for creation of an index.
- ``Mutect2.variants_for_contamination`` -- (optional, but required if ``Mutect2.run_contamination_model`` is ``true``) vcf containing population allele frequencies (AF) of common SNPs. If omitted, cross-sample contamination will not be calculated and contamination filtering will not be applied. This can be generated from a gnomAD vcf using the GATK4 tool ``SelectVariants`` with the argument ``-select 'AF > 0.05'``. For speed, one can get very good results using only SNPs on chromosome 1. For example, ``java -jar $gatk SelectVariants -V gnomad.vcf -L 1 -select 'AF > 0.05' -O variants_for_contamination.vcf``.
- ``Mutect2.variants_for_contamination_index`` -- (optional, but required if ``Mutect2.variants_for_contamination`` is specified) VCF index for contamination variants. Please see GATK4 tool ``IndexFeatureFile`` for creation of an index.
- ``Mutect2.bwa_mem_index_image`` -- (optional, but required if ``Mutect2.run_realignment_filter`` is ``true``) resource for the FilterAlignmentArtifacts task. Generated by ``BwaMemIndexImageCreator``.
- ``Mutect2.funcotator_transcript_list`` -- Canonical transcript file, used when multiple possible transcripts exist.
- ``Mutect2.funcotator_data_sources_tar_gz`` -- Resource for Funcotator. If not specified, the most recent one is downloaded and used. (Specifying a resource is significantly faster.)

Extra arguments:
- ``Mutect2.split_intervals_extra_args`` -- (optional) a string of additional command line arguments of the form "-argument1 value1 -argument2 value2" for the SpitIntervals task. Most users will not need this.
- ``Mutect2.m2_extra_args`` -- (optional) a string of additional command line arguments of the form "-argument1 value1 -argument2 value2" for Mutect 2. Most users will not need this.
- ``Mutect2.m2_filter_extra_args`` -- (optional) a string of additional command line arguments of the form "-argument1 value1 -argument2 value2" for filtering of Mutect 2 calls. Most users will not need this.
- ``Mutect2.select_variants_extra_args`` -- (optional) a string of additional command line arguments of the form "-argument1 value1 -argument2 value2" for the selection of filtered Mutect 2 calls. Most users will not need this.
- ``Mutect2.select_low_conficence_variants_jexl_arg`` -- (optional) a JEXL filtering expression to select low confidence somatic variants. Only relevant if ``Mutect2.run_realignment_filter_only_on_high_confidence_variants`` is ``true``. This is useful in the tumor-only case, where we expect to call ~50k variants, most of which are rare germline variants or sequencing artifacts that were missed by the PoN or previous filtering. The default is ``"'(vc.getAttribute(\"GERMQ\") < 30) || (vc.getAttribute(\"DP\") < 4) || (vc.getAttribute(\"MBQ\").0 == 0) || (vc.getAttribute(\"MFRL\").0 == 0)'"``, filtering likely germline variants (quality score < 30), low-coverage variants, and variants without a read supporting the reference allele (likely artifacts). Those variants are not passed on to the realignment filter.
- ``Mutect2.realignment_extra_args`` -- (optional) a string of additional command line arguments of the form "-argument1 value1 -argument2 value2" for filtering realignment artifacts.  Most users will not need this.
- ``Mutect2.funcotate_extra_args`` -- (optional) a string of additional command line arguments of the form "-argument1 value1 -argument2 value2" for functional annotation using Funcotator. Most users will not need this.

Runtime options:
- ``Mutect2.scatter_count`` -- Number of parallel jobs to generate when scattering over intervals. Low scatter counts have long runtimes for each shard and increase the chance for a preemptible to fail prematurely, thus increasing cost. On the other hand, each shard has an overhead of ~5 minutes for spinup and spindown, which is the main source of compute cost for large scatter counts. The main resource difference between WES and WGS is in the filtering alignment artifacts task. Same for paired tumor-normal vs tumor-only calling. For multi-sample calling, the runtimes significantly increase. With 50 shards on 25 samples, the VariantCall task runs at least 4 hours, some shards up to 18 hours. Tumor-only cases generate ~50k variants, which, if split over only less than the default number of shards (42), require the alignment artifact filter task to use more memory than the default. TL;DR: When in doubt, use a larger scatter_count.
- ``Mutect2.gatk_docker`` -- Docker image to use for the Mutect2 tasks. This is only used for backends configured to use docker.
- ``Mutect2.gatk4_jar_override`` -- (optional) A GATK4 jar file to be used instead of the jar file in the docker image. This can be very useful for developers. Please note that you need to be careful that the docker image you use is compatible with the GATK4 jar file given here.
- ``Mutect2.preemptible`` -- Number of times to attempt running a task on a preemptible VM. This is only used for cloud backends in cromwell and is ignored for local and SGE backends. Preemptibles are very likely to fail for long tasks like VariantCall (Mutect2) or FilterAlignmentArtifacts if they take more than ~24 hours per shard. Increasing the ``scatter_count`` for those cases actually reduces costs.
- ``Mutect2.max_retries`` -- Number of times to retry failed tasks -- very important on the cloud when there are transient errors

Each task has an exposed memory argument (in MB), which can be set if the task fails due to lack of memory. The default values are optimized to work for all tested cases while keeping computing cost at a minimum.


### Functional annotation (Funcotator)

Funcotator (**FUNC**tional ann**OTATOR**) is a functional annotation tool in the core GATK toolset and was designed to handle both somatic and germline use cases. It analyzes given variants for their function (as retrieved from a set of data sources) and produces the analysis in a specified output file.  Funcotator reads in a VCF file, labels each variant with one of twenty-three distinct variant classifications, produces gene information (e.g. affected gene, predicted variant amino acid sequence, etc.), and associations to information in datasources. Default supported datasources include GENCODE (gene information and protein change prediction), dbSNP, gnomAD, and COSMIC (among others). The corpus of datasources is extensible and user-configurable and includes cloud-based datasources supported with Google Cloud Storage. Funcotator produces either a Variant Call Format (VCF) file (with annotations in the INFO field) or a Mutation Annotation Format (MAF) file.

Funcotator allows the user to add their own annotations to variants based on a set of data sources.  Each data source can be customized to annotate a variant based on several matching criteria.  This allows a user to create their own custom annotations easily, without modifying any Java code.

By default the M2 WDL runs Funcotator for functional annotation and produce a TCGA MAF from the M2 VCF.  There are several notes and caveats
- Several parameters should be passed in to populate the TCGA MAF metadata fields.  Default values are provided, though we recommend that you specify the values.  These parameters are ignored if you do not run Funcotator.
- Several fields in a TCGA MAF cannot be generated by M2 and Funcotator, such as all fields relating to validation alleles.  These will need to be populated by a downstream process created by the user.
- Funcotator does not enforce the TCGA MAF controlled vocabulary, since it is often too restrictive for general use.  This is up to the user to specify correctly.
  *Therefore, we cannot guarantee that a TCGA MAF generated here will pass the TCGA Validator*.  If you are unsure about the ramifications of this statement, then it probably does not concern you.
- More information about Funcotator can be found at: https://gatkforums.broadinstitute.org/dsde/discussion/11193/funcotator-information-and-tutorial/ 

### Important Notes :
- Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
- The provided JSON is a ready to use example JSON template of the workflow. Users are responsible for reviewing the [GATK Tool and Tutorial Documentations](https://gatk.broadinstitute.org/hc/en-us/categories/360002310591) to properly set the reference and resource variables. 
- For help running workflows on the Google Cloud Platform or locally please
view the following tutorial [(How to) Execute Workflows from the gatk-workflows Git Organization](https://gatk.broadinstitute.org/hc/en-us/articles/360035530952).
- Relevant reference and resources bundles can be accessed in [Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652).

### Contact Us :
- The following material is provided by the Data Science Platforum group at the Broad Institute. Please direct any questions or concerns to one of our forum sites : [GATK](https://gatk.broadinstitute.org/hc/en-us/community/topics) or [Terra](https://support.terra.bio/hc/en-us/community/topics/360000500432).

### LICENSING :
Copyright Broad Institute, 2020 | BSD-3

This script is released under the WDL open source code license (BSD-3) (full license text at https://github.com/openwdl/wdl/blob/master/LICENSE). Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script.
