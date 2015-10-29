Tools supplied in this Docker image:

STAR - https://code.google.com/p/rna-star/ - GNU General Public License v3
Cufflinks - http://cole-trapnell-lab.github.io/cufflinks/ - OSI-approved Boost License
FASTQC - http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ - GNU General Public License v3 or later
Picard - https://github.com/broadinstitute/picard - MIT License
STAR pipeline from ICGC - https://github.com/akahles/icgc_rnaseq_align - MIT License
Biobambam - https://github.com/gt1/biobambam - GNU General Public License v3
Samtools - https://github.com/samtools/samtools - MIT License

This docker image takes as input a compressed tar archive of two fastq files and provides an aligned BAM file, QC statistics and gene expression quantification.
Using the Docker image

1.  Load the Docker image. Must be done as root. 

Command: docker load pipeline_star_cuff.tar

2.  Run the Docker image as root.

Command: docker run -v /home/ubuntu/:/host/home -v /etc:/host/etc -v /mnt/SCRATCH:/home/ubuntu/SCRATCH -i -t <docker_id> /usr/bin/python /home/ubuntu/expression/pipeline_elastic_cluster_new.py 
--analysis_id <analysis_id> 
--gtf <path/gencode.v19.annotation.hs37d5_chr.gtf> 
--p 8 
--star_pipeline /home/ubuntu/expression/icgc_rnaseq_align/star_align.py 
--input_dir <path/to/input/directory> 
--genome_fasta_file <path/hs37d5.fa>  
--genome_dir <path/to/star_genome_build> 
--quantMode TranscriptomeSAM 
--cufflinks_pipeline /home/ubuntu/expression/compute_expression.py
--ref_flat <path/to/refFlat.txt>

The successful completion of the Docker will create the following directory structure in the directory containing the input file:

Assuming the initial input was:
 </path/sample.tar>

 The output is expected to be:
 <path/sample.tar>
 <path/sample_star.log>
 <path/star_2_pass/cufflinks_sample.log>
 <path/star_2_pass/sample_star.bam>
 <path/star_2_pass/genes.fpkm_tracking>
 <path/star_2_pass/isoforms.fpkm_tracking>
 <path/star_2_pass/skipped.gtf>
 <path/star_2_pass/transcripts.gtf>
 <path/qc/>
 <path/qc/fastqc_results>
 <path/qc/sample.validate>
 <path/qc/sample.rna_seq_metrics.txt>


Supporting files to run the docker:
1. Reference Genome: reference.fa
2. Gene Annotation File: annotation.gtf
4. STAR reference genome build
5. Picard sequence dictionary: Command: java –jar </path/picard.jar> --CreateSequenceDictionary R=<path/reference.fa> O=reference.fa.dict
6. Samtools reference index (.faidx): Command: samtools faidx <path/hs37d5.fa>
7. RefFlat file: Convert annotation.gtf to refFlat.txt.
    -Download: gtfToGenePred: http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/gtfToGenePred
    -Convert: gtfToGenePred -genePredExt -geneNameAsName2 genes.gtf refFlat.tmp.txt
    paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > refFlat.txt
    gzip refFlat.txt
    (Source: https://gist.github.com/igordot/4467f1b02234ff864e61)


