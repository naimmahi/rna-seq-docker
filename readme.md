Tools supplied in this Docker image:

(a) STAR - https://code.google.com/p/rna-star/ - GNU General Public License v3

(b) Cufflinks - http://cole-trapnell-lab.github.io/cufflinks/ - OSI-approved Boost License

(c) FASTQC - http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ - GNU General Public License v3 or later

(d) Picard - https://github.com/broadinstitute/picard - MIT License

(e) STAR pipeline from ICGC - https://github.com/akahles/icgc_rnaseq_align - MIT License

(f) Biobambam - https://github.com/gt1/biobambam - GNU General Public License v3

(g) Samtools - https://github.com/samtools/samtools - MIT License


This docker image takes as input a compressed tar archive of two fastq files and provides an aligned BAM file, QC statistics and gene expression quantification.

****Install Docker****
Instructions for the installation of Docker are available on the Docker website here: https://docs.docker.com/installation/
The Docker version currently being used is: version 1.9.0, build 76d6bc9. 


****Using the Docker image****

*   Load the Docker image. Must be done as root. 

        Command: docker load pipeline_star_cuff.tar

*   Run the Docker image as root.

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

Assuming the initial input was: /path/sample.tar
*   The output is expected to be:
        
        path/sample.tar
        path/sample_star.log
        path/star_2_pass/cufflinks_sample.log
        path/star_2_pass/sample_star.bam
        path/star_2_pass/genes.fpkm_tracking
        path/star_2_pass/isoforms.fpkm_tracking
        path/star_2_pass/skipped.gtf
        path/star_2_pass/transcripts.gtf
        path/qc/
        path/qc/fastqc_results
        path/qc/sample.validate
        path/qc/sample.rna_seq_metrics.txt


Supporting files to run the docker:

*   Reference Genome

        Download the reference genome used in the 1000 Genomes Project and decompress the reference genome.
        Link:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

        Command: gunzip <path/hs37d5.fa.gz> 


*   Gene Annotation File: annotation.gtf

        Download the Gencode Annotation file for the above reference and decompress the annotation file.
        Link:ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

        Command: gunzip <path/gencode.v19.annotation.gtf.gz>

*   Convert the chromosome names in the GTF file to reflect those in the reference genome

        Command: tail -n +6 gencode.v19.annotation.gtf | sed –e  "s/^chrM/MT/g;s/^chr//g" > gencode.v19.annotation.hs37d5_chr.gtf
        Reference: Pancancer-PCAWG Wiki

*   Build the STAR reference genome:

        STAR
        --runMode genomeGenerate 
        --genomeDir /path/star_genome/ 
        --genomeFastaFiles /path/hs37d5.fa
        --sjdbOverhang 100 
        --sjdbGTFfile /path/gencode.v19.annotation.gtf> 
        --runThreadN runThreadN

        Reference: Pancancer-PCAWG Wiki


*   Picard sequence dictionary:

        Command:
        java –jar /path/picard.jar 
        --CreateSequenceDictionary 
        R=path/hs37d5.fa 
        O=hs37d5.fa.dict


*   Samtools reference index (.faidx): 

        Command: 
        samtools faidx path/hs37d5.fa

*   RefFlat file: Convert annotation.gtf to refFlat.txt.
        
        -Download: gtfToGenePred: http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/gtfToGenePred
        -Convert: gtfToGenePred -genePredExt -geneNameAsName2 genes.gtf refFlat.tmp.txt
        paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > refFlat.txt
        gzip refFlat.txt
    
        (Source: https://gist.github.com/igordot/4467f1b02234ff864e61)


