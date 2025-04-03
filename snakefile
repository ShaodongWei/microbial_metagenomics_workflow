configfile: "config.yaml"  # Path to config file

if config["assembler"] == "spades":
    assembly = config["output_directory"] + "/.spades_assembly.done"
elif config["assembler"] == "megahit":
    assembly = config["output_directory"] + "/.megahit_assembly.done"
else:
    raise ValueError("Unknown assembler: " + config["assembler"])


rule all:
    input:
        config["output_directory"] + "/.fastp_trim.done",
        config["output_directory"] + "/.remove_host.done",
        config["output_directory"] + "/.multiqc.done",
        assembly,
        config["output_directory"] + "/.binning.done",
        config["output_directory"] + "/.metaphlan.db.done",
        config["output_directory"] + "/.metaphlan.taxonomy.done"



rule fastp_trim:
    input:
        input=config['input_directory']
    output:
        config["output_directory"] + "/.fastp_trim.done"
    params:
        threads = config["threads"],
        output = config["output_directory"]
    conda:
        "env/fastp.yaml"
    shell:
        """
        mkdir -p {params.output}/fastp_trim
        find {input} -name '*1.fastq' | while read -r r1; do
            r2=$(echo "$r1" | sed 's/1.fastq$/2.fastq/')
            basename=$(basename "$r1")
            sample_prefix=$(echo "$basename" | sed 's/_1.fastq$//')

            r1_trimmed={params.output}/fastp_trim/${{sample_prefix}}_1.trimmed.fastq
            r2_trimmed={params.output}/fastp_trim/${{sample_prefix}}_2.trimmed.fastq

            fastp -i "$r1" -I "$r2" -o "$r1_trimmed" -O "$r2_trimmed" --thread {params.threads} --json {params.output}/fastp_trim/"$sample_prefix"_fastp.json

        done
        touch {output}
        """

rule remove_host:
    input:
        config["output_directory"] + "/.fastp_trim.done"
    output:
        config["output_directory"] + "/.remove_host.done"
    params:
        threads = config["threads"],
        input = config['input_directory'],
        output = config["output_directory"],
        host_genome = config["host_genome"]
    conda:
        "env/remove_host.yaml"
    shell:
        """
        mkdir -p {params.output}/clean_fastq
        mkdir -p {params.output}/clean_fastq/bowtie2_log

        # Build Bowtie2 index if not already present
        if [ ! -f {params.host_genome}.1.bt2 ]; then
            bowtie2-build {params.host_genome} {params.host_genome}
        fi
      
        find {params.output}/fastp_trim -name '*1.trimmed.fastq' | while read -r r1; do
            r2=$(echo "$r1" | sed 's/1.trimmed.fastq$/2.trimmed.fastq/')
            basename=$(basename "$r1")
            sample_prefix=$(echo "$basename" | sed 's/_1.trimmed.fastq$//')

            bowtie2 -x {params.host_genome} -1 $r1 -2 $r2 -p {params.threads} --very-sensitive\
            --un-conc {params.output}/clean_fastq/"$sample_prefix"_unmapped.fastq -S /dev/null \
            2> {params.output}/clean_fastq/bowtie2_log/"$sample_prefix".bowtie2.log

            mv {params.output}/clean_fastq/"$sample_prefix"_unmapped.1.fastq {params.output}/clean_fastq/"$sample_prefix"_1.fastq
            mv {params.output}/clean_fastq/"$sample_prefix"_unmapped.2.fastq {params.output}/clean_fastq/"$sample_prefix"_2.fastq

        done
        touch {output}
        """

rule multiqc:
    input:
        config["output_directory"] + "/.remove_host.done"
    output:
        config["output_directory"] + "/.multiqc.done"
    params:
        output = config["output_directory"],
    conda:
        "env/multiqc.yaml"
    shell:
        """
        mkdir -p {params.output}/multiqc

        multiqc {params.output} -o {params.output}/multiqc

        touch {output}
        """


rule spades_assembly:
    input:
        done = config["output_directory"] + "/.multiqc.done"
    output:
        config["output_directory"] + "/.spades_assembly.done"
    params:
        output = config["output_directory"],
        assembler = config["assembler"],
        threads = config["threads"],
        assembler_mode = config['assembler_mode']
    conda:
        "env/spades.yaml"

    shell:
        """
        mkdir -p {params.output}/assembly
        find {params.output}/clean_fastq -name '*1.fastq' | while read -r r1; do
            r2=$(echo "$r1" | sed 's/1.fastq$/2.fastq/')
            basename=$(basename "$r1")
            sample_prefix=$(echo "$basename" | sed 's/_1.fastq$//')

            spades.py -1 $r1 -2 $r2 -o {params.output}/assembly/"$sample_prefix" -t {params.threads} --{params.assembler_mode}

        done

        touch {output}
        """

rule megahit_assembly:
    input:
        done = config["output_directory"] + "/.multiqc.done"
    output:
        config["output_directory"] + "/.megahit_assembly.done"
    params:
        output = config["output_directory"],
        assembler = config["assembler"],
        threads = config["threads"],
    conda:
        "env/megahit.yaml"

    shell:
        """
        mkdir -p {params.output}/assembly
        find {params.output}/clean_fastq -name '*1.fastq' | while read -r r1; do
            r2=$(echo "$r1" | sed 's/1.fastq$/2.fastq/')
            basename=$(basename "$r1")
            sample_prefix=$(echo "$basename" | sed 's/_1.fastq$//')
            megahit -1 $r1 -2 $r2 -o {params.output}/assembly/"$sample_prefix" -t {params.threads} > {params.output}/assembly/"$sample_prefix".megahit.log 2>&1
        done

        touch {output}
        """

rule binning:
    input:
        assembly
    output:
        config["output_directory"] + "/.binning.done"
    params:
        output = config["output_directory"],
        threads = config["threads"],
        assembler = config["assembler"]
    conda:
        "env/comebin.yaml"
    shell:
        """
        mkdir -p {params.output}/binning

        if [[ -f {params.output}/samples.txt ]]; then rm {params.output}/samples.txt; fi
        find {params.output}/clean_fastq -name "*1.fastq" | while read -r file; do basename $file | cut -d_ -f1 | sort | uniq >> {params.output}/samples.txt;done 


        while read sample; do
            echo "Processing $sample"

            if [[ {params.assembler} == 'spades' ]]; then

                # Alignment
                minimap2 -ax sr -t {params.threads} \
                    {params.output}/assembly/$sample/scaffolds.fasta \
                    {params.output}/clean_fastq/"$sample"_1.fastq \
                    {params.output}/clean_fastq/"$sample"_2.fastq \
                    > {params.output}/binning/"$sample".alignment.sam

                # convert sam to bam 
                samtools view -hb {params.output}/binning/"$sample".alignment.sam | \
                samtools sort -o {params.output}/binning/"$sample".alignment.sorted.bam

                # delete sam file 
                rm {params.output}/binning/"$sample".alignment.sam

                # combebin -p is directory,
                run_comebin.sh -a {params.output}/assembly/$sample/scaffolds.fasta \
                -o  {params.output}/binning/$sample\
                -p {params.output}/binning  \
                -t {params.threads}
            else
                # Alignment
                minimap2 -ax sr -t {params.threads} \
                    {params.output}/assembly/$sample/final.contigs.fa \
                    {params.output}/clean_fastq/"$sample"_1.fastq \
                    {params.output}/clean_fastq/"$sample"_2.fastq \
                    > {params.output}/binning/"$sample".alignment.sam

                # convert sam to bam 
                samtools view -hb {params.output}/binning/"$sample".alignment.sam | \
                samtools sort -o {params.output}/binning/"$sample".alignment.sorted.bam

                # delete sam file 
                rm {params.output}/binning/"$sample".alignment.sam

                # combebin -p is directory,
                run_comebin.sh -a {params.output}/assembly/$sample/final.contigs.fa \
                -o  {params.output}/binning/$sample\
                -p {params.output}/binning  \
                -t {params.threads}
            fi

        done < {params.output}/samples.txt

        # rm -f samples.txt
        touch {output}
        """



# rule binning:
#     input:
#         assembly
#     output:
#         config["output_directory"] + "/.binning.done"
#     params:
#         output = config["output_directory"],
#         threads = config["threads"],
#         assembler = config["assembler"]
#     conda:
#         "env/metadecoder.yaml"
#     shell:
#         """
#         mkdir -p {params.output}/binning

#         if [[ -f {params.output}/samples.txt ]]; then rm {params.output}/samples.txt; fi
#         find {params.output}/clean_fastq -name "*1.fastq" | while read -r file; do basename $file | cut -d_ -f1 | sort | uniq >> {params.output}/samples.txt;done 

#         # not knowing why the dependency fraggenescan has no permission 
#         #chmod +x $(pwd)/.snakemake/conda/fccdad0af2863466fb6cf6d8a303b12d_/lib/python3.13/site-packages/metadecoder/fraggenescan
#         chmod +x $(pwd)/.snakemake/conda/8b37a24dbb6a344e5a1b463c89fe6e71_/lib/python3.13/site-packages/metadecoder/fraggenescan 
        
#         while read sample; do
#             echo "Processing $sample"
#             mkdir -p {params.output}/binning/$sample

#             if [[ {params.assembler} == 'spades' ]]; then
#                 # Alignment
#                 minimap2 -ax sr -t {params.threads} \
#                     {params.output}/assembly/$sample/scaffolds.fasta \
#                     {params.output}/clean_fastq/"$sample"_1.fastq \
#                     {params.output}/clean_fastq/"$sample"_2.fastq \
#                     > {params.output}/binning/"$sample"/"$sample".alignment.sam

#                 # convert sam to bam 
#                 samtools view -hb {params.output}/binning/"$sample"/"$sample".alignment.sam | \
#                 samtools sort -o {params.output}/binning/"$sample"/"$sample".alignment.sorted.bam

#                 # delete sam file 
#                 rm {params.output}/binning/"$sample"/"$sample".alignment.sam

#                 # Coverage
#                 metadecoder coverage \
#                     -b {params.output}/binning/"$sample"/"$sample".alignment.sorted.bam \
#                     -o {params.output}/binning/"$sample"/"$sample".metadecoder.coverage \
#                     --threads {params.threads}

#                 # Seed
#                 metadecoder seed \
#                     -f {params.output}/assembly/"$sample"/scaffolds.fasta \
#                     -o {params.output}/binning/"$sample"/"$sample".metadecoder.seed \
#                     --threads {params.threads}
                
#                 # Cluster
#                     mkdir -p {params.output}/binning/bins
#                     metadecoder cluster \
#                         -f {params.output}/assembly/"$sample"/scaffolds.fasta \
#                         -c {params.output}/binning/"$sample"/"$sample".metadecoder.coverage \
#                         -s {params.output}/binning/"$sample"/"$sample".metadecoder.seed \
#                         -o {params.output}/binning/bins/"$sample".bins.tsv

#             else
#                 # Alignment
#                 # minimap2 -ax sr -t {params.threads} \
#                 #     {params.output}/assembly/$sample/final.contigs.fa \
#                 #     {params.output}/clean_fastq/"$sample"_1.fastq \
#                 #     {params.output}/clean_fastq/"$sample"_2.fastq \
#                 #     > {params.output}/binning/"$sample"/"$sample".alignment.sam

#                 # # convert sam to bam 
#                 # samtools view -hb {params.output}/binning/"$sample"/"$sample".alignment.sam | \
#                 # samtools sort -o {params.output}/binning/"$sample"/"$sample".alignment.sorted.bam

#                 # # delete sam file 
#                 # rm {params.output}/binning/"$sample"/"$sample".alignment.sam

#                 # # Coverage
#                 # metadecoder coverage \
#                 #     -b {params.output}/binning/"$sample"/"$sample".alignment.sorted.bam \
#                 #     -o {params.output}/binning/"$sample"/"$sample".metadecoder.coverage \
#                 #     --threads {params.threads}

#                 # Seed
#                 # which FragGeneScan
#                 # which fraggenescan
#                 # which metadecoder
#                 # cd $(dirname $(which FragGeneScan))
#                 # FragGeneScan -s /data/output/assembly/a1/final.contigs.fa -o test_fraggenescan -w 0 -t complete
#                 # ln -s $(dirname $(which FragGeneScan))/train {params.output}/binning/$sample/train

#                 # metadecoder seed \
#                 #     -f {params.output}/assembly/"$sample"/final.contigs.fa \
#                 #     -o {params.output}/binning/"$sample"/"$sample".metadecoder.seed \
#                 #     --threads {params.threads}
                
#                 # # Cluster
#                 # mkdir -p {params.output}/binning/bins
#                 # metadecoder cluster \
#                 #     -f {params.output}/assembly/"$sample"/final.contigs.fa \
#                 #     -c {params.output}/binning/"$sample"/"$sample".metadecoder.coverage \
#                 #     -s {params.output}/binning/"$sample"/"$sample".metadecoder.seed \
#                 #     -o {params.output}/binning/bins/"$sample".bins.tsv
                    
#             fi
        
#         done < {params.output}/samples.txt

#         # rm -f samples.txt
#         touch {output}
#         """

rule download_metaphlan_db:
    output:
        config["output_directory"] + "/.metaphlan.db.done"
    conda:
        "env/metaphlan.yaml"
    params:
        metaphlan_db = config['metaphlan_db']
    shell:
        """
        mkdir -p {params.metaphlan_db}
        if [[ -z $(ls -A {params.metaphlan_db}) ]]; then 
            metaphlan --install --bowtie2db {params.metaphlan_db}
        fi
        touch {output}
        """

rule taxonomy:
    input:
        done = config["output_directory"] + "/.megahit_assembly.done"
    output:
        config["output_directory"] + "/.metaphlan.taxonomy.done"
    params:
        output = config["output_directory"],
        threads = config["threads"],
        metaphlan_db = config['metaphlan_db']
    conda:
        "env/metaphlan.yaml"
    shell:
        """
        mkdir -p {params.output}/taxonomy

        find {params.output}/clean_fastq -name '*1.fastq' | while read -r r1; do
            r2=$(echo "$r1" | sed 's/1.fastq$/2.fastq/')
            basename=$(basename "$r1")
            sample_prefix=$(echo "$basename" | sed 's/_1.fastq$//')
            
            # seems there is memory problem
            metaphlan $r1,$r2 \
                --input_type fastq \
                --bowtie2db {params.metaphlan_db} \
                --bowtie2out {params.output}/taxonomy/"$sample_prefix".bowtie2.bz2 \
                --nproc {params.threads} \
                -o {params.output}/taxonomy/"$sample_prefix".profile.txt \
                > {params.output}/taxonomy/"$sample_prefix".metaphlan.log 2>&1
        done

        touch {output}
        """


