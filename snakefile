configfile: "config.yaml"  # Path to config file

if config["assembler"] == "spades":
    assembly = config["output_directory"] + "/.spades_assembly.done"
elif config["assembler"] == "megahit":
    assembly = config["output_directory"] + "/.megahit_assembly.done"
else:
    raise ValueError("Unknown assembler: " + config["assembler"])

if config['run_binning']:
    binning = config["output_directory"] + "/.semibin.done"
else: 
    binning = None

if config['functional_profiling']:
    functional_profiling = config["output_directory"] + "/.humann.done"
    human_db = config["output_directory"] + "/.humann.db.done"
else: 
    functional_profiling = None
    human_db = None

rule all:
    input:
        config["output_directory"] + "/.fastp_trim.done",
        config["output_directory"] + "/.remove_host.done",
        config["output_directory"] + "/.multiqc.done",
        assembly,
        binning if config['run_binning'] else [],
        config["output_directory"] + "/.metaphlan.db.done",
        config["output_directory"] + "/.metaphlan.taxonomy.done",
        functional_profiling if config['functional_profiling'] else [],
        human_db if config['functional_profiling'] else []
        
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
            mv {params.output}/assembly/"$sample_prefix".megahit.log {params.output}/assembly/"$sample_prefix"
        done

        touch {output}
        """

rule binning:
    input:
        assembly
    output:
        config["output_directory"] + "/.semibin.done"
    params:
        output = config["output_directory"],
        threads = config["threads"],
        assembler = config["assembler"]
    conda:
        "env/semibin.yaml"
    shell:
        """
        mkdir -p {params.output}/binning

        if [[ -f {params.output}/samples.txt ]]; then rm {params.output}/samples.txt; fi
        find {params.output}/clean_fastq -name "*1.fastq" | while read -r file; do basename $file | cut -d_ -f1 | sort | uniq >> {params.output}/samples.txt;done 

        while read sample; do

            if [[ {params.assembler} == 'spades' ]]; then
                
                mkdir -p {params.output}/binning/$sample

                # Alignment
                minimap2 -ax sr -t {params.threads} \
                    {params.output}/assembly/$sample/scaffolds.fasta \
                    {params.output}/clean_fastq/"$sample"_1.fastq \
                    {params.output}/clean_fastq/"$sample"_2.fastq \
                    > {params.output}/binning/"$sample"/"$sample".alignment.sam
                
                # convert sam to bam 
                samtools view -hb {params.output}/binning/"$sample"/"$sample".alignment.sam | \
                samtools sort -o {params.output}/binning/"$sample"/"$sample".alignment.sorted.bam

                # delete sam file 
                rm {params.output}/binning/"$sample"/"$sample".alignment.sam

                
                SemiBin2 single_easy_bin --input-fasta {params.output}/assembly/$sample/scaffolds.fasta --input-bam {params.output}/binning/"$sample"/"$sample".alignment.sorted.bam --environment human_gut --output {params.output}/binning/$sample

            else

                mkdir -p {params.output}/binning/$sample

                # Alignment
                minimap2 -ax sr -t {params.threads} \
                    {params.output}/assembly/$sample/final.contigs.fa \
                    {params.output}/clean_fastq/"$sample"_1.fastq \
                    {params.output}/clean_fastq/"$sample"_2.fastq \
                    > {params.output}/binning/"$sample"/"$sample".alignment.sam

                # convert sam to bam 
                samtools view -hb {params.output}/binning/"$sample"/"$sample".alignment.sam | \
                samtools sort -o {params.output}/binning/"$sample"/"$sample".alignment.sorted.bam

                # delete sam file 
                #rm {params.output}/binning/"$sample"/"$sample".alignment.sam

                
                SemiBin2 single_easy_bin --input-fasta {params.output}/assembly/$sample/final.contigs.fa --input-bam {params.output}/binning/"$sample"/"$sample".alignment.sorted.bam --environment human_gut --output {params.output}/binning/$sample

            fi
            
            
            
        done < {params.output}/samples.txt

        rm -f {params.output}/samples.txt
        touch {output}
        """

rule download_metaphlan_db:
    input:
        assembly
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
        config["output_directory"] + "/.metaphlan.db.done"
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
            
            mkdir -p {params.output}/taxonomy/$sample_prefix
            # seems there is memory problem
            metaphlan $r1,$r2 \
                --input_type fastq \
                --bowtie2db {params.metaphlan_db} \
                --bowtie2out {params.output}/taxonomy/$sample_prefix/"$sample_prefix".bowtie2.bz2 \
                --nproc {params.threads} \
                -o {params.output}/taxonomy/$sample_prefix/"$sample_prefix".profile.txt \
                > {params.output}/taxonomy/$sample_prefix/"$sample_prefix".metaphlan.log 2>&1
        done

        touch {output}
        """


rule download_humann_db:
    input:
        config["output_directory"] + "/.metaphlan.taxonomy.done"
    output:
        config["output_directory"] + "/.humann.db.done"
    conda:
        "env/humann.yaml"
    params:
        humann_db = config['humann_db']
    shell:
        """
        mkdir -p {params.chocophlan_db}
        if [[ -z $(ls -A {params.chocophlan_db}) ]]; then 
            humann_databases --download chocophlan full {params.humann_db}
            humann_databases --download uniref uniref90_diamond {params.humann_db}

        fi
        touch {output}
        """


rule functional_profiling:
    input:
        config["output_directory"] + "/.humann.chocophlan.done"
    output:
        config["output_directory"] + "/.humann.done"
    params:
        output = config["output_directory"],
        threads = config["threads"],
        chocophlan_db = config['chocophlan_db']
    conda:
        "env/humann.yaml"
    shell:
        """
        mkdir -p {params.output}/humann

        find {params.output}/clean_fastq -name '*1.fastq' | while read -r r1; do
            r2=$(echo "$r1" | sed 's/1.fastq$/2.fastq/')
            basename=$(basename "$r1")
            sample_prefix=$(echo "$basename" | sed 's/_1.fastq$//')

            mkdir -p {params.output}/humann/$sample_prefix
            cat "$r1" "$r2" > {params.output}/humann/$sample_prefix/concatenated_fastq
            humann \
                --input {params.output}/humann/$sample_prefix/concatenated_fastq \
                --output {params.output}/humann/$sample_prefix \
                --threads {params.threads} \
                --nucleotide-database {params.humann_db}/chocophlan \
                --protein-database {params.humann_db}/uniref \
                > {params.output}/humann/$sample_prefix/"$sample_prefix".humann.log 2>&1
            rm {params.output}/humann/$sample_prefix/concatenated_fastq
        done

        touch {output}
        """
