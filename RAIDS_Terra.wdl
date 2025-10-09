version 1.0

# Workflow that preprocesses VCFs and then runs RAIDS
workflow RAIDS {
    input {
        File vcf_in
        String sample_id
        String normal_sample_id = "NORMAL"
        File ref_genotype = "gs://fc-738f2fe5-f3c7-4c79-a381-b88c5cb91cda/refGDS/matGeno1000g.gds"
        File ref_annotation = "gs://fc-738f2fe5-f3c7-4c79-a381-b88c5cb91cda/refGDS/matAnnot1000g.gds"
        File ref_fai = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        Int nb_profiles = 30

        String docker_bcftools = "us-central1-docker.pkg.dev/image-rchiv/gdc-raids/bcftools:1.22--h3a4d415_1"
        String docker_raids = "us-central1-docker.pkg.dev/image-rchiv/gdc-raids/raids:generic"

        Int preemptible = 1
        Int raids_cpu = 1
        String raids_mem = "25G"
        String raids_disk = "local-disk 50 HDD"
    }

    call VcfPreprocessing {
        input:
            vcf_in=vcf_in,
            sample_id=sample_id,
            normal_sample_id=normal_sample_id,
        
            docker=docker_bcftools,
            preemptible=preemptible
    }

    call RunRAIDS {
        input:
            input_csv=VcfPreprocessing.csv,
            ref_genotype=ref_genotype,
            ref_annotation=ref_annotation,
            ref_fai=ref_fai,
            nb_profiles=nb_profiles,

            docker=docker_raids,
            preemptible=preemptible,
            raids_cpu=raids_cpu,
            raids_mem=raids_mem,
            raids_disk=raids_disk
    }

    output {
        File ancestry_csv = RunRAIDS.ancestry_csv
        File AUROC_png = RunRAIDS.AUROC_png
    }
}

task VcfPreprocessing {
    input {
        File vcf_in
        String sample_id
        String normal_sample_id
        String docker
        Int preemptible
    }

    command 
    <<<
        # create index
        bcftools index --tbi ~{vcf_in}

        # filter for multiallelic non-SNPs
        bcftools view -v snps -m2 -M2 \
        -s ~{normal_sample_id} \
        -Oz -o tmp.vcf.gz ~{vcf_in} && \
        bcftools index -t tmp.vcf.gz

        # building generic format file from VCF
        bcftools query -s ~{normal_sample_id} -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' tmp.vcf.gz \
        | awk 'BEGIN{
                OFS=","; print "Chromosome,Position,Ref,Alt,Count,File1R,File1A"
            }
            {
                split($5,a,","); Chrom=$1; Pos0=$2-1; Ref=$3; Alt=$4;
                File1R=a[1]; File1A=a[2]; Count=File1R+File1A;
                print Chrom,Pos0,Ref,Alt,Count,File1R,File1A
            }' > ~{sample_id}.csv


    >>>
    runtime {
        docker: docker
        memory: "2G"
        cpu: 1
        preemptible: preemptible
    }
    output {
        File csv = "~{sample_id}.csv"
    }
}

task RunRAIDS {
    input {
        File input_csv
        File ref_genotype
        File ref_annotation
        File ref_fai
        Int nb_profiles

        String docker
        Int preemptible
        Int raids_cpu
        String raids_mem
        String raids_disk
        
        String sample_name = basename(input_csv, ".csv")
    }
    command {
        Rscript /opt/raids/RAIDS_script.R \
            ~{ref_genotype} \
            ~{ref_annotation} \
            ~{ref_fai} \
            ~{input_csv} \
            ~{nb_profiles}
    }
    runtime {
        docker: docker
        cpu: raids_cpu
        preemptible: preemptible
        memory: raids_mem
        disks: raids_disk
    }
    output {
        File ancestry_csv = "~{sample_name}_ancestry.csv"
        File AUROC_png = "~{sample_name}_AUROC.png"
    }
}