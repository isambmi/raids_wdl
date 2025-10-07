version 1.0

# Workflow that preprocesses VCFs and then runs RAIDS
workflow RAIDS {
    input {
        File vcf_in
        String sample_id
        File ref_genotype = "gs://fc-738f2fe5-f3c7-4c79-a381-b88c5cb91cda/refGDS/matGeno1000g.gds"
        File ref_annotation = "gs://fc-738f2fe5-f3c7-4c79-a381-b88c5cb91cda/refGDS/matAnnot1000g.gds"
        File ref_fai = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        Int nb_profiles = 30

        String docker_bcftools = "quay.io/biocontainers/bcftools:1.22--h3a4d415_1"
        String docker_raids = "isabmi/raids:latest"

        Int preemptible = 1
    }

    call VcfPreprocessing {
        input:
            vcf_in=vcf_in,
            sample_id=sample_id,
        
            docker=docker_bcftools,
            preemptible=preemptible
    }

    call RunRAIDS {
        input:
            final_vcf=VcfPreprocessing.final_vcf,
            final_vcf_idx=VcfPreprocessing.final_vcf_idx,
            ref_genotype=ref_genotype,
            ref_annotation=ref_annotation,
            ref_fai=ref_fai,
            nb_profiles=nb_profiles,

            docker=docker_raids,
            preemptible=preemptible
    }

    output {
        File ancestry_csv = RunRAIDS.ancestry_csv
        File AUROC_png = RunRAIDS.AUROC_png
    }
}

task VcfPreprocessing {
    input {
        String sample_id
        File vcf_in
        String docker
        Int preemptible
    }

    command 
    <<<
        bcftools index --tbi ~{vcf_in}

        bcftools view -v snps -M2 \
        -Oz -o ~{sample_id}.vcf.gz ~{vcf_in} && \
        bcftools index -t ~{sample_id}.vcf.gz
    >>>
    runtime {
        docker: docker
        memory: "2G"
        cpu: 1
        preemptible: preemptible
    }
    output {
        File final_vcf = "~{sample_id}.vcf.gz"
        File final_vcf_idx = "~{sample_id}.vcf.gz.tbi"
    }
}

task RunRAIDS {
    input {
        File final_vcf
        File final_vcf_idx
        File ref_genotype
        File ref_annotation
        File ref_fai
        Int nb_profiles

        String docker
        Int preemptible
        
        String sample_name = basename(final_vcf, ".vcf.gz")
    }
    command {
        Rscript /opt/raids/RAIDS_script.R \
            ~{ref_genotype} \
            ~{ref_annotation} \
            ~{ref_fai} \
            ~{final_vcf} \
            ~{nb_profiles}
    }
    runtime {
        docker: docker
        cpu: 1
        preemptible: preemptible
        memory: "10G"
    }
    output {
        File ancestry_csv = "~{sample_name}_ancestry.csv"
        File AUROC_png = "~{sample_name}_AUROC.png"
    }
}