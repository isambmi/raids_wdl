version 1.0

# Workflow that processes normal BAMs into RAIDS input and runs RAIDS

workflow RAIDSFromBam {
    input {
        
        File bam
        File bam_bai

        File bam_fasta
        File bam_fasta_fai
        File bam_fasta_dict
        
        File lift_fasta
        File lift_fasta_fai
        File lift_fasta_dict
        
        String sample_id
        File intervals

        Boolean liftover = true
        File lift_chain

        File ref_genotype
        File ref_annotation
        File ref_2bit
        Int nb_profiles = 30

        File fasta = if(liftover) then lift_fasta else bam_fasta
        File fasta_fai = if(liftover) then lift_fasta_fai else bam_fasta_fai
        File fasta_dict = if(liftover) then lift_fasta_dict else bam_fasta_dict

        String docker_bcftools = "bcftools"
        String docker_gatk = "arashi-gatk-426"
        String docker_raids = "raids"
        String docker_crossmap

    }

    call generateVcf {
        input:
            fasta = bam_fasta,
            fasta_fai = bam_fasta_fai,
            fasta_dict = bam_fasta_dict,
            bam = bam,
            bam_bai = bam_bai,
            sample_id = sample_id,
            intervals = intervals,
            docker = docker_gatk
    }

    if (liftover) {
        call liftOverVcf {
            input:
                lift_chain = lift_chain,
                sample_id = sample_id,
                vcf = generateVcf.out_vcf,
                fasta = lift_fasta,
                fasta_fai = lift_fasta_fai,
                fasta_dict = lift_fasta_dict,
                docker = docker_crossmap
        }

    }

    call processVcf {
        
        input:
            sample_id = sample_id,
            vcf = select_first([liftOverVcf.lifted_vcf, generateVcf.out_vcf]),
            fasta = fasta,
            fasta_fai = fasta_fai,
            fasta_dict = fasta_dict,
            docker = docker_bcftools

    }

    call runRAIDS {
        input:
            vcf = processVcf.final_vcf,
            vcf_idx = processVcf.final_vcf_idx,
            ref_genotype = ref_genotype,
            ref_annotation = ref_annotation,
            ref_2bit = ref_2bit,
            nb_profiles = nb_profiles,
            docker = docker_raids
    }

    output {
        File final_vcf = processVcf.final_vcf
        File final_vcf_idx = processVcf.final_vcf_idx

        File ancestry_csv = runRAIDS.ancestry_csv
        File AUROC_png = runRAIDS.AUROC_png

    }
}

task generateVcf {
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File bam
        File bam_bai
        String sample_id
        File intervals
        String docker
    }
    command {
        gatk --java-options "-Xmx16g" HaplotypeCaller  \
            -R ~{fasta} \
            -I ~{bam} \
            -O ~{sample_id}.vcf.gz \
            -L ~{intervals} \
            -ERC NONE \
            --output-mode EMIT_ALL_CONFIDENT_SITES
    
    }
    runtime {
        docker: docker
    }
    output {
        File out_vcf = "~{sample_id}.vcf.gz"
    }
}

task liftOverVcf {
    input {
        File lift_chain
        String sample_id
        File vcf
        File fasta
        File fasta_fai
        File fasta_dict
        String docker
    }
    command {
        CrossMap vcf --chromid l \
            ~{lift_chain} \
            ~{vcf} \
            ~{fasta} \
            ~{sample_id}.lifted.vcf
    }
    runtime {
        docker: docker
    }
    output {
        File lifted_vcf = "~{sample_id}.lifted.vcf"
    }
}

task processVcf {
    input {
        String sample_id
        File vcf
        File fasta
        File fasta_fai
        File fasta_dict
        String docker
    }
    command {
        bcftools sort ~{vcf} -Oz -o  ~{sample_id}.sorted.vcf.gz
        bcftools index -t ~{sample_id}.sorted.vcf.gz 

        bcftools norm -m -any -f ~{fasta} \
            -Oz -o ~{sample_id}.sorted.norm.vcf.gz ~{sample_id}.sorted.vcf.gz
        bcftools index -t ~{sample_id}.sorted.norm.vcf.gz

        bcftools view -m2 -M2 \
        -i 'strlen(REF)=1 && strlen(ALT)=1 && REF~"^[ACGT]$" && ALT~"^[ACGT]$" && FMT/DP>0' \
        -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
        -Oz -o  ~{sample_id}.vcf.gz ~{sample_id}.sorted.norm.vcf.gz
        bcftools index -t ~{sample_id}.vcf.gz 
    }
    runtime {
        docker: docker
    }
    output {
        File final_vcf = "~{sample_id}.vcf.gz"
        File final_vcf_idx = "~{sample_id}.vcf.gz.tbi"
    }
}

task runRAIDS {
    input {
        File vcf
        File vcf_idx
        File ref_genotype
        File ref_annotation
        File ref_2bit
        Int nb_profiles

        String docker
        String sample_name = basename(vcf, ".vcf.gz")
    }
    command {
        Rscript /opt/raids/RAIDS_script.R \
            ~{ref_genotype} \
            ~{ref_annotation} \
            ~{ref_2bit} \
            ~{vcf} \
            ~{nb_profiles}
    }
    runtime {
        docker: docker
    }
    output {
        File ancestry_csv = "~{sample_name}_ancestry.csv"
        File AUROC_png = "~{sample_name}_AUROC.png"
    }
}
