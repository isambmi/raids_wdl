version 1.0

# Workflow that preprocesses VCFs and then runs RAIDS
workflow RAIDS {
    input {
        String vcf_ref = "hg38"
        File? chr_map
        File vcf_in
        File vcf_in_idx
        String sample_id
        File? lift_chain
        File? target_ref_fa
        File? target_ref_fai
        File? target_ref_dict
        File ref_genotype = "gs://fc-738f2fe5-f3c7-4c79-a381-b88c5cb91cda/refGDS/matGeno1000g.gds"
        File ref_annotation = "gs://fc-738f2fe5-f3c7-4c79-a381-b88c5cb91cda/refGDS/matAnnot1000g.gds"
        File ref_fai = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        Int nb_profiles = 30

        String docker_bcftools = "biocontainers/bcftools:v1.9-1-deb_cv1"
        String docker_gatk = "broadinstitute/gatk:4.6.2.0"
        String docker_raids = "isabmi/raids:latest"

        Boolean run_TN = false
    }

    if (vcf_ref != "hg38") {
        call ReannotateVcf {
            input:
                chr_map=chr_map,
                sample_id=sample_id,
                vcf_in=vcf_in,
                vcf_in_idx=vcf_in_idx,
                
                docker=docker_bcftools
        }

        call LiftOverVcf {
            input:
                input_vcf = ReannotateVcf.reannotated_vcf,
                input_vcf_idx = ReannotateVcf.reannotated_vcf_idx,
                sample_id = sample_id,
                lift_chain = lift_chain,
                target_ref_fa = target_ref_fa,
                target_ref_fai = target_ref_fai,
                target_ref_dict = target_ref_dict,
                
                docker = docker_gatk
        }
    }

    File pre_raids_vcf = select_first([LiftOverVcf.lifted_vcf, vcf_in])
    File pre_raids_idx = select_first([LiftOverVcf.lifted_idx, vcf_in_idx])

    call VcfPreprocessing {
        input:
            sample_id = sample_id,
            lifted_vcf = pre_raids_vcf,
            lifted_idx = pre_raids_idx,
        
            docker = docker_bcftools
    }

    call GetTNFromHeader {
        input:
            vcf_gz = VcfPreprocessing.final_vcf,
            docker = docker_bcftools
    }

    Boolean haveTN = (GetTNFromHeader.tumor_name != "") && (GetTNFromHeader.normal_name != "")

    if (haveTN && run_TN) {
        Array[String] sub_ids = read_lines(VcfPreprocessing.samples_list)

        scatter (sid in sub_ids) {
                String label = if (sid == GetTNFromHeader.tumor_name) then "tumor"
                    else if (sid == GetTNFromHeader.normal_name) then "normal"
                    else "unknown"

            call SubSample as SubSample_TN {
                input:
                    label = label,
                    sub_id = sid,
                    sample_id = sample_id,
                    final_vcf = VcfPreprocessing.final_vcf,
                    final_vcf_idx = VcfPreprocessing.final_vcf_idx,
                    docker = docker_bcftools
            }
            call RunRAIDS as RAIDS_TN {
                input:
                    final_vcf = SubSample_TN.sub_vcf,
                    final_vcf_idx = SubSample_TN.sub_vcf_idx,
                    ref_genotype = ref_genotype,
                    ref_annotation = ref_annotation,
                    ref_fai = ref_fai,
                    nb_profiles = nb_profiles,

                    docker = docker_raids
            }
        }
    }

    call RunRAIDS {
        input:
            final_vcf = VcfPreprocessing.final_vcf,
            final_vcf_idx = VcfPreprocessing.final_vcf_idx,
            ref_genotype = ref_genotype,
            ref_annotation = ref_annotation,
            ref_fai = ref_fai,
            nb_profiles = nb_profiles,

            docker = docker_raids
    }

    output {
        File germline_sub30 = VcfPreprocessing.germline_sub30
        File germline_filter = VcfPreprocessing.germline_filter
        File samples_list = VcfPreprocessing.samples_list

        File ancestry_csv = RunRAIDS.ancestry_csv
        File AUROC_png = RunRAIDS.AUROC_png

        Array[File]? TN_ancestry_csv = if (haveTN) then RAIDS_TN.ancestry_csv else []
        Array[File]? TN_AUROC_png = if (haveTN) then RAIDS_TN.AUROC_png else []
    }
}

task ReannotateVcf {
    input {
        File? chr_map
        String sample_id
        File vcf_in
        File vcf_in_idx
        String docker

    }
    command {
        bcftools annotate --rename-chrs ~{chr_map} \
            -Oz -o ~{sample_id}.hg19.vcf.gz ~{vcf_in} && \
        bcftools index -t ~{sample_id}.hg19.vcf.gz
    }
    runtime {
        docker: docker
        memory: "100M"
        cpu: 1
        preemptible: 1
    }
    output {
        File reannotated_vcf = "~{sample_id}.hg19.vcf.gz"
        File reannotated_vcf_idx = "~{sample_id}.hg19.vcf.gz.tbi"
    }

}

task LiftOverVcf {
    input {
        File input_vcf
        File input_vcf_idx
        String sample_id
        File? lift_chain
        File? target_ref_fa
        File? target_ref_fai
        File? target_ref_dict

        String docker
        
    }
    command {
        gatk --java-options "-Xmx4500m" LiftoverVcf \
            --I ~{input_vcf} \
            -O ~{sample_id}.hg38.vcf.gz\
            -CHAIN ~{lift_chain} \
            --REJECT ~{sample_id}.rejects.hg19tohg38.vcf.gz \
            -R ~{target_ref_fa}
    }
    runtime {
        docker: docker
        memory: "10G"
        cpu: 1
        preemptible: 1
    }
    output {
        File lifted_vcf = "~{sample_id}.hg38.vcf.gz"
        File lifted_idx = "~{sample_id}.hg38.vcf.gz.tbi"
        File liftover_rejects = "~{sample_id}.rejects.hg19tohg38.vcf.gz"
    }
}

task VcfPreprocessing {
    input {
        String sample_id
        File lifted_vcf
        File lifted_idx
        
        String docker
        
    }
    command {
        bcftools view -v snps -M2 \
        -Oz -o ~{sample_id}.vcf.gz ~{lifted_vcf} && \
        bcftools index -t ~{sample_id}.vcf.gz

        # germline sanity check
        bcftools view -i 'INFO/GERMQ < 30' ~{sample_id}.vcf.gz -H > ~{sample_id}_sub30germq.txt
        bcftools view -f germline -H ~{sample_id}.vcf.gz > ~{sample_id}_germline.txt

        # creating T/N VCFs too just in case
        bcftools query -l ~{sample_id}.vcf.gz > ~{sample_id}_samples.txt

    }
    runtime {
        docker: docker
        memory: "100M"
        cpu: 1
        preemptible: 1
    }
    output {
        File final_vcf = "~{sample_id}.vcf.gz"
        File final_vcf_idx = "~{sample_id}.vcf.gz.tbi"
        File germline_sub30 = "~{sample_id}_sub30germq.txt"
        File germline_filter = "~{sample_id}_germline.txt"
        File samples_list = "~{sample_id}_samples.txt"
    }
}

task GetTNFromHeader {
  input {
    File vcf_gz
    String docker = "bcftools"
    
  }
  command <<<
    set -euo pipefail

    bcftools view -h ~{vcf_gz} > header.txt

    awk -F'[=,>[:space:]]+' '/^##tumor_sample=/{print $2; exit}'  header.txt > tumor_name.txt || true
    awk -F'[=,>[:space:]]+' '/^##normal_sample=/{print $2; exit}' header.txt > normal_name.txt || true

    [ -s tumor_name.txt ]  || : > tumor_name.txt
    [ -s normal_name.txt ] || : > normal_name.txt
  >>>
  runtime { 
        docker: docker 
        cpu: 1
        preemptible: 1
        memory: "100M"
    }
  output {
    String tumor_name  = read_string("tumor_name.txt")
    String normal_name = read_string("normal_name.txt")
  }
}

task SubSample {
    input {
        String label
        String sub_id
        String sample_id
        File final_vcf
        File final_vcf_idx

        String docker
        
    }
    command {
        bcftools view -s ~{sub_id} -Oz -o ~{sample_id}.~{label}.vcf.gz ~{final_vcf} && \
        bcftools index -t ~{sample_id}.~{label}.vcf.gz
    }
    runtime {
        docker: docker
        cpu: 1
        preemptible: 1
        memory: "100M"
    }
    output {
        File sub_vcf = "~{sample_id}.~{label}.vcf.gz"
        File sub_vcf_idx = "~{sample_id}.~{label}.vcf.gz.tbi"
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
        preemptible: 1
        memory: "10G"
    }
    output {
        File ancestry_csv = "~{sample_name}_ancestry.csv"
        File AUROC_png = "~{sample_name}_AUROC.png"
    }
}
