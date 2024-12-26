import pandas as pd
import os

smp = pd.read_csv(
    "<dataDir>/rasqual/harmonised_donor_3study.txt", header=None
)[0]

chrm = sorted([str(x) for x in range(1, 23)])


rule all:
    input:
        "<dataDir>/glimpse_3study_donor/merged_chr1-22.merged.phased.vcf.gz",
        "<dataDir>/glimpse_3study_donor/merged_chr1-22.merged.1kgfmt.vcf.gz",
        "<dataDir>/glimpse_3study_donor/merged_chr1-22.merged.phased.rename.filtered.vcf.gz",
        expand(
            "<dataDir>/glimpse_3study_donor/dosage/merged_chr{chrm}.merged.phased.rename.filtered.DS.txt.gz",
            chrm=chrm,
        ),
        "<dataDir>/scatac/inferred_gender.txt",


rule merge_GL:
    input:
        covid="<dataDir>/glimpse_scATAC_donor/mpileup/merged_chr{chrm}.vcf.gz",
        benaglio="<BenaglioDataDir>/glimpse_scATAC/mpileup/merged_chr{chrm}.vcf.gz",
        you="<YouDataDir>/glimpse_scATAC/mpileup/merged_chr{chrm}.vcf.gz",
    output:
        "<dataDir>/glimpse_3study_donor/mpileup/merged_chr{chrm}.vcf.gz",
    resources:
        runtime=100,
        mem_mb=24000,
        cpus_per_task=8,
    shell:
        "bcftools merge --threads 5 -m none -Ou {input.covid} {input.benaglio} {input.you} | "
        "bcftools view --threads 5 -Oz -o {output} -S <dataDir>/rasqual/harmonised_donor_3study.txt; "
        "bcftools index -f {output}"


rule chunk:
    input:
        vcf="<dataDir>/glimpse/ref/ALL.chr{chrm}.sites.vcf.gz",
    output:
        "<dataDir>/glimpse/chunks/chunks.chr{chrm}.txt",
    resources:
        runtime=60,
        mem_mb=24000,
        cpus_per_task=8,
    shell:
        "GLIMPSE_chunk_static --thread 4 --input {input.vcf} --region {wildcards.chrm} "
        "--window-size 2000000 --buffer-size 200000 --output {output}"


rule impute:
    input:
        vcf=rules.merge_GL.output,
        ref="<dataDir>/glimpse/ref/ALL.chr{chrm}.bcf",
        gmap="/project/yangili1/zpmu/1000GP_Phase3/genetic_map_chr{chrm}_combined_b37.txt",
        chunk=rules.chunk.output,
    output:
        "<dataDir>/glimpse_3study_donor/out/merged_chr{chrm}.done.txt",
    resources:
        runtime=500,
        mem_mb=70000,
        cpus_per_task=30,
    shell:
        "while IFS='' read -r LINE || [ -n \"${{LINE}}\" ]; do "
        "printf -v ID '%02d' $(echo ${{LINE}} | cut -d' ' -f1); "
        "IRG=$(echo ${{LINE}} | cut -d' ' -f3); "
        "ORG=$(echo ${{LINE}} | cut -d' ' -f4); "
        "GLIMPSE_phase_static --thread 25 "
        "--input {input.vcf} --reference {input.ref} --map {input.gmap} "
        "--input-region ${{IRG}} --output-region ${{ORG}} "
        "--output <dataDir>/glimpse_3study_donor/out/merged_chr{wildcards.chrm}.${{ID}}.bcf; "
        "bcftools index -f <dataDir>/glimpse_3study_donor/out/merged_chr{wildcards.chrm}.${{ID}}.bcf; "
        "done < {input.chunk}; "
        "echo DONE > {output}"


rule ligate:
    input:
        rules.impute.output,
    output:
        "<dataDir>/glimpse_3study_donor/ligated/merged_chr{chrm}.merged.vcf.gz",
    resources:
        runtime=60,
        mem_mb=24000,
        cpus_per_task=10,
    shell:
        "ls <dataDir>/glimpse_3study_donor/out/merged_chr{wildcards.chrm}.*.bcf > "
        "<dataDir>/glimpse_3study_donor/ligated/merged.list.chr{wildcards.chrm}.txt; "
        "GLIMPSE_ligate_static --thread 8 "
        "--input <dataDir>/glimpse_3study_donor/ligated/merged.list.chr{wildcards.chrm}.txt "
        "--output {output}; tabix -fp vcf {output}"


rule eagle:
    input:
        vcf=rules.ligate.output,
        hap="/home/zepengmu/tools/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz",
        ref="<dataDir>/glimpse/ref/ALL.chr{chrm}.bcf",
    output:
        "<dataDir>/glimpse_3study_donor/phased/merged_chr{chrm}.merged.phased.vcf.gz",
    resources:
        runtime=180,
        mem_mb=100000,
        cpus_per_task=8,
    params:
        outputPrefix="<dataDir>/glimpse_3study_donor/phased/merged_chr{chrm}.merged.phased",
    shell:
        "eagle --geneticMapFile {input.hap} --vcfRef {input.ref} --vcfTarget {input.vcf} "
        "--outPrefix {params.outputPrefix} --allowRefAltSwap --numThreads {resources.cpus_per_task} --vcfOutFormat z; "
        "bcftools index -f {output} --thread {resources.cpus_per_task}"


rule concat:
    input:
        expand(
            "<dataDir>/glimpse_3study_donor/phased/merged_chr{chrm}.merged.phased.vcf.gz",
            chrm=chrm,
        ),
    output:
        "<dataDir>/glimpse_3study_donor/merged_chr1-22.merged.phased.vcf.gz",
    resources:
        runtime=30,
        mem_mb=80000,
        cpus_per_task=8,
    shell:
        "bcftools concat --threads 8 -Oz -o {output} {input}; tabix -fp vcf {output}"


rule rename_chr:
    input:
        rules.concat.output,
    output:
        "<dataDir>/glimpse_3study_donor/merged_chr1-22.merged.rename.vcf.gz",
    resources:
        runtime=100,
        mem_mb=24000,
        cpus_per_task=8,
    shell:
        "bcftools annotate --rename-chrs rename_chr.txt "
        "-Oz -o {output} {input}; tabix -fp vcf {output}"


rule fmt_1kg:
    input:
        rules.concat.output,
    output:
        "<dataDir>/glimpse_3study_donor/merged_chr1-22.merged.1kgfmt.vcf.gz",
    resources:
        runtime=120,
        mem_mb=24000,
        cpus_per_task=3,
    params:
        nobgzip="<dataDir>/glimpse_3study_donor/merged_chr1-22.merged.1kgfmt.vcf",
    shell:
        "python3 1kg_format.py {input} {params.nobgzip}; bgzip {params.nobgzip}; "
        "tabix -fp vcf {output}"


rule filter_vcf:
    input:
        rules.rename_chr.output,
    output:
        "<dataDir>/glimpse_3study_donor/merged_chr1-22.merged.phased.rename.filtered.vcf.gz",
    resources:
        runtime=120,
        mem_mb=30000,
        cpus_per_task=12,
    shell:
        "vcftools --gzvcf {input} "
        "--exclude-bed /project/yangili1/zpmu/ENCODE_blacklist.bed "
        "--remove-indels --recode --recode-INFO-all --stdout | "
        "bcftools view --threads 10 -S <dataDir>/rasqual/harmonised_donor_3study.txt -Ou | "
        "bcftools plugin fill-tags --threads 10 -Ou -- -t MAF,AF | "
        "bcftools filter --threads 10 --include 'INFO/INFO>=0.7 & INFO/RAF>0.05 & INFO/RAF<0.95 & MAF>0' -Oz -o {output}; "
        "tabix -fp vcf {output}"


rule get_DS:
    input:
        rules.filter_vcf.output,
    output:
        "<dataDir>/glimpse_3study_donor/dosage/merged_chr{chrm}.merged.phased.rename.filtered.DS.txt.gz",
    resources:
        runtime=120,
        mem_mb=30000,
        cpus_per_task=8,
    shell:
        "bcftools annotate -r chr{wildcards.chrm} -Ou --rename-chrs remove_chr.txt {input} | "
        "bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%DS]\\n' | bgzip -c > {output}; "
        "tabix -f -s 1 -b 2 -e 2 {output}"
