#!/usr/bin/env nextflow

// aux
params.gt = "$projectDir/lofreq_gt.sh"
params.threads = 80

// FASTQs
params.f1 = " "
params.f2 = " "

// Reference genomes
params.genome = " "
params.ugenome = " "
params.hg38 = " "

// Known variants/databases
params.known = " "
params.cosmic = " "
params.dbsnp = " "


// ML predictions
params.model = " "
params.output = " "
params.threshold = 0.75






// Step 1: Data preprocessing
process bwaMem {
    input:
    tuple val(sampleId), file(reads)

    output:
    path "${sampleId}.sort.bam" 

    script:
    """
    python $projectDir/code/bwa.py -r ${params.hg38} -t ${params.threads} -f1 ${params.f1} -f2 ${params.f2}
    """
}

process indexVCFs{
    publishDir "$projectDir/databases", mode: 'symlink'

    input:
    file(vcf)

    output:
    path "${vcf}.tbi"

    script:
    """
    tabix -p vcf ${vcf}
    """
}



process bamPreprocessing {
    publishDir "$projectDir/code/", mode: 'symlink'
    input:
    file(bam)


    output:
    path "*.rg.bam"

    script :
    """
    python $projectDir/code/gatk.py -f ${params.genome} -k ${params.known}
    """
}


process bqsrProcess {
    publishDir "$projectDir/code/", mode: 'symlink'
    input:
    file(bam)

    output:
    path "*.table"

    script :
    """
    python $projectDir/code/bqsr.py -f ${params.genome} -k ${params.known} -b ${bam}
    """
}


process applyBQSR {
    publishDir "$projectDir/code/", mode: 'symlink'
    input:
    file(bam)
    file(table)

    output:
    path "*.bam"

    script :
    """
    python $projectDir/code/apply.py -f ${params.genome} -k ${params.known} -b ${bam} -t ${table}
    """
}


// Step 2: Variant calling
process bcftoolsCall {
    publishDir "$projectDir/code/", mode: 'symlink'

    input:
    file(bam)

    output:
    path "bcftools.*.vcf.gz"

    script:
    """
    python $projectDir/code/bcftools.py -f ${params.genome} -b ${bam}
    """
}

process lofreqCall {
    publishDir "$projectDir/code/", mode: 'symlink'

    input:
    file(bam)

    output:
    path "*.vcf.gz"

    script:
    """
    python $projectDir/code/lofreq.py -f ${params.genome} -b ${bam}
    """
}

process freebayesCall {
    publishDir "$projectDir/code/", mode: 'symlink'
    input:
    file(bam)

    output:
    path "fbayes.*.vcf.gz"

    script:
    """
    python $projectDir/code/freebayes.py -f ${params.ugenome} -b ${bam}
    """
}

process mutectCall {
    publishDir "$projectDir/code/", mode: 'symlink'

    input:
    file(bam)

    output:
    path "mutect2*.vcf.gz"

    script:
    """
    python $projectDir/code/mutect2.py -R ${params.genome} -b ${bam}
    """
}

process preprocessBCF {
    publishDir "$projectDir/code/", mode: 'symlink'

    input:
    file(vcf) 
    file (bam)

    output:
    path "gatk*.vcf.gz"

    script:
    """
    python $projectDir/code/annotate_serum.py -f ${params.genome} -v ${vcf} -b ${bam}
    """
}

process preprocessFreebayes {
    publishDir "$projectDir/code/", mode: 'symlink'

    input:
    file(vcf) 
    file (bam)

    output:
    path "gatk*.vcf.gz"

    script:
    """
    python $projectDir/code/annotate_serum.py -f ${params.genome} -v ${vcf} -b ${bam}
    
    """
}


// Step 3: Variant 'filtering'
process annotateVCFs {
    publishDir "$projectDir/code/", mode: 'symlink'

    input:
    file(bcftools)
    file(lofreq)
    file(fbayes)
    file(mutect2)


    output:
    path "COSMIC_*.vcf.gz"

    script:
    """
    $projectDir/code/serum_filter.sh -c ${params.dbsnp} -d ${params.cosmic} -b ${bcftools} -l ${lofreq} -f ${fbayes} -m ${mutect2} -p $projectDir/
    """
}

// Step 4: Merging VCFs
process mergeVCFs {
    input:
    file(x)

    output:
    file ("*dedup.snps.vcf.gz")

    script:
    """
    samtools index $projectDir/code/*bqsr.bam
    python $projectDir/code/merge.py -f ${params.genome}/ -d $projectDir/code/ -x ${x}
    """
}

// Step 5: Dataframe creation
process generateDataframe {
    publishDir "$projectDir/code/", mode: 'symlink'
    
    input:
    file(bcftools)
    file(lofreq)
    file(fbayes)
    file(mutect2)
    file(merged)
    
    
    output:
    file ("dataframe.csv")

    script:
    """    
    python $projectDir/code/dataframe.py -f ${params.genome} -b ${bcftools} -r ${fbayes} -l ${lofreq} -m ${mutect2} -x ${merged}
    """
}

// Step 6: Prediction
process predictVariants {
    publishDir "$projectDir/code/", mode: 'move'
    input:
    file(merged)
    file(dataframe)
    
    output:

    script:
    """    
    python $projectDir/code/predictor.py -o ${params.output} -x ${params.threshold} -d ${merged} -m ${params.model} -q ${dataframe}
    rm -f $projectDir/code/*.bai $projectDir/code/COSMOS_* $projectDir/code/COSMIC_* $projectDir/code/gatkbcftools*vcf* $projectDir/code/gatkfbayes*vcf* $projectDir/code/lofreq.*vcf* $projectDir/code/mutect2*vcf* $projectDir/code/bcftools*vcf* $projectDir/code/*.bqsr.bam* $projectDir/code/fbayes.*vcf* $projectDir/code/*.dedup.snps.vcf.gz* $projectDir/code/*_vcfs.list $projectDir/code/dataframe.csv $projectDir/code/*table $projectDir/code/*rg.bam
    """
}


// Workflow
workflow {
    bam_ch = Channel.fromFilePairs( '*{1,2}.fastq.gz') | bwaMem
    indexVCF = Channel.fromPath("databases/*vcf.gz") | indexVCFs
    gatk_ch = bamPreprocessing(bam_ch)
    bqsr_ch = bqsrProcess(gatk_ch)
    apply_ch = applyBQSR(gatk_ch, bqsr_ch)
    bcf_ch  = bcftoolsCall(apply_ch)
    lof_ch  = lofreqCall(apply_ch)
    bay_ch  = freebayesCall(apply_ch)
    mut_ch  = mutectCall(apply_ch)
    ann1_ch = preprocessBCF(bcf_ch, apply_ch)
    ann2_ch = preprocessFreebayes(bay_ch, apply_ch)
    text_ch = annotateVCFs(ann1_ch, lof_ch, ann2_ch, mut_ch)
    merge_ch = mergeVCFs(text_ch)
    datafr_ch = generateDataframe(ann1_ch, ann2_ch, lof_ch, mut_ch, merge_ch)
    predictVariants(merge_ch, datafr_ch)
}
