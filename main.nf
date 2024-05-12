#!/usr/bin/env nextflow

// aux
params.gt = "$projectDir/Code/lofreq_gt.sh"
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
    python $projectDir/Code/bwa.py -r ${params.hg38} -t ${params.threads} -f1 ${params.f1} -f2 ${params.f2}
    """
}

process indexVCFs{
    publishDir "$projectDir/", mode: 'symlink'

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
    publishDir "$projectDir/", mode: 'symlink'
    input:
    file(bam)

    output:
    path "*.bqsr.bam"

    script :
    """
    python $projectDir/Code/gatk.py -f ${params.genome} -k ${params.known}
    """
}


// Step 2: Variant calling
process bcftoolsCall {
    publishDir "$projectDir/", mode: 'symlink'

    input:
    file(bam)

    output:
    path "bcftools.*.vcf.gz"

    script:
    """
    python $projectDir/Code/bcftools.py -f ${params.genome} -b ${bam}
    """
}

process lofreqCall {
    publishDir "$projectDir/", mode: 'symlink'

    input:
    file(bam)

    output:
    path "*.vcf.gz"

    script:
    """
    python $projectDir/Code/lofreq.py -f ${params.genome} -b ${bam}
    """
}

process freebayesCall {
    publishDir "$projectDir/", mode: 'symlink'
    input:
    file(bam)

    output:
    path "fbayes.*.vcf.gz"

    script:
    """
    python $projectDir/Code/freebayes.py -f ${params.ugenome} -b ${bam}
    """
}

process mutectCall {
    publishDir "$projectDir/", mode: 'symlink'

    input:
    file(bam)

    output:
    path "mutect2*.vcf.gz"

    script:
    """
    python $projectDir/Code/mutect2.py -R ${params.genome} -b ${bam}
    """
}

process preprocessBCF {
    publishDir "$projectDir/", mode: 'symlink'

    input:
    file(vcf) 
    file (bam)

    output:
    path "gatk*.vcf.gz"

    script:
    """
    python $projectDir/Code/annotate_serum.py -f ${params.genome} -v ${vcf} -b ${bam}
    """
}

process preprocessFreebayes {
    publishDir "$projectDir/", mode: 'symlink'

    input:
    file(vcf) 
    file (bam)

    output:
    path "gatk*.vcf.gz"

    script:
    """
    python $projectDir/Code/annotate_serum.py -f ${params.genome} -v ${vcf} -b ${bam}
    
    """
}


// Step 3: Variant 'filtering'
process annotateVCFs {
    publishDir "$projectDir/", mode: 'symlink'

    input:
    file(bcftools)
    file(lofreq)
    file(fbayes)
    file(mutect2)


    output:
    path "COSMIC_*.vcf.gz"

    script:
    """
    $projectDir/Code/serum_filter.sh -c ${params.dbsnp} -d ${params.cosmic} -b ${bcftools} -l ${lofreq} -f ${fbayes} -m ${mutect2} -p $projectDir/
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
    samtools index $projectDir/*bam
    python $projectDir/Code/merge.py -f ${params.genome}/ -d $projectDir/ -x ${x}
    """
}

// Step 5: Dataframe creation
process generateDataframe {
    publishDir "$projectDir/", mode: 'symlink'
    
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
    python $projectDir/Code/dataframe.py -f ${params.genome} -b ${bcftools} -r ${fbayes} -l ${lofreq} -m ${mutect2} -x ${merged}
    """
}

// Step 6: Prediction
process predictVariants {
    publishDir "$projectDir/", mode: 'move'
    input:
    file(merged)
    file(dataframe)
    
    output:

    script:
    """    
    python $projectDir/Code/predictor.py -o ${params.output} -x ${params.threshold} -d ${merged} -m ${params.model} -q ${dataframe}
    rm -f $projectDir/*.bai $projectDir/COSMOS_* $projectDir/COSMIC_* $projectDir/gatkbcftools*vcf* $projectDir/gatkfbayes*vcf* $projectDir/lofreq.*vcf* $projectDir/mutect2*vcf* $projectDir/bcftools*vcf* $projectDir/*.bqsr.bam* $projectDir/fbayes.*vcf* $projectDir/*.dedup.snps.vcf.gz* $projectDir/*_vcfs.list $projectDir/dataframe.csv
    """
}


// Workflow
workflow {
    bam_ch = Channel.fromFilePairs( '*{1,2}.fastq.gz') | bwaMem
    Channel.fromPath("*vcf.gz") | indexVCFs
    gatk_ch = bamPreprocessing(bam_ch)
    bcf_ch  = bcftoolsCall(gatk_ch)
    lof_ch  = lofreqCall(gatk_ch)
    bay_ch  = freebayesCall(gatk_ch)
    mut_ch  = mutectCall(gatk_ch)
    ann1_ch = preprocessBCF(bcf_ch, gatk_ch)
    ann2_ch = preprocessFreebayes(bay_ch, gatk_ch)
    text_ch = annotateVCFs(ann1_ch, lof_ch, ann2_ch, mut_ch)
    merge_ch = mergeVCFs(text_ch)
    datafr_ch = generateDataframe(ann1_ch, ann2_ch, lof_ch, mut_ch, merge_ch)
    predictVariants(merge_ch, datafr_ch)
}
