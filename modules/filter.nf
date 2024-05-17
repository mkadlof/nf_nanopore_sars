process filter {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId),  path('sorted.bam'), path('sorted.bam.bai')
    path(primers)

    output:
    tuple val(sampleId), path('forvariants.bam'), path('forvariants.bam.bai')

    script:
    """
    mask_rubbish=`echo "${params.mask} + 10" | bc -l`
    simple_filter_nanopore.py sorted.bam ${primers} \
                                         ${params.bed_offset} \
                                         ${params.max_depth} \
                                         ${params.length} \
                                         ${params.mapq} \
                                         ${params.extra_bed_offset} \
                                         \${mask_rubbish}
                                         
    # Step 4. Laczenie, sortowanie i indeksowanie plikow ze skryptu wyzej
    samtools sort -@ ${params.threads} -o to_amplicon_clip.bam reads_inner_strict.bam
    samtools index to_amplicon_clip.bam

    ## Laczenie czastkowych plikow bam zaweirajacych ready obejmujace polaczone amplikony
    TO_MERGE_OVER=`ls -l  reads_two_amplicons_*sorted.bam | tr -s " " |  cut -d " " -f9 | tr "\n" " "`
    TO_MERGE_OVER_COUNT=`ls -l  reads_two_amplicons*sorted.bam | wc -l`

    if [ \${TO_MERGE_OVER_COUNT} -gt 0 ];then
            echo "Merging reads from ampicons fusion"
            echo "samtools merge -o tmp_bis.bam \${TO_MERGE_OVER}"
        samtools merge -o tmp_bis.bam \${TO_MERGE_OVER}
        samtools sort -@ ${params.threads} -o two_amplicon_sorted.bam tmp_bis.bam
        samtools index two_amplicon_sorted.bam

        for P in `ls Statystyki_two_amplicons_*`; do if [ -s \${P} ]; then cat \${P} >> Statystyki.txt ; fi ; done
    fi
    
    # Step 5 Maskowanie primerow
    length_filtering=`echo "${params.length} - 40 " | bc -l`
    bed_offset_filtering=`echo "${params.bed_offset} + ${params.extra_bed_offset}" | bc -l`

    ## Maskowanie ze standardowo tolerancja
    samtools ampliconclip --filter-len \${length_filtering} --both-ends -b ${primers} --tolerance ${params.bed_offset} -o clip_part1.bam -O bam to_amplicon_clip.bam

    ## Maskowanie z ekstra tolerancja dla readow overshoot
    if [ -e reads_overshot.bam ]; then
        samtools ampliconclip --filter-len \${length_filtering} --both-ends -b ${primers} --tolerance \${bed_offset_filtering} -o clip_part2.bam -O bam reads_overshot.bam
    fi

    ## Maskowanie z ekstra tolerancja dwa amplikony
    if [ -e two_amplicon_sorted.bam ];then
        samtools ampliconclip --filter-len \${length_filtering} --both-ends -b ${primers} --tolerance \${bed_offset_filtering} -o clip_part3.bam -O bam two_amplicon_sorted.bam
        rm tmp_bis.bam
    fi

    #Maskowanie odczytow obejmujacych jeden amplikon ale tylko jeden z primerow
    if [ -e reads_partial_strict_sort.bam ]; then
        samtools ampliconclip --filter-len \${length_filtering} --both-ends -b ${primers} --tolerance ${params.bed_offset} -o clip_part4.bam -O bam reads_partial_strict_sort.bam
    fi

    #Maskowanie odczytow-smieci z midnight
    if [ -e reads_smieci_sorted.bam ]; then
        samtools ampliconclip --filter-len 1 --both-ends -b ${primers} --tolerance ${params.bed_offset} -o clip_part5.bam -O bam reads_smieci_sorted.bam

    fi

    ## Laczenie plikow z 3 niezalenzych filtrowań
    if [ -e clip_part1.bam ]; then
            PART1="clip_part1.bam"
    else
            PART1=''
    fi

    if [ -e clip_part2.bam ]; then
            PART2="clip_part2.bam"
    else
            PART2=''
    fi

    if [ -e clip_part3.bam ]; then
            PART3="clip_part3.bam"
    else
            PART3=''
    fi

    if [ -e clip_part4.bam ]; then
            PART4="clip_part4.bam"
    else
            PART4=''
    fi

    if [ -e clip_part5.bam ]; then
            PART5="clip_part5.bam"
    else
            PART5=''
    fi

    samtools merge -o tmp.bam \${PART1} \${PART2} \${PART3} \${PART4} \${PART5}
    samtools sort -@ ${params.threads} -o forvariants.bam tmp.bam
    samtools index forvariants.bam
    """
}