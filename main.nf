nextflow.enable.dsl=2

include {FASTQC} from './modules/fastqc'
include {TRIM} from './modules/trimmomatic'
include {MULTIQC} from './modules/multiqc'
include {BOWTIE2_BUILD} from './modules/bowtie2_build'
include {BOWTIE2_ALIGN} from './modules/bowtie2_align'
include {SAMTOOLS_SORT} from './modules/samtools_sort'
include {SAMTOOLS_IDX} from './modules/samtools_idx'
include {SAMTOOLS_FLAGSTAT} from './modules/samtools_flagstat'
include {BAMCOVERAGE} from './modules/deeptools_bamcoverage'
include {MULTIBWSUMMARY} from './modules/deeptools_multibwsummary'
include {PLOTCORRELATION} from './modules/deeptools_plotcorrelation'
include {COMPUTEMATRIX} from './modules/deeptools_computematrix'
include {PLOTPROFILE} from './modules/deeptools_plotprofile'
include {CALLPEAKS} from './modules/macs3_callpeak'
include {INTERSECT} from './modules/bedtools_intersect'
include {REMOVE} from './modules/bedtools_remove'
include {ANNOTATE} from './modules/homer_annotatepeaks'
include {FINDPEAKS} from './modules/homer_findpeaks'
include {FIND_MOTIFS_GENOME} from './modules/homer_findmotifsgenome'
//include {TAGDIR} from './modules/homer_maketagdir'

workflow {

    // Week1
    Channel.fromPath(params.samplesheet)
    | splitCsv( header: true )
    | map{row -> tuple(row.name, row.path) }
    | set { read_ch }
    
    FASTQC(read_ch)
    BOWTIE2_BUILD(params.genome)
    TRIM(read_ch, params.adapter_fa)
    
    BOWTIE2_ALIGN(TRIM.out.trimmed, BOWTIE2_BUILD.out.index, BOWTIE2_BUILD.out.name)
    SAMTOOLS_FLAGSTAT(BOWTIE2_ALIGN.out.bam)
    
    TRIM.out.log.concat(FASTQC.out.zip, SAMTOOLS_FLAGSTAT.out.flagstat).collect()
    | set { multiqc_ch }
    
    MULTIQC(multiqc_ch)
    
    SAMTOOLS_SORT(BOWTIE2_ALIGN.out)
    SAMTOOLS_IDX(SAMTOOLS_SORT.out)
    BAMCOVERAGE(SAMTOOLS_IDX.out.index)

    BAMCOVERAGE.out.bigwig.collect{ it[1] }
    | set { bws_ch }
    
    // Week2
    // Perform Multibwsummary and plotcorrelation
    MULTIBWSUMMARY(bws_ch)
    PLOTCORRELATION(MULTIBWSUMMARY.out.multibwsummary, params.cortype)

    // //This outputs two tuples containing the paired files from replicates
    // // The IP and Control are stored in maps for named access during peak calling
    // BOWTIE2_ALIGN.out
    // | map { name, path -> tuple(name.split('_')[1], [(path.baseName.split('_')[0]): path]) }
    // | groupTuple(by: 0)
    // | map { rep, maps -> tuple(rep, maps[0] + maps[1])}
    // | map { rep, samples -> tuple(rep, samples.IP, samples.INPUT)}
    // | set { peakcalling_ch }

    // // Peak calling should work with real files, not enough reads to build model
    // CALLPEAKS(peakcalling_ch, params.macs3_genome)

    // // Make reproducible peaks, filter blacklisted regions
    // CALLPEAKS.out.peaks.collect{ it[1] }
    // | map { files -> tuple('repr_peaks', files[0], files[1]) }
    // | set { intersect_ch }
    Channel.of(params.rep1_peaks, params.rep2_peaks)
    .collect()
    .map { files -> tuple('repr_peaks', files[0], files[1]) }
    .set { intersect_ch }

    INTERSECT(intersect_ch)
    REMOVE(INTERSECT.out.intersect, params.blacklist)

    //Annotate peaks and perform motif finding
    ANNOTATE(REMOVE.out.filtered, params.genome, params.gtf)

// Step: Filter only IP samples from BAMCOVERAGE results
BAMCOVERAGE.out.bigwig
| filter { name, file -> name.toLowerCase().contains('ip') }
| map { name, bw -> tuple(name, bw) }
| set { ip_bw_ch }

// Step: Compute matrix with scale-regions, 2kb upstream/downstream
ip_bw_ch
| map { meta, bw -> tuple(meta, bw, file(params.genes_bed), 2000) }
| set { compute_matrix_input_ch }

COMPUTEMATRIX(compute_matrix_input_ch)

// Step: Plot profile based on matrix output
COMPUTEMATRIX.out.matrix
| PLOTPROFILE

    FIND_MOTIFS_GENOME(REMOVE.out.filtered, params.genome)

}
