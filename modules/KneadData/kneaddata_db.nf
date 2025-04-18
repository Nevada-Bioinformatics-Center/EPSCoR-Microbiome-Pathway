
//process kneaddata_db_human_genome {

//    publishDir "${params.output}/kneaddata_db", mode: 'copy'

//    output:
//        path ("human_genome/*"), emit: human_genome_db

//    script:
//    """
//    kneaddata_database \
//    --download human_genome \
//    --database-location ${params.output}/kneaddata_db/human_genome
//    """
//}