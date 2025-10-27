import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 


    # Emily uniprot mapping
em_mapping_path = 'gs://genetics-gym-not-public/Trisha/Emily_uniprot_mapping/'
ht_mis_path = 'gs://genetics-gym-not-public/Trisha/linker_ens_no_uni.ht'


chroms = [str(x) for x in range(1, 23, 1)] + ['X']
print(7)
x = 0
for c in chroms: 
    print('x: ', x)
    # format uniprot mapping file
    fname = f'chr{c}_ensembl_tid_to_uniprot.tsv'
    ht_uni = hl.import_table(em_mapping_path + fname, delimiter='\t', 
                                 types = {'transcript_mane_select':hl.tbool},
                                 force=True,
                                 missing='')
    print('a')
    ht_uni = ht_uni.annotate(
        uniprot_em_id = hl.coalesce(ht_uni.uniprot_isoform, ht_uni.uniprotswissprot, ht_uni.uniprotsptrembl)
    )
    
    print('b')
    ht_uni = ht_uni.key_by('ensembl_transcript_id')
    ht_uni = ht_uni.select('uniprot_em_id', 'transcript_mane_select')
    print('c')
    # filter linker
    ht_mis = hl.read_table(ht_mis_path)
    ht_mis = ht_mis.filter(ht_mis.locus.contig == f'chr{c}')
    ht_mis = ht_mis.key_by('enst')
    print('c')

    # join
    ht_mis_chr = ht_mis.join(ht_uni, how='left')
    print('d')

    ht_mis_chr.write(f'gs://trisha-tmp/make_linker/linker_temp_chr_{c}.ht', overwrite=True)
