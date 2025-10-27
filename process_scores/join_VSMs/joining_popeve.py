import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 

method = 'popeve'

popeve_path = 'gs://missense-scoring/popeve_2025.ht'
popeve_scores = ['popEVE', 'EVE', 'ESM-1v']

def check_score_completenes(ht, score_name): 
    check_ht = ht.key_by('locus', 'alleles')
    check_ht = check_ht.collect_by_key()
    check_ht = check_ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/{score_name}_collected.ht', overwrite=True)
    check_ht = check_ht.annotate(
        score_mane = hl.find(lambda x: x.mane_select == True, check_ht.values)[score_name],
        is_mane = hl.find(lambda x: x.mane_select == True, check_ht.values)['mane_select'],
        score_canon =  hl.find(lambda x: x.canonical == True, check_ht.values)[score_name],
        is_canonical = hl.find(lambda x: x.canonical == True, check_ht.values)['canonical'],
        score_any = hl.find(lambda x: hl.is_defined(x[score_name]) == True, check_ht.values)[score_name]

    )  
    check_htf1 = check_ht.filter(hl.is_defined(check_ht.score_mane))
    check_htf2 = check_ht.filter(hl.is_defined(check_ht.score_canon))
    check_htf3 = check_ht.filter(hl.is_defined(check_ht.score_any))
    mane_count = check_htf1.count()
    canon_count = check_htf2.count()
    any_count = check_htf3.count()
    count_results = {'mane': mane_count, 'canon': canon_count, 'any': any_count}
    return count_results



vsm_ht = hl.read_table('gs://genetics-gym-not-public/Trisha/linker_rs_no_uni.ht')

vsm_ht = vsm_ht.key_by('NP', 'aa_pos', 'aa_ref', 'aa_alt')

ht = hl.read_table(popeve_path)
ht = ht.key_by('NP', 'aa_pos', 'aa_ref', 'aa_alt')
ht = ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/popeve_2025_keyed.ht', overwrite=True)
vsm_ht = vsm_ht.join(ht, how='left')
vsm_ht = vsm_ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/vsm_{method}.ht', overwrite=True)

vsm_ht = hl.read_table(f'gs://trisha-tmp/new_VSM_temp/vsm_{method}.ht')
for s in popeve_scores:
    results = check_score_completenes(vsm_ht, s)
    print('Method: ', method)
    print('Score: ', s)
    print(results)
    print('--------------------------------')








