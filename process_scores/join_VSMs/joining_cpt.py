import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 

method = 'cpt'

vsm_ens = {
    'esm1b': {'path': 'gs://genetics-gym-not-public/Trisha/updated_data_tables/esm1b_scores.ht', 
              'score_cols': ['esm1b']
             },
    'misfit': {'path': 'gs://genetics-gym-not-public/Trisha/updated_data_tables/misfit_scores.ht', 
               'score_cols': ['MisFit_D', 'MisFit_S']
              },
    'AM': {'path': 'gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/AM_fixed.ht', 
           'score_cols': ['AM']
           }, 
    'RASP': {'path': 'gs://missense-scoring/mutation/rasp_scores.ht', 
             'score_cols': ['rasp_score']
            },
    'MPC': {'path': 'gs://missense-scoring/mpc_grch38_deduped_with_outliers_2024-04-30.ht', 
             'score_cols': ['mpc']
            },
    'ProteinMPNN': {'path': 'gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/protein_mpnn.ht', 
             'score_cols': ['proteinmpnn_llr']
            }, 
    'AM_canon': {'path': 'gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/AM_all_aa_sub_uniprot_canon.ht',
                 'score_cols': ['AM']

    },
    'AM_noncanon': {'path': 'gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/AM_aa_noncanon_enst.ht', 
                    'score_cols': ['AM']

    },
    'PAI3D': {'path': 'gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/PAI3D.ht', 
          'score_cols': ['score_PAI3D']
             }, 
    'revel': {'path': 'gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/revel_enst.ht',
            'score_cols': ['revel']

    },
    'cpt': {'path': 'gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/cpt.ht',
            'score_cols': ['cpt1_score']

    }
}

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



# vsm_ht = hl.read_table('gs://trisha-tmp/make_linker/linker2_temp_union_final.ht')
# vsm_ht = vsm_ht.annotate(
#     uniprot_base = hl.if_else(
#         vsm_ht.uniprot_em_id.contains('-'), vsm_ht.uniprot_em_id.split('-')[0], vsm_ht.uniprot_em_id
#     )
# )

# vsm_ht = vsm_ht.key_by('uniprot_base', 'aa_pos', 'aa_ref', 'aa_alt')

# ht = hl.read_table(vsm_ens[method]['path'])
# # ht = ht.key_by('uniprot_id', 'aa_pos', 'aa_ref', 'aa_alt') # already keyed
# vsm_ht = vsm_ht.join(ht, how='left')
# vsm_ht = vsm_ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/vsm_{method}.ht', overwrite=True)

vsm_ht = hl.read_table(f'gs://trisha-tmp/new_VSM_temp/vsm_{method}.ht')
results = check_score_completenes(vsm_ht, vsm_ens[method]['score_cols'][0])
print('Method: ', method)
print(results)








