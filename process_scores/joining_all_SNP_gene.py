import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 

# all are keyed by locus, alleles
collected_paths = {
    'proteinmpnn': {'path': 'gs://trisha-tmp/new_VSM_temp/proteinmpnn_llr_collected.ht', 
             'score_cols': ['proteinmpnn_llr'], 
             'higher_is_less_deleterious' : {'proteinmpnn_llr': True}
            }, 
    'esm1b': {'path': 'gs://trisha-tmp/new_VSM_temp/esm1b_lj2_collected.ht', 
              'score_cols': ['esm_comb'],
              'higher_is_less_deleterious' : {'esm_comb': True}
             },
    'pai3d': {'path': 'gs://trisha-tmp/new_VSM_temp/score_PAI3D_collected.ht', 
              'score_cols': ['score_PAI3D'], 
              'higher_is_less_deleterious' : {'score_PAI3D': False}
            },
    'revel': {'path': 'gs://trisha-tmp/new_VSM_temp/revel_collected.ht',
            'score_cols': ['revel'], 
            'higher_is_less_deleterious' : {'revel': False}
            },
    'rasp': {'path': 'gs://trisha-tmp/new_VSM_temp/rasp_score_collected.ht', 
             'score_cols': ['rasp_score'], 
             'higher_is_less_deleterious' : {'rasp_score': False}
            },
    'AM': {'path': 'gs://trisha-tmp/new_VSM_temp/AM_comb_collected.ht', 
            'score_cols': ['AM_comb'],
            'higher_is_less_deleterious' : {'AM_comb': False}
            },
    'misfit_D': {'path': 'gs://trisha-tmp/new_VSM_temp/MisFit_D_collected.ht', 
            'score_cols': ['MisFit_D'],
            'higher_is_less_deleterious' : {'MisFit_D': False}
            },
    'misfit_S': {'path': 'gs://trisha-tmp/new_VSM_temp/MisFit_S_collected.ht', 
            'score_cols': ['MisFit_S'],
            'higher_is_less_deleterious' : {'MisFit_S': False}
            },
    'mpc': {'path': 'gs://missense-scoring/mpc_grch38_deduped_with_outliers_2024-04-30.ht',
            'score_cols': ['mpc'],
            'higher_is_less_deleterious' : {'mpc': False}
        }
}

joined_paths = {
    'proteinmpnn': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_proteinmpnn.ht', 
             'score_cols': ['proteinmpnn_llr'], 
             'higher_is_less_deleterious' : {'proteinmpnn_llr': True}
            }, 
    'esm1b': {'path': 'gs://trisha-tmp/new_VSM_temp/esm1b_lj2.ht', 
              'score_cols': ['esm1b'],
              'higher_is_less_deleterious' : {'esm1b': True}
             },
    'pai3d': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_pai3d.ht', 
              'score_cols': ['score_PAI3D'], 
              'higher_is_less_deleterious' : {'score_PAI3D': False}
            },
    'revel': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_revel.ht',
            'score_cols': ['revel'], 
            'higher_is_less_deleterious' : {'revel': False}
            },
    'rasp': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_rasp.ht', 
             'score_cols': ['rasp_score'], 
             'higher_is_less_deleterious' : {'rasp_score': False}
            },
    'AM': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_am.ht', 
            'score_cols': ['AM_score'],
            'higher_is_less_deleterious' : {'AM_score': False}
            },
    'misfit': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_misfit.ht', 
            'score_cols': ['MisFit_D', 'MisFit_S'],
            'higher_is_less_deleterious' : {'MisFit_D': False, 'MisFit_S': False}
            },
    'popeve': {
        'path': 'gs://trisha-tmp/new_VSM_temp/vsm_popeve.ht',
        'score_cols': ['popEVE', 'EVE', 'ESM_1v'],
        'higher_is_less_deleterious' : {'popEVE': True, 'EVE': False, 'ESM_1v': True}
    },
}
joined2_paths = {
        'polyphen': {
        'path': 'gs://trisha-tmp/new_VSM_temp/vsm_polyphen.ht',
        'score_cols': ['polyphen_score'],
        'higher_is_less_deleterious' : {'polyphen_score': False}
    },
    'cpt': {
        'path': 'gs://trisha-tmp/new_VSM_temp/vsm_cpt.ht',
        'score_cols': ['cpt1_score'],
        'higher_is_less_deleterious' : {'cpt1_score': False}
    }

}

no_transcript_tables = { 
    'mpc': {'path': 'gs://missense-scoring/mpc_grch38_deduped_with_outliers_2024-04-30.ht',
            'score_cols': ['mpc'],
            'higher_is_less_deleterious' : {'mpc': False}
        },
    'cadd': {
        'path': 'gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/cadd.ht',
        'score_cols': ['cadd_score'],
        'higher_is_less_deleterious' : {'cadd_score': False}
    },
    'gpn_msa': {
        'path': 'gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/gpn_msa.ht', 
        'score_cols': ['gpn_msa_score'],
        'higher_is_less_deleterious' : {'gpn_msa_score': False}
    }
}
    # just need to run this the first time
# vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/vsm_base_em_gene.ht')
# vsm_ht = vsm_ht.key_by('locus', 'alleles', 'ensg')
# vsm_ht = vsm_ht.checkpoint('gs://trisha-tmp/new_VSM_temp/vsm_base_em_gene_k_lag.ht', overwrite=True)
# vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/vsm_base_em_gene_k_lag.ht')

# For adding ANY 
# all_score_cols = []
# x = 0
# x = 6
# vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_misfit_5.ht')
# for method, info in joined2_paths.items(): 
#     print(method)
#     # if method == 'popeve': 
#     #     continue
#     ht = hl.read_table(f'gs://trisha-tmp/new_VSM_temp/per_SNP_gene/max_per_SNP_gene_{method}.ht')
#     final_score_names = []
#     for s in info['score_cols']:
#         if info['higher_is_less_deleterious'][s]:
#             final_score_names.append(f'{s}_neg')
#         else:
#             final_score_names.append(s)
#     ht = ht.select(*final_score_names)
#     vsm_ht = vsm_ht.join(ht, how='left')
#     vsm_ht = vsm_ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_{method}_{x}.ht', overwrite=True)
#     x += 1
#     vsm_ht.show()

# x = 8
# method = 'popeve'
# info = joined_paths[method]
# vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_cpt_7.ht')
# vsm_ht = vsm_ht.key_by('locus', 'alleles', 'gene_symbol')
# ht = hl.read_table(f'gs://trisha-tmp/new_VSM_temp/per_SNP_gene/max_per_SNP_gene_{method}.ht')
# final_score_names = []
# for s in info['score_cols']:
#     if info['higher_is_less_deleterious'][s]:
#         final_score_names.append(f'{s}_neg')
#     else:
#         final_score_names.append(s)
# ht = ht.select(*final_score_names)
# vsm_ht = vsm_ht.join(ht, how='left')
# vsm_ht = vsm_ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_{method}_{x}.ht', overwrite=True)
# x += 1
# vsm_ht.show()

# method = 'popeve'
# x = 8
# vsm_ht = hl.read_table(f'gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_{method}_{x}.ht')
# vsm_ht = vsm_ht.key_by('locus', 'alleles')

# for method, info in no_transcript_tables.items(): 
#     print(method)
#     final_score_names = []
#     for s in info['score_cols']:
#         if info['higher_is_less_deleterious'][s]:
#             final_score_names.append(f'{s}_neg')
#         else:
#             final_score_names.append(s)
#     ht = hl.read_table(info['path'])
#     ht = ht.select(*final_score_names)
#     vsm_ht = vsm_ht.join(ht, how='left')
#     vsm_ht = vsm_ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_{method}_{x}.ht', overwrite=True)
#     x += 1
#     vsm_ht.show()

# method = 'gpn_msa'
# vsm_ht = hl.read_table(f'gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_cadd_9.ht')
# x = 10
# for method, info in no_transcript_tables.items(): 
#     if method != 'gpn_msa':
#         continue
#     print(method)
#     final_score_names = []
#     for s in info['score_cols']:
#         if info['higher_is_less_deleterious'][s]:
#             final_score_names.append(f'{s}_neg')
#         else:
#             final_score_names.append(s)
#     ht = hl.read_table(info['path'])
#     ht = ht.select(*final_score_names)
#     vsm_ht = vsm_ht.join(ht, how='left')
#     vsm_ht = vsm_ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_{method}_{x}.ht', overwrite=True)
#     x += 1
#     vsm_ht.show()

# Format 
vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_gpn_msa_10.ht')
vsm_ht = vsm_ht.key_by('locus', 'alleles', 'ensg')
vsm_ht = vsm_ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_all_k_lag.ht', overwrite=True)
vsm_ht = vsm_ht.distinct()
vsm_ht = vsm_ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_all_k_lag_distinct.ht', overwrite=True)

# fix AM 
vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_all_k_lag_distinct.ht')
method = 'AM'
vsm_ht = vsm_ht.drop('AM_score')
print(method)
ht = hl.read_table(f'gs://trisha-tmp/new_VSM_temp/per_SNP_gene/max_per_SNP_gene_{method}_f.ht')
ht = ht.select('AM')
vsm_ht = vsm_ht.join(ht, how='left')
vsm_ht = vsm_ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_all_k_lag.ht', overwrite=True)

# vsm_ht.show()

# Format
vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_all_k_lag.ht')
score_columns = ['proteinmpnn_llr_neg', 'esm1b_neg', 
            'score_PAI3D', 'revel', 'rasp_score', 'AM', 
            'MisFit_D', 'MisFit_S', 'polyphen_score', 'cpt1_score', 
            'popEVE_neg', 'EVE', 'ESM_1v_neg', 'cadd_score', 'gpn_msa_score']
other_columns = ['aa_pos', 'aa_ref', 'aa_alt', 'gene_symbol',
                'ensp','uniprot_base', 'em_gene']
columns_to_specify = ['enst', 'mane_NM', 'mane_select', 'canonical', 'uniprot_em_id',  'transcript_mane_select']
vsm_ht = vsm_ht.collect_by_key()              
for s in score_columns:
    vsm_ht = vsm_ht.annotate(
        **{f'{s}':  hl.find(lambda x: hl.is_defined(x[s]) == True, vsm_ht.values)[s]}
    )
for c in other_columns:
    vsm_ht = vsm_ht.annotate(
        **{f'{c}':  hl.find(lambda x: hl.is_defined(x[c]) == True, vsm_ht.values)[c]}
    )
vsm_ht = vsm_ht.annotate(
    mane_enst = hl.find(lambda x: x.mane_select == True, vsm_ht.values)['enst'], 
    canon_enst = hl.find(lambda x: x.canonical == True, vsm_ht.values)['enst'], 
    **{'all_enst':  hl.filter(lambda x: hl.is_defined(x.enst) == True, vsm_ht.values).map(lambda x: x['enst'])},

    mane_NM = hl.find(lambda x: x.mane_select == True, vsm_ht.values)['mane_NM'], 
    mane_select = hl.if_else(hl.is_defined(hl.find(lambda x: x.mane_select == True, vsm_ht.values)['mane_select']), True, False),
    canonical = hl.if_else(hl.is_defined(hl.find(lambda x: x.canonical == True, vsm_ht.values)['canonical']), True, False),

    mane_uniprot_em_id = hl.find(lambda x: x.mane_select == True, vsm_ht.values)['uniprot_em_id'],
    canon_uniprot_em_id = hl.find(lambda x: x.canonical == True, vsm_ht.values)['uniprot_em_id'],
    **{'all_uniprot_em_id':  hl.filter(lambda x: hl.is_defined(x.uniprot_em_id) == True, vsm_ht.values).map(lambda x: x['uniprot_em_id'])},

    transcript_mane_select = hl.find(lambda x: x.transcript_mane_select == True, vsm_ht.values)['transcript_mane_select'],

)
vsm_ht = vsm_ht.checkpoint('gs://genetics-gym-not-public/Trisha/vsm_all_tables/vsm_all_SNP_gene.ht', overwrite=True)

vsm_ht = vsm_ht.drop('values')
vsm_ht.write('gs://genetics-gym-not-public/Trisha/vsm_all_tables/vsm_all_SNP_gene_cleaned.ht', overwrite=True)
