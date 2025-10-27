import hail as hl
import sys
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
    'esm1b': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_esm1b.ht', 
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
    'AM': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_am_f.ht', 
            'score_cols': ['AM'],
            'higher_is_less_deleterious' : {'AM': False}
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
    }
}
# vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/vsm_base_em_gene.ht')
# vsm_ht = vsm_ht.key_by('locus', 'alleles', 'ensg')
# vsm_ht = vsm_ht.checkpoint('gs://trisha-tmp/new_VSM_temp/vsm_base_em_gene_k_lag.ht', overwrite=True)
vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/vsm_base_em_gene_k_lag.ht')
print('done')
# For adding ANY 

import sys

if len(sys.argv) < 2:
    print("Usage: python coalesce_vsm.py <method>")
    sys.exit(1)

method = sys.argv[1]
if method not in joined_paths:
    print(f"Method '{method}' not found in joined_paths.")
    sys.exit(1)

info = joined_paths[method]
ht = hl.read_table(info['path'])
print(ht.describe())
ht = ht.key_by('locus', 'alleles', 'ensg')
ht = ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/per_SNP_gene/vsm_k_lag_{method}.ht', overwrite=True)
ht = ht.collect_by_key()
ht = ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/per_SNP_gene/vsm_k_lag_{method}_collected.ht', overwrite=True)
ht.describe()
score_cols = info['score_cols']
higher_is_less_deleterious = info['higher_is_less_deleterious']
final_scores = []
for s in score_cols:
    ht = ht.annotate(
        **{f'{s}_mane_scores':  hl.filter(lambda x: x.mane_select == True, ht.values).map(lambda x: x[s])},
        **{f'{s}_canon_scores':  hl.filter(lambda x: x.canonical == True, ht.values).map(lambda x: x[s])},
        **{f'{s}_any_scores':  hl.filter(lambda x: hl.is_defined(x[s]), ht.values).map(lambda x: x[s])},
    )
    if higher_is_less_deleterious[s]:
        ht = ht.annotate(
            **{f'{s}_mane_scores_neg': ht[f'{s}_mane_scores'].map(lambda x: 1 - x)},
            **{f'{s}_canon_scores_neg': ht[f'{s}_canon_scores'].map(lambda x: 1 - x)},
            **{f'{s}_any_scores_neg': ht[f'{s}_any_scores'].map(lambda x: 1 - x)},
        )
        ht = ht.annotate(
            **{f'{s}_mane_max': hl.max(ht[f'{s}_mane_scores_neg'])},
            **{f'{s}_canon_max': hl.max(ht[f'{s}_canon_scores_neg'])},
            **{f'{s}_any_max': hl.max(ht[f'{s}_any_scores_neg'])},
        )
        ht = ht.annotate(
            **{f'{s}_neg': hl.coalesce(ht[f'{s}_mane_max'], ht[f'{s}_canon_max'], ht[f'{s}_any_max'])}
        )
        to_drop = [f'{s}_mane_scores_neg', f'{s}_canon_scores_neg', f'{s}_any_scores_neg', f'{s}_mane_max', f'{s}_canon_max', f'{s}_any_max']
        ht = ht.drop(*to_drop)
        final_scores.append(f'{s}_neg')
    else: 
        # append
        ht = ht.annotate(
            **{f'{s}_mane_max': hl.max(ht[f'{s}_mane_scores'])},
            **{f'{s}_canon_max': hl.max(ht[f'{s}_canon_scores'])},
            **{f'{s}_any_max': hl.max(ht[f'{s}_any_scores'])},
        )
        ht = ht.annotate(
            **{f'{s}': hl.coalesce(ht[f'{s}_mane_max'], ht[f'{s}_canon_max'], ht[f'{s}_any_max'])}
        )
        to_drop = [f'{s}_mane_max', f'{s}_canon_max', f'{s}_any_max']
        ht = ht.drop(*to_drop)
        final_scores.append(f'{s}')
    print(final_scores)
    
ht = ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/per_SNP_gene/max_per_SNP_gene_{method}.ht', overwrite=True)


