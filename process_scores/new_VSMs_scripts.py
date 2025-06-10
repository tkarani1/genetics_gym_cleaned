import hail as hl
hl.init()

# ESM1b
esm1b_ht = hl.import_table(
        'gs://missense-scoring/esm1b-650m-brandes/proc/*.txt.bgz',
        types = {
            'aa_pos':hl.tint,
            'aa_ref':hl.tstr,
            'aa_alt':hl.tstr,
            'esm_score':hl.tfloat,
            'uniprot':hl.tstr,
        },
        force=True,
        missing='',
    )
esm1b_ht = esm1b_ht.rename({'uniprot': 'uniprot_isoform'})
esm1b_ht.write('gs://genetics-gym-not-public/Trisha/updated_data_tables/esm1b_scores.ht')

# MisFit
misfit_ht = hl.import_table(
        'gs://missense-scoring/misfit/raw/MisFit_by_Uniprot/*.txt.gz',
        source_file_field='filename',
        types = {
            'Uniprot_position':hl.tint,
            'AA_alt':hl.tstr,
            'MisFit_D':hl.tfloat,
            'MisFit_S':hl.tfloat,
        },
        force=True,
        missing='',
    )
misfit_ht = misfit_ht.annotate(
    uniprot_id = misfit_ht.filename.split('/')[-1].split('\.')[0],
)
misfit_ht = misfit_ht.rename({'Uniprot_position': 'aa_pos', 'AA_alt': 'aa_alt'})
misfit_ht = misfit_ht.drop('filename')

mf_mapping = hl.import_table('gs://missense-scoring/misfit/raw/geneset_s_gene.txt')
mf_mapping = mf_mapping.key_by('UniprotID')
misfit_ht = misfit_ht.annotate(
    enst = mf_mapping[misfit_ht.uniprot_id].TranscriptID, 
    ensg = mf_mapping[misfit_ht.uniprot_id].GeneID, 
    ensp = mf_mapping[misfit_ht.uniprot_id].ProteinID, 
    gene_symbol = mf_mapping[misfit_ht.uniprot_id].Symbol
)
misfit_ht.write('gs://genetics-gym-not-public/Trisha/updated_data_tables/misfit_scores.ht')

# PopEVE
popeve_raw_ht = hl.import_table(
    'gs://missense-scoring/popEVE_ukbb_20250312/*.csv', # 2025 popeve scores
    delimiter=',',
    source_file_field='filename',
    types = {
        'popEVE':hl.tfloat,
        'popped EVE':hl.tfloat,
        'popped ESM-1v':hl.tfloat,
        'EVE':hl.tfloat,
        'ESM-1v':hl.tfloat,
        
    },
    force=True,
    missing='',
)

popeve_ht = popeve_raw_ht.select(
    aa_pos = hl.int(popeve_raw_ht.mutant[1:-1]),
    aa_ref = popeve_raw_ht.mutant[0],
    aa_alt = popeve_raw_ht.mutant[-1],
    popEVE = popeve_raw_ht['popEVE'], 
    popped_EVE = popeve_raw_ht['popped EVE'], 
    popped_ESM_1v = popeve_raw_ht['popped ESM-1v'], 
    EVE = popeve_raw_ht['EVE'], 
    ESM_1v = popeve_raw_ht['ESM-1v'], 
    NP = popeve_raw_ht.filename.split('/')[-1][:-4]
    )
popeve_ht.write('gs://missense-scoring/popeve_2025.ht')
popeve_ht.write('gs://genetics-gym-not-public/Trisha/updated_data_tables/popeve_2025.ht')

# MPC
# original: 
ht = hl.read_table('gs://asc-v17/mpc_gnomad_2.1.1/mpc_grch38_deduped_with_outliers_2024-04-30.ht')
ht = ht.key_by('locus','alleles')
ht.write('gs://missense-scoring/mpc_grch38_deduped_with_outliers_2024-04-30.ht')

mpc_ht = hl.read_table('gs://missense-scoring/mpc_grch38_deduped_with_outliers_2024-04-30.ht')

# RASP
rasp_table = 'gs://nnfc-fdp-konrad-public/RaSP/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2.ht'
ht = hl.read_table(rasp_table)
ht = ht.naive_coalesce(50)
ht = ht.annotate(
    uniprot_id = ht.uniprot,
    uniprot_aa_pos = hl.int(ht.variant[1:-1]),
    uniprot_aa_ref = ht.variant[0],
    uniprot_aa_alt = ht.variant[-1],
    rasp_score = ht.score_ml,
)
ht = ht.key_by(
    'uniprot_id',
    'uniprot_aa_pos',
    'uniprot_aa_ref',
    'uniprot_aa_alt',
)
ht = ht.select(ht.rasp_score)
ht.write('gs://missense-scoring/mutation/rasp_scores.ht')

# AM
am_table = hl.import_table('gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz', 
                           delimiter='\t', force_bgz=True, no_header=True, comment='#')
am_table_fixed = am_table.annotate(
    locus = hl.locus(am_table.f0, hl.int(am_table.f1), reference_genome='GRCh38'),
    alleles = [am_table.f2, am_table.f3], 
    AM_score = hl.float(am_table.f8), 
    enst_orig = am_table.f6, 
    consequence = am_table.f9, 
    genome = am_table.f4
    )
am_table = am_table_fixed.annotate(enst = am_table_fixed.enst_orig.split('\.')[0])
am_table_fixed = am_table_fixed.drop('f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9')
am_table_fixed_k = am_table_fixed.key_by('locus', 'alleles')
am_table_fixed_k.write('gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/AM_fixed.ht')
