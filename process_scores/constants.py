"""Constants for processing missense scores."""

########################################################################################
# Top-level fields
########################################################################################
ENSEMBL_TRANSCRIPT_ID_FIELD = "enst"
"""The Ensembl transcript ID field."""
ENSEMBL_GENE_ID_FIELD = "ensg"
"""The Ensembl gene ID field."""
ENSEMBL_PROTEIN_ID_FIELD = "ensp"
"""The Ensembl protein ID field."""
GENE_SYMBOL_FIELD = "gene_symbol"
"""The gene symbol field."""
AA_POS_FIELD = "aa_pos"
"""The amino acid position field."""
AA_REF_FIELD = "aa_ref"
"""The amino acid reference field."""
AA_ALT_FIELD = "aa_alt"
"""The amino acid alternative field."""
UNIPROT_ID_FIELD = "uniprot_id"
"""The UniProt ID field."""
UNIPROT_ISOFORM_FIELD = "uniprot_isoform"
"""The UniProt isoform field."""

########################################################################################
# Key field groups
########################################################################################
VARIANT_KEY_FIELDS = ["locus", "alleles"]
"""The variant key fields."""
VARIANT_TRANSCRIPT_KEY_FIELDS = [*VARIANT_KEY_FIELDS, ENSEMBL_TRANSCRIPT_ID_FIELD]
"""The variant transcript key fields."""
AA_SUBSTITUTION_KEY_FIELDS = [AA_POS_FIELD, AA_REF_FIELD, AA_ALT_FIELD]
"""The amino acid substitution key fields."""
AA_SUBSTITUTION_UNIPROT_KEY_FIELDS = [*AA_SUBSTITUTION_KEY_FIELDS, UNIPROT_ID_FIELD]
"""The amino acid substitution UniProt key fields."""
AA_SUBSTITUTION_TRANSCRIPT_KEY_FIELDS = [
    *AA_SUBSTITUTION_KEY_FIELDS,
    ENSEMBL_TRANSCRIPT_ID_FIELD,
]
"""The amino acid substitution UniProt transcript key fields."""
VARIANT_UNIPROT_TRANSCRIPT_KEY_FIELDS = [
    *VARIANT_KEY_FIELDS,
    UNIPROT_ID_FIELD,
    ENSEMBL_TRANSCRIPT_ID_FIELD,
]
"""The variant UniProt transcript key fields."""
AA_ALT_UNIPROT_TRANSCRIPT_KEY_FIELDS = [
    AA_POS_FIELD,
    AA_ALT_FIELD,
    UNIPROT_ID_FIELD,
    ENSEMBL_TRANSCRIPT_ID_FIELD,
]
"""The amino acid substitution alt UniProt and transcript key fields."""
LINKER_KEY_FIELDS = [
    *VARIANT_KEY_FIELDS,
    *AA_SUBSTITUTION_UNIPROT_KEY_FIELDS,
    ENSEMBL_TRANSCRIPT_ID_FIELD,
]
"""The linker key fields."""

KEY_GROUPS = {
    "variant": VARIANT_KEY_FIELDS,
    "variant_transcript": VARIANT_TRANSCRIPT_KEY_FIELDS,
    "variant_uniprot_transcript": VARIANT_UNIPROT_TRANSCRIPT_KEY_FIELDS,
    "aa_substitution": AA_SUBSTITUTION_KEY_FIELDS,
    "aa_substitution_transcript": AA_SUBSTITUTION_TRANSCRIPT_KEY_FIELDS,
    "aa_substitution_uniprot": AA_SUBSTITUTION_UNIPROT_KEY_FIELDS,
    "aa_substitution_alt_uniprot_transcript": AA_ALT_UNIPROT_TRANSCRIPT_KEY_FIELDS,
}
"""The table key groups."""

########################################################################################
# Score fields
########################################################################################

ESM1B_SCORE_FIELD = "esm1b"
"""The ESM1B score field."""
MISFIT_D_SCORE_FIELD = "MisFit_D"
"""The MisFit D score field."""
MISFIT_S_SCORE_FIELD = "MisFit_S"
"""The MisFit S score field."""
POPEVE_SCORE_FIELD = "popEVE"
"""The popEVE score field."""
EVE_SCORE_FIELD = "EVE"
"""The EVE score field."""
ESM_1V_SCORE_FIELD = "ESM_1v"
"""The ESM-1v score field."""
MPC_SCORE_FIELD = "mpc"
"""The MPC score field."""
RASP_SCORE_FIELD = "rasp_score"
"""The RaSP score field."""
REVEL_SCORE_FIELD = "revel"
"""The REVEL score field."""
CPT1_SCORE_FIELD = "cpt1_score"
"""The CPT1 score field."""
PROTEINMPNN_LLR_SCORE_FIELD = "proteinmpnn_llr"
"""The ProteinMPNN LLR score field."""
PAI3D_SCORE_FIELD = "score_PAI3D"
"""The PAI3D score field."""
POLYPHEN_SCORE_FIELD = "polyphen_score"
"""The PolyPhen score field."""
CADD_SCORE_FIELD = "cadd_score"
"""The CADD score field."""
GPN_MSA_SCORE_FIELD = "gpn_msa_score"
"""The GPN-MSA score field."""
AM_SCORE_FIELD = "AM_score"
"""The AlphaMissense score field."""
AM_CANONICAL_SCORE_FIELD = "AM_score"
"""The AlphaMissense canonical score field."""

SCORE_FIELDS = {
    "esm1b": [ESM1B_SCORE_FIELD],
    "misfit": [MISFIT_D_SCORE_FIELD, MISFIT_S_SCORE_FIELD],
    "popeve": [POPEVE_SCORE_FIELD, EVE_SCORE_FIELD, ESM_1V_SCORE_FIELD],
    "mpc": [MPC_SCORE_FIELD],
    "rasp": [RASP_SCORE_FIELD],
    "revel": [REVEL_SCORE_FIELD],
    "cpt": [CPT1_SCORE_FIELD],
    "proteinmpnn": [PROTEINMPNN_LLR_SCORE_FIELD],
    "pai3d": [PAI3D_SCORE_FIELD],
    "polyphen": [POLYPHEN_SCORE_FIELD],
    "cadd": [CADD_SCORE_FIELD],
    "gpn_msa": [GPN_MSA_SCORE_FIELD],
    "am_isos": [AM_SCORE_FIELD],
    "am_canonical": [AM_SCORE_FIELD],
}
"""The score fields in each score resource."""

SCORE_KEY_GROUPS = {
    "esm1b": "aa_substitution_uniprot",
    "misfit": "aa_substitution_alt_uniprot_transcript",
    "popeve": "aa_substitution",
    # "mpc": "variant", # TODO: Uncomment this when MPC scores are processed.
    "rasp": "aa_substitution_uniprot",
    "revel": "variant_transcript",
    "cpt": "aa_substitution",
    "proteinmpnn": "aa_substitution_uniprot",
    "pai3d": "variant_transcript",
    "polyphen": "variant_transcript",
    "cadd": "variant",
    "gpn_msa": "variant",
    "am_isos": "aa_substitution_transcript",
    "am_canonical": "variant_uniprot_transcript",
}

HIGHER_IS_LESS_DELETERIOUS = {
    "esm1b": {ESM1B_SCORE_FIELD: True},
    "misfit": {MISFIT_D_SCORE_FIELD: False, MISFIT_S_SCORE_FIELD: False},
    "popeve": {
        POPEVE_SCORE_FIELD: True,
        EVE_SCORE_FIELD: False,
        ESM_1V_SCORE_FIELD: True,
    },
    "mpc": {MPC_SCORE_FIELD: False},
    "rasp": {RASP_SCORE_FIELD: False},
    "revel": {REVEL_SCORE_FIELD: False},
    "cpt": {CPT1_SCORE_FIELD: False},
    "proteinmpnn": {PROTEINMPNN_LLR_SCORE_FIELD: True},
    "pai3d": {PAI3D_SCORE_FIELD: False},
    "polyphen": {POLYPHEN_SCORE_FIELD: False},
    "cadd": {CADD_SCORE_FIELD: False},
    "gpn_msa": {GPN_MSA_SCORE_FIELD: False},
    "am_isos": {AM_SCORE_FIELD: False},
    "am_canonical": {AM_SCORE_FIELD: False},
}
"""Whether higher scores are less deleterious for each score."""
