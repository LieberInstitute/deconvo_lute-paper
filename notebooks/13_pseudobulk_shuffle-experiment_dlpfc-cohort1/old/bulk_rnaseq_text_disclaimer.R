#--------------------
# Experiment Manifest
#--------------------
# IMPORTANT: this is your sample identifier manifest
#
# DATA DICTIONARY
# * sample: Brnum, the brain identifier, specific to donor
# * position: DLPFC region, anatomic region within DLPFC region, (mid, ant, post)
# * cell.compartment: cellular compartment, (nuc, cyto, bulk)
# * library.prep: DNA library preparation type (polyA, ribozero)
#
# COLUMN IDENTIFIERS
# * base identifer pattern:
#  [sample.id]-[anatomic.position.id]-[cell.compartment.id]-[library.prep.id]
#
# * primary identifier consolidation
# [BLOCK_ID]-[cell.compartment.id]-[library.prep.id]
#
# * secondary identifier consolidation
# [BLOCK_ID]-[EXPERIMENT_CONDITION]
