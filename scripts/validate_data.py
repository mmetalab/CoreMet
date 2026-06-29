#!/usr/bin/env python3
"""Validate data sources, PMID coverage, and evidence quality across all CoreMet databases."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd

from app.config import Config
cfg = Config()

# MPI
mpi = pd.read_csv(cfg.MPI_DB_PATH, dtype=str).fillna('')
print('=== MPI ===')
print(f'  Rows: {len(mpi)}')
src_cols = [c for c in mpi.columns if 'source' in c.lower() or 'pmid' in c.lower() or 'evidence' in c.lower()]
print(f'  Source cols: {src_cols}')
for c in src_cols:
    non_empty = (mpi[c] != '').sum()
    print(f'    {c}: {non_empty}/{len(mpi)} ({100*non_empty/len(mpi):.1f}%)')

# MDI
from app.services.mdi_service import get_mdi_db
mdi = get_mdi_db()
print('\n=== MDI ===')
print(f'  Rows: {len(mdi)}')
src_cols = [c for c in mdi.columns if 'source' in c.lower() or 'pmid' in c.lower() or 'evidence' in c.lower()]
print(f'  Source cols: {src_cols}')
for c in src_cols:
    vals = mdi[c].astype(str)
    non_empty = ((vals != '') & (vals != 'nan')).sum()
    print(f'    {c}: {non_empty}/{len(mdi)} ({100*non_empty/len(mdi):.1f}%)')
if 'Source' in mdi.columns:
    print('  Source dist:')
    for v, cnt in mdi['Source'].value_counts().head(10).items():
        print(f'    {v}: {cnt}')

# MMI
from app.services.mmi_service import get_mmi_db
mmi = get_mmi_db()
print('\n=== MMI ===')
print(f'  Rows: {len(mmi)}')
src_cols = [c for c in mmi.columns if 'source' in c.lower() or 'pmid' in c.lower() or 'evidence' in c.lower()]
print(f'  Source cols: {src_cols}')
for c in src_cols:
    vals = mmi[c].astype(str)
    non_empty = ((vals != '') & (vals != 'nan')).sum()
    print(f'    {c}: {non_empty}/{len(mmi)} ({100*non_empty/len(mmi):.1f}%)')

# MDrI
from app.services.mdri_service import get_mdri_db
mdri = get_mdri_db()
print('\n=== MDrI ===')
print(f'  Rows: {len(mdri)}')
src_cols = [c for c in mdri.columns if 'source' in c.lower() or 'pmid' in c.lower() or 'evidence' in c.lower()]
print(f'  Source cols: {src_cols}')
for c in src_cols:
    vals = mdri[c].astype(str)
    non_empty = ((vals != '') & (vals != 'nan')).sum()
    print(f'    {c}: {non_empty}/{len(mdri)} ({100*non_empty/len(mdri):.1f}%)')

# MGI
from app.services.mgi_service import get_mgi_db
mgi = get_mgi_db()
print('\n=== MGI ===')
print(f'  Rows: {len(mgi)}')
src_cols = [c for c in mgi.columns if 'source' in c.lower() or 'pmid' in c.lower() or 'evidence' in c.lower()]
print(f'  Source cols: {src_cols}')
for c in src_cols:
    vals = mgi[c].astype(str)
    non_empty = ((vals != '') & (vals != 'nan')).sum()
    print(f'    {c}: {non_empty}/{len(mgi)} ({100*non_empty/len(mgi):.1f}%)')

# mGWAS
from app.services.mgwas_service import get_mgwas_db
mgwas = get_mgwas_db()
print('\n=== mGWAS ===')
print(f'  Rows: {len(mgwas)}')
src_cols = [c for c in mgwas.columns if 'source' in c.lower() or 'pmid' in c.lower() or 'evidence' in c.lower()]
print(f'  Source cols: {src_cols}')
for c in src_cols:
    vals = mgwas[c].astype(str)
    non_empty = ((vals != '') & (vals != 'nan')).sum()
    print(f'    {c}: {non_empty}/{len(mgwas)} ({100*non_empty/len(mgwas):.1f}%)')
if 'Source' in mgwas.columns:
    print('  Source dist:')
    for v, cnt in mgwas['Source'].value_counts().items():
        print(f'    {v}: {cnt}')

# Summary
print('\n=== DATA QUALITY SUMMARY ===')
total = len(mpi) + len(mdi) + len(mmi) + len(mdri) + len(mgi) + len(mgwas)
print(f'Total edges: {total:,}')

# PMID coverage
pmid_total = 0
pmid_has = 0
for name, df in [('MPI', mpi), ('MDI', mdi), ('MMI', mmi), ('MDrI', mdri), ('MGI', mgi), ('mGWAS', mgwas)]:
    pc = [c for c in df.columns if 'pmid' in c.lower()]
    pmid_total += len(df)
    if pc:
        vals = df[pc[0]].astype(str)
        has = ((vals != '') & (vals != 'nan') & (vals != '0')).sum()
        pmid_has += has
        print(f'  {name}: PMID {has}/{len(df)} ({100*has/len(df):.1f}%)')
    else:
        print(f'  {name}: No PMID column')

print(f'\nOverall PMID coverage: {pmid_has}/{pmid_total} ({100*pmid_has/pmid_total:.1f}%)')
