""" Sample the repackaged cosmoDC2 catalog from step 00 and
hexgrid the galaxies onto tract 4430.
"""
import os
from lsst.daf.butler import Butler
import numpy as np
from ssi_tools.layout_utils import make_hexgrid_for_tract
import astropy.io.fits as fits
from astropy.table import Table
from tqdm import tqdm

repo = "/global/cfs/cdirs/lsst/production/gen3/DC2/Run2.2i/repo"
collections = "u/descdm/coadds_Y1_4430"
butler = Butler(repo, collections=collections)
skymap = butler.get("skyMap")
tract = skymap[4430]
grid = make_hexgrid_for_tract(tract, rng=57721)

hdul = fits.open(
    os.path.join(
        "/global/cfs/cdirs/lsst/groups/fake-source-injection/ssi_vs_resim/",
        "ssi_catalog.fits"
    )
)
cat = hdul[1].data

# randomly subsample cat with replacement
rng = np.random.default_rng(57721)
indices = rng.choice(len(cat), size=len(grid), replace=True)

out = Table()
out['ra'] = np.deg2rad(grid['ra'])
out['dec'] = np.deg2rad(grid['dec'])
out['x'] = grid['x']
out['y'] = grid['y']
out['original_ra'] = np.deg2rad(cat['ra'][indices])
out['original_dec'] = np.deg2rad(cat['dec'][indices])

for col in tqdm(cat.dtype.names):
    if col in ['ra', 'dec']:
        continue
    out[col] = cat[col][indices]

out.write(
    "/global/cfs/cdirs/lsst/groups/fake-source-injection/ssi_vs_resim/"
    "tract_4430.fits",
    overwrite=True
)
