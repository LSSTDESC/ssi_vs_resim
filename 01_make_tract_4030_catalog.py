""" Sample the repackaged cosmoDC2 catalog from step 00 and
hexgrid the galaxies onto tract 4030.
"""
import os
from desc_dc2_dm_data import REPOS
from lsst.daf.persistence import Butler
import numpy as np
from ssi_tools.layout_utils import make_hexgrid_for_tract
import astropy.io.fits as fits
from astropy.table import Table

butler = Butler(REPOS['2.2i_dr6_wfd'])
skymap = butler.get("deepCoadd_skyMap")
tract = skymap[4030]
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

for col in cat.dtype.names:
    if col in ['ra', 'dec']:
        continue
    out[col] = cat[col][indices]

out.write(
    "/global/cfs/cdirs/lsst/groups/fake-source-injection/ssi_vs_resim/"
    "tract_4030.fits",
    overwrite=True
)
