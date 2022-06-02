"""This script repackages a subset of the DC2 data into a file that
corresponds to 4 square degrees of objects brighter than i = 28.

The resulting file is ~500 MB which is a useful size in that it is not too big
but not too small either.
"""
import os
import tqdm
import numpy as np
import fitsio

import GCRCatalogs

# numpy structured array dtype for the outputs
dtype = [
    ("ra", "f8"),
    ("dec", "f8"),
    ("redshift", "f4"),
    ("u_mag", "f4"),
    ("g_mag", "f4"),
    ("r_mag", "f4"),
    ("i_mag", "f4"),
    ("z_mag", "f4"),
    ("y_mag", "f4"),
    ("sourceType", "U6"),
    ("bulge_semimajor", "f4"),
    ("bulge_axis_ratio", "f4"),
    ("bulge_pa", "f4"),
    ("bulge_n", "f4"),
    ("disk_semimajor", "f4"),
    ("disk_axis_ratio", "f4"),
    ("disk_pa", "f4"),
    ("disk_n", "f4"),
    ("u_bulge_disk_flux_ratio", "f4"),
    ("g_bulge_disk_flux_ratio", "f4"),
    ("r_bulge_disk_flux_ratio", "f4"),
    ("i_bulge_disk_flux_ratio", "f4"),
    ("z_bulge_disk_flux_ratio", "f4"),
    ("y_bulge_disk_flux_ratio", "f4"),
    ("bulge_disk_flux_ratio", "f4"),
]

# these are mappings from the column names expected by the ssi code and
# the DC2 column names.
dc2_name_map = {
    "ra": "ra",
    "dec": "dec",
    "redshift": "redshift",
    "u_mag": "mag_true_u_lsst",
    "g_mag": "mag_true_g_lsst",
    "r_mag": "mag_true_r_lsst",
    "i_mag": "mag_true_i_lsst",
    "z_mag": "mag_true_z_lsst",
    "y_mag": "mag_true_Y_lsst",
    "bulge_semimajor": "morphology/spheroidMajorAxisArcsec",
    "bulge_axis_ratio": "morphology/spheroidAxisRatio",
    "bulge_pa": "morphology/positionAngle",
    "bulge_n": "morphology/spheroidSersicIndex",
    "disk_semimajor": "morphology/diskMajorAxisArcsec",
    "disk_axis_ratio": "morphology/diskAxisRatio",
    "disk_pa": "morphology/positionAngle",
    "disk_n": "morphology/diskSersicIndex",
}

cat = GCRCatalogs.load_catalog("cosmoDC2_v1.1.4_small")

# we keep a catalog corresponding to objects from 4 square degrees
# of everything brighter than i = 28
rng = np.random.RandomState(seed=3463)
frac = 4.0/cat.sky_area

# we have to make the random selection after the i mag cut, so we do that column
# first
col_d = cat.get_quantities("mag_true_i_lsst")["mag_true_i_lsst"]
imsk = (col_d <= 28.0)
tot = int(np.sum(imsk))
num = int(tot * frac)
inds = rng.choice(tot, size=num, replace=False)
d = np.zeros(num, dtype=dtype)

for col, dc2_col in tqdm.tqdm(dc2_name_map.items()):
    d[col] = cat.get_quantities(dc2_col)[dc2_col][imsk][inds]
d["sourceType"] = "galaxy"

# Fill in b/d flux ratios manually
cols = []
for band in tqdm.tqdm('ugrizy'):
    bname = f"LSST_filters/spheroidLuminositiesStellar:LSST_{band}:observed"
    dname = f"LSST_filters/diskLuminositiesStellar:LSST_{band}:observed"
    bflux = cat.get_quantities(bname)[bname][imsk][inds]
    dflux = cat.get_quantities(dname)[dname][imsk][inds]
    d[f"{band}_bulge_disk_flux_ratio"] = bflux/dflux

# Add a placeholder "bulge_disk_flux_ratio" column filled with nans so we don't accidentally
# use the default bdfr of 1.0
d["bulge_disk_flux_ratio"] = np.nan

os.makedirs(
    "/global/cfs/cdirs/lsst/groups/fake-source-injection/ssi_vs_resim/",
    exist_ok=True,
)
fitsio.write(
    "/global/cfs/cdirs/lsst/groups/fake-source-injection/ssi_vs_resim/"
    "ssi_catalog.fits",
    d,
    clobber=True,
)
