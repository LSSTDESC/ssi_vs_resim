import numpy as np
from astropy.io import fits
import pandas as pd

data = fits.getdata(
    "/global/cfs/cdirs/lsst/groups/fake-source-injection/ssi_vs_resim/"
    "tract_4430.fits"
)

data = {name: data[name].tolist() for name in data.dtype.names}  # byteswapping here
df = pd.DataFrame(data=data)
df.to_parquet(
    "/global/cfs/cdirs/lsst/groups/fake-source-injection/ssi_vs_resim/"
    "tract_4430.parq"
)
