import numpy as np
from desc_dc2_dm_data import REPOS
from lsst.daf.persistence import Butler
import lsst.geom as geom
from astropy.table import vstack
from tqdm.notebook import tqdm


repo = REPOS['2.2i_dr6_wfd']
butler = Butler(repo)

skymap = butler.get("deepCoadd_skyMap")
tractInfo = skymap[4030]

ccds = []
for ipatch in tqdm(range(49)):
    patchInfo = tractInfo.getPatchInfo(ipatch)
    for filt in 'ugrizy':
        bbox = butler.get(
            "deepCoadd_bbox",
            tract=tractInfo.getId(),
            patch="%d,%d"%patchInfo.getIndex(),
            filter=filt
        )
        bbox = geom.Box2I(bbox.getMin(), bbox.getMin()+geom.Extent2I(1, 1))
        coadd = butler.get(
            "deepCoadd_sub",
            bbox=bbox,
            tract=tractInfo.getId(),
            patch="%d,%d"%patchInfo.getIndex(),
            filter=filt
        )
        ccds.append(coadd.getInfo().getCoaddInputs().ccds.asAstropy())

table = vstack(ccds)

print(table)
