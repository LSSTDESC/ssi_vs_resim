""" Read the tract 4030 SSI catalog and create a phosim-style
instance catalog from it.  In particular, need to 'deprecess' the
ra/dec coords into the phosim celestial reference frame.  Also
handle ellipticity and unit conversions.
"""
import numpy as np
from tqdm import tqdm
import erfa
from astropy.time import Time
from astropy.table import Table


def ICRS_to_observed(
    ra, dec,
    mjd_tai, lon, lat, height
):
    obstime = Time(mjd_tai, format='mjd', scale='tai')
    utc1, utc2 = erfa.dtf2d("UTC", *obstime.utc.ymdhms)
    aob, zob, hob, dob, rob, eo = erfa.atco13(
        ra, dec,  # ICRF radec
        0.0, 0.0,  # proper motion
        0.0, 0.0,  # parallax, radial velocity
        utc1, utc2,  # [julian years]
        obstime.delta_ut1_utc,  # [sec]
        lon, lat, height, # [observatory location]
        0.0, 0.0,  # polar motion [rad]
        # Set pressure to 0 to omit refraction
        0.0,  # pressure [hPa = mbar]
        # temp, humid, wave don't matter for zero refraction
        0.0,  # temperature [C]
        0.0,  # relative humidity [0-1]
        1.0  # wavelength [micron]
    )
    return rob - eo, dob


def sphereToCart(ra, dec):
    return np.array([
        np.cos(ra)*np.cos(dec),
        np.sin(ra)*np.cos(dec),
        np.sin(dec)
    ])


def cartToSphere(x, y, z):
    return np.array([
        np.arctan2(y, x),
        np.arcsin(z)
    ])

# Sims has a non-robust version of this, I'm coding up a robust version
# for now, even though it might be slightly non-reproducing.
def rotatePointToPoint(v1, v2):
    axis = np.cross(v1, v2)
    axis /= np.sqrt(np.sum(np.square(axis)))
    #angle = np.arccos(np.dot(v1, v2))
    dsq = np.sum((v1-v2)**2)
    angle = 2. * np.arcsin(0.5 * np.sqrt(dsq))
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    R = np.eye(3) + np.sin(angle)*K + (1-np.cos(angle))*(K@K)
    return R


# Second part of transformation is the field rotator nonsense...
class Deprecessor:
    def __init__(self, orig_ra, orig_dec, mjd_tai, lon, lat, height):
        precessed_ra, precessed_dec = ICRS_to_observed(
            orig_ra, orig_dec,
            mjd_tai, lon, lat, height
        )
        self.precessed_ra = precessed_ra
        self.precessed_dec = precessed_dec

        # find the rotation that carries the original field center
        # to the new field center
        xyz0 = sphereToCart(orig_ra, orig_dec)
        xyz1 = sphereToCart(precessed_ra, precessed_dec)
        first_rotation = rotatePointToPoint(xyz0, xyz1)

        # create a basis set in which the unit vector
        # defining the new field center is the x axis
        xx = np.dot(first_rotation, xyz0)
        rng = np.random.default_rng(99)
        mag = np.NaN
        while np.abs(mag)<1.0e-10 or np.isnan(mag):
            random_vec = rng.normal(size=3)
            comp = np.dot(random_vec, xx)
            yy = random_vec - comp*xx
            mag = np.sqrt(np.sum(np.square(yy)))
            yy /= mag

        zz = np.cross(xx, yy)

        to_self_bases = np.array([xx,
                                  yy,
                                  zz])

        out_of_self_bases = to_self_bases.transpose()

        # Take a point due north of the original field
        # center.  Apply first_rotation to carry it to
        # the new field.  Transform it to the [xx, yy, zz]
        # bases and find the rotation about xx that will
        # make it due north of the new field center.
        # Finally, transform back to the original bases.
        d_dec = np.deg2rad(0.1)
        north = sphereToCart(orig_ra, orig_dec+d_dec)

        north = np.dot(first_rotation, north)

        #print(np.degrees(sphericalFromCartesian(north)))

        north_true = sphereToCart(precessed_ra, precessed_dec+d_dec)

        north = np.dot(to_self_bases, north)
        north_true = np.dot(to_self_bases, north_true)
        north = np.array([north[1], north[2]])
        north /= np.sqrt((north**2).sum())
        north_true = np.array([north_true[1], north_true[2]])
        north_true /= np.sqrt((north_true**2).sum())

        c = north_true[0]*north[0]+north_true[1]*north[1]
        s = north[0]*north_true[1]-north[1]*north_true[0]
        norm = np.sqrt(c*c+s*s)
        c = c/norm
        s = s/norm

        nprime = np.array([c*north[0]-s*north[1],
                           s*north[0]+c*north[1]])

        yz_rotation = np.array([[1.0, 0.0, 0.0],
                                [0.0, c, -s],
                                [0.0, s, c]])

        second_rotation = np.dot(out_of_self_bases,
                                 np.dot(yz_rotation,
                                        to_self_bases))

        self._R = np.dot(second_rotation, first_rotation)

    def __call__(self, ra, dec):
        return cartToSphere(*(self._R.T @ sphereToCart(ra, dec)))


table = Table.read(
    "/global/cfs/cdirs/lsst/groups/fake-source-injection/ssi_vs_resim/"
    "tract_4030.fits"
)

lat = -np.deg2rad(30.2444)
lon = -np.deg2rad(70.7494)
height = 2650.0
pt_ra = np.deg2rad(63.69371156414971580)
pt_dec = np.deg2rad(-34.56012665293710029)
mjd = 59634.08601661110878922  # might be off by 15 sec
obstime = Time(mjd, format='mjd', scale='tai')

raObs, decObs = ICRS_to_observed(
    table['ra'], table['dec'],
    mjd, lon, lat, height
)

deprecessor = Deprecessor(pt_ra, pt_dec, mjd, lon, lat, height)
table['ra_phosim'], table['dec_phosim'] = deprecessor(raObs, decObs)

with open(
     "/global/cfs/cdirs/lsst/groups/fake-source-injection/ssi_vs_resim/"
     "tract_4030_inject.txt",
    'w'
) as fout:
    for irow, row in enumerate(tqdm(table)):
        if row['sourceType'] == 'galaxy':
            bdrat = row['r_bulge_disk_flux_ratio']
            if bdrat == 0.:
                bmag = 50.0
                dmag = row['r_mag']
            else:
                bmag = row['r_mag'] - 2.5*np.log10(bdrat/(1+bdrat))
                dmag = row['r_mag'] - 2.5*np.log10(1/(1+bdrat))

            outstr = "object"
            outstr += f" inject_{irow:012d}_b"
            outstr += f" {np.rad2deg(row['ra_phosim']):.18f}"
            outstr += f" {np.rad2deg(row['dec_phosim']):.18f}"
            outstr += f" {bmag:.10f}"
            outstr += f" flatSED/sed_flat.txt.gz"
            outstr += " 0"  # redshift
            outstr += " 0 0 0"  # gamma1, gamma2, kappa
            outstr += " 0 0"  # dra, ddec
            outstr += " sersic2d"
            outstr += f" {row['bulge_semimajor']:.10f}"
            outstr += f" {row['bulge_semimajor']*row['bulge_axis_ratio']:.10f}"
            outstr += f" {180-row['bulge_pa']:.10f}"
            outstr += f" {row['bulge_n']:.10f}"
            outstr += f" CCM 0.0 3.1 CCM 0.0 3.1\n"
            fout.write(outstr)

            outstr = "object"
            outstr += f" inject_{irow:012d}_d"
            outstr += f" {np.rad2deg(row['ra_phosim']):.18f}"
            outstr += f" {np.rad2deg(row['dec_phosim']):.18f}"
            outstr += f" {dmag:.10f}"
            outstr += f" flatSED/sed_flat.txt.gz"
            outstr += " 0"  # redshift
            outstr += " 0 0 0"  # gamma1, gamma2, kappa
            outstr += " 0 0"  # dra, ddec
            outstr += " sersic2d"
            outstr += f" {row['disk_semimajor']:.10f}"
            outstr += f" {row['disk_semimajor']*row['disk_axis_ratio']:.10f}"
            outstr += f" {180-row['disk_pa']:.10f}"
            outstr += f" {row['disk_n']:.10f}"
            outstr += f" CCM 0.0 3.1 CCM 0.0 3.1\n"
            fout.write(outstr)
