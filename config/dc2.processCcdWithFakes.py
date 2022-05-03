import lsst.meas.extensions.shapeHSM
import lsst.meas.extensions.photometryKron

config.calibrate.measurement.plugins.names=['base_SdssCentroid', 'base_Blendedness', 'base_LocalWcs', 'base_LocalPhotoCalib', 'base_Variance', 'base_CircularApertureFlux', 'base_SkyCoord', 'base_SdssShape', 'base_PsfFlux', 'base_GaussianFlux', 'base_NaiveCentroid', 'base_LocalBackground', 'base_PixelFlags', 'base_FPPosition', 'base_Jacobian', 'ext_photometryKron_KronFlux', 'ext_shapeHSM_HsmSourceMoments', 'ext_shapeHSM_HsmPsfMoments', 'ext_shapeHSM_HsmShapeRegauss']

config.calibrate.astrometry.referenceSelector.doUnresolved = True
config.calibrate.astrometry.referenceSelector.unresolved.name = 'resolved'
config.calibrate.astrometry.referenceSelector.unresolved.minimum = None
config.calibrate.astrometry.referenceSelector.unresolved.maximum = 0.5

# Make sure galaxies are not used for zero-point calculation.
config.calibrate.photoCal.match.referenceSelection.doUnresolved = True
config.calibrate.photoCal.match.referenceSelection.unresolved.name = 'resolved'
config.calibrate.photoCal.match.referenceSelection.unresolved.minimum = None
config.calibrate.photoCal.match.referenceSelection.unresolved.maximum = 0.5

# S/N cuts for zero-point calculation
config.calibrate.photoCal.match.sourceSelection.doSignalToNoise = True
config.calibrate.photoCal.match.sourceSelection.signalToNoise.minimum = 150
config.calibrate.photoCal.match.sourceSelection.signalToNoise.fluxField = 'base_PsfFlux_instFlux'
config.calibrate.photoCal.match.sourceSelection.signalToNoise.errField = 'base_PsfFlux_instFluxErr'
