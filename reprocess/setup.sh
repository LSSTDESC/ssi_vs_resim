weekly_version=w_2022_24
LSST_DISTRIB=/cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/${weekly_version}
source "${LSST_DISTRIB}/loadLSST-ext.bash"
setup lsst_distrib
setup -k gen3_workflow
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
wq_env=/global/cfs/cdirs/desc-co/jchiang8/wq_env_py_3.10
export PYTHONPATH=${wq_env}/lib/python3.10/site-packages:${PYTHONPATH}
export PATH=${wq_env}/bin:${PATH}
