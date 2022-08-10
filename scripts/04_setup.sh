weekly_version=w_2022_21
LSST_DISTRIB=/cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/${weekly_version}
source "${LSST_DISTRIB}/loadLSST-ext.bash"
setup lsst_distrib
setup -k -r /global/homes/j/jemeyers/src/gen3_workflow
setup -k -r ~/src/pipe_tasks
export OMP_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
wq_env=/global/cfs/cdirs/desc-co/jchiang8/wq_env
export PYTHONPATH=${wq_env}/lib/python3.8/site-packages:${PYTHONPATH}
export PATH=${wq_env}/bin:${PATH}
