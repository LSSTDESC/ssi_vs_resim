import subprocess
from astropy.table import Table
from tqdm import tqdm
from collections import defaultdict

subdirs = []

cmd = "hsi ls -l /home/projects/desc/DC2_ImSim/Run2.2i/instCat/"
out = subprocess.run(cmd, shell=True, capture_output=True)
result = out.stderr.decode("utf-8")
result = result.split('\n')
result = result[31:]
for line in result[:-1]:
    subdirs.append(line[60:])

index = []

for subdir in tqdm(subdirs):
    fn_root = f"/home/projects/desc/DC2_ImSim/Run2.2i/instCat/{subdir}/"
    cmd = f"hsi ls -l {fn_root}/*.idx"
    out = subprocess.run(cmd, shell=True, capture_output=True)
    result = out.stderr.decode("utf-8")
    result = result.split('\n')
    result = result[30:-1]
    for line in result:
        begin = int(line[60:68])
        end = int(line[70:78])
        fn = fn_root+line[60:]
        index.append((begin, end, fn))

table = Table.read("data/tract4030ccds.fits")

# Make a big table, 1 row per ccd, with original coaddInputs columns
# and a new column for NERSC HPSS instance catalog
instCats = []
for row in tqdm(table):
    visit = row['visit']
    for min, max, name in index:
        if min <= visit and visit <= max:
            instCats.append(name)
            break
    else:
        raise ValueError(f"Could not find visit {visit}")

table['instCat'] = instCats

table.write("data/tract4030_instCatMap.fits")
