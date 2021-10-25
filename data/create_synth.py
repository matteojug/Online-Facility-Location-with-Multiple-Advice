import os
from tqdm import tqdm

# For uniform facility cost, using unit cost for the file to avoid storing multiple copies of the same dataset. Can be changed once loaded in the testing code
f = 1
seeds = range(42,42+50) # [42,43,44,45,46]

cmds = []
def exec(cmd):
    cmds.append(cmd)

#synth
out_dir = "final_synth"
os.system("mkdir %s" % out_dir)
for n in [1000,5000,10000]:
    for d in ["mixture"]:
        for mix_cnt in [5,10,25,50,100]:
            for s in seeds:
                fname = f"{out_dir}/{d}-m{mix_cnt}-{n}-f{f}-s{s}"
                exec(f"python3 gen_synth.py {n} {fname} -d {d} -mixtures {mix_cnt} -f {f} -s {s} --save-plot")

for cmd in tqdm(cmds):
    os.system(cmd)