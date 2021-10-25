import pathlib, subprocess, itertools, argparse, re
from tqdm import tqdm
import multiprocessing.dummy

parser = argparse.ArgumentParser()
parser.add_argument("-m", help="Multi thread", type=int, default=1)
parser.add_argument("-online", default=False, action="store_true")
parser.add_argument("-offline", default=False, action="store_true")
args = parser.parse_args()

def process_online(args):
    fout, fin, f, seed, sorting = args
    extra = get_extra_solutions(fin.as_posix())
    # cmd = f"./solve_online {fin.as_posix()} {f} {seed} {sorting} \"{extra}\" 2>/dev/null" % ()
    cmd = f"./solve_synth {fin.as_posix()} {f} {seed} {sorting} 2>/dev/null" % ()
    out = subprocess.check_output(cmd, shell=True).decode()
    return args, out

def process_offline(args):
    fout, fin, f, seed, sorting = args
    out = subprocess.check_output(f"./solve_offline {fin.as_posix()} {f} {seed} 2>/dev/null", shell=True).decode()
    return args, out

def skip_fname(fname):
    return fname.suffix not in (".png",".sol")

def sort_exp(exp):
    d = int(re.search("/(?:[^-]+-)+([\d]+)-",exp[1].as_posix()).group(1))
    return d

def get_extra_solutions(fin):
    return ""

# seeds = [99,100,101,102,103]
seeds = range(99,99+10)

exps = []
exps.extend(itertools.product(["logs/synth.log"], filter(skip_fname, pathlib.Path("data/final_synth/").glob("*")), [2**i for i in range(8)], seeds, ["none","rand"]))
exps.sort(key=sort_exp)
print("Exps:", len(exps))

done = set()
for p in set(map(lambda x: x[0], exps)):
    with open(p,"a") as _: pass
    with open(p) as fin:
        for l in fin:
            if l.startswith("#"):
                done.add((p, *l.strip().split("\t")[1:5]))

for i,e in enumerate(exps):
    fout, fin, f, seed, sorting = e
    if (fout,fin.as_posix(),str(f),str(seed),sorting) in done:
        exps[i] = None
exps = [*filter(lambda x: x is not None, exps)]
print("Exps after filtering:", len(exps))

p = multiprocessing.dummy.Pool(args.m)
print("Using %d procs" % args.m)
if args.offline:
    exps_single = [*filter(lambda x: x[3] == seeds[0] and x[4].startswith("none"), exps)]
    print("Singleton:", len(exps_single))
    pbi = tqdm(p.imap_unordered(process_offline, exps_single), total=len(exps_single))
    for args_, out in pbi:
        fnout, fin, f, seed, sorting = args_
        with open(fnout,"a") as fout:
            fout.write("#\t%s\t%s\t%s\t%s\n" % (fin.as_posix(), f, seed, sorting))
            fout.write(out)
if args.online:
    pbi = tqdm(p.imap_unordered(process_online, exps), total=len(exps))
    for args_, out in pbi:
        fnout, fin, f, seed, sorting = args_
        with open(fnout,"a") as fout:
            fout.write("#\t%s\t%s\t%s\t%s\n" % (fin.as_posix(), f, seed, sorting))
            fout.write(out)
