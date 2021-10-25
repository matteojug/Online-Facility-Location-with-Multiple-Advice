import pathlib, subprocess, itertools, argparse, re
from tqdm import tqdm
import multiprocessing.dummy

parser = argparse.ArgumentParser()
parser.add_argument("-m", help="Multi thread", type=int, default=1)
parser.add_argument("-online", default=False, action="store_true")
parser.add_argument("-offline", default=False, action="store_true")
parser.add_argument("-redo-skipped", default=False, action="store_true")
args = parser.parse_args()

assert not (args.offline and args.online), "Run first -offline, then -online"

def process_online(args):
    fout, fin, f, seed, sorting = args
    extra = get_extra_solutions(fin.as_posix())
    if extra is None: return args, None
    cmd = f"./solve_online {fin.as_posix()} {f} {seed} {sorting} \"{extra}\" 2>/dev/null" % ()
    # print(cmd); return args, None
    out = subprocess.check_output(cmd, shell=True).decode()
    return args, out

def process_offline(args):
    fout, fin, f, seed, sorting = args
    out = subprocess.check_output(f"./solve_offline {fin.as_posix()} {f} {seed} 2>/dev/null", shell=True).decode()
    return args, out


def skip_fname_ext(fname):
    return fname.suffix not in (".png",".sol")

def skip_fname(fname):
    g = re.search(r"-D(\d+)(_\d+)?[^/]*$", fname.name)
    # if g.group(2) is None: return False # This
    if g.group(2) is not None: return False # This
    return skip_fname_ext(fname)

def get_extra_solutions(fin):
    # -1d=data/gen_pts_day_all/gw-us-D98;-7d=data/gen_pts_day_all/gw-us-D92;-30d=data/gen_pts_day_all/gw-us-D0;
    g = re.search(r"-D(\d+)(_\d+)?[^/]*$", fin)
    day, timew = g.group(1), g.group(2)
    s = []
    if timew is None:
        timew = 0
        for delta in [-1, -2, -3, -4, -5, -6, -7, -8]:
        # for delta in [-4,-3,-2,-1]:
            new_day = int(day)
            new_time = timew + delta
            while new_time < 0:
                new_day -= 1
                new_time += 4 # TODO: this
            if new_day < 0: continue
            ss = "%sstp=%s" % (delta, fin.replace(f"-D{day}", f"-D{new_day}_{new_time}"))
            s.append(ss)
    else:
        return None

    return ";".join(s)

seeds = range(99,99+10)

exps = []
if args.offline:
    exps.extend(itertools.product(["logs/days.log"], sorted(filter(skip_fname_ext, pathlib.Path("data/final_days/").glob("*"))), [100,500,1000,5000,10000,50000], seeds, ["none"]))
    exps.extend(itertools.product(["logs/days.log"], sorted(filter(skip_fname_ext, pathlib.Path("data/final_uber/").glob("*"))), [100,500,1000,5000,10000,50000], seeds, ["none"]))

if args.online:
    exps.extend(itertools.product(["logs/days.log"], sorted(filter(skip_fname, pathlib.Path("data/final_days/").glob("*"))), [100,500,1000,5000,10000,50000], seeds, ["none"]))
    exps.extend(itertools.product(["logs/days.log"], sorted(filter(skip_fname, pathlib.Path("data/final_uber/").glob("*"))), [100,500,1000,5000,10000,50000], seeds, ["none"])) 

print("Exps:", len(exps))

done = set()
for p in set(map(lambda x: x[0], exps)):
    with open(p,"a") as _: pass
    with open(p) as fin:
        lines = fin.readlines()
        for i,l in enumerate(lines):
            if not l.startswith("#"): continue
            if i+1 < len(lines) and "SKIPPED_TOO_BIG" in lines[i+1] and args.redo_skipped: continue
            done.add((p, *l.strip().split("\t")[1:5]))

for i,e in enumerate(exps):
    fout, fin, f, seed, sorting = e
    if (fout,fin.as_posix(),str(f),str(seed),sorting) in done:
        # pass
        exps[i] = None
exps = [*filter(lambda x: x is not None, exps)]
print("Exps after filtering:", len(exps))

p = multiprocessing.dummy.Pool(args.m)
print("Using %d procs" % args.m)
if args.offline:
    exps_single = [*filter(lambda x: x[3] == seeds[0] and x[4].startswith("none"), exps)]
    print("Singleton:", len(exps_single))
    skipped_big = 0
    pbi = tqdm(p.imap_unordered(process_offline, exps_single), total=len(exps_single))
    for args_, out in pbi:
        fnout, fin, f, seed, sorting = args_
        with open(fnout,"a") as fout:
            fout.write("#\t%s\t%s\t%s\t%s\n" % (fin.as_posix(), f, seed, sorting))
            fout.write(out)
        if "SKIPPED_TOO_BIG" in out: skipped_big += 1
    print("Skipped:",skipped_big)
if args.online:
    skipped_big = 0
    pbi = tqdm(p.imap_unordered(process_online, exps), total=len(exps))
    for args_, out in pbi:
        if out is None:
            print(args_, "got None")
            continue
        
        fnout, fin, f, seed, sorting = args_
        with open(fnout,"a") as fout:
            fout.write("#\t%s\t%s\t%s\t%s\n" % (fin.as_posix(), f, seed, sorting))
            fout.write(out)
        if "SKIPPED_TOO_BIG" in out: skipped_big += 1
    print("Skipped:",skipped_big)
