import argparse, numpy as np, sys
from common import *

parser = argparse.ArgumentParser()
parser.add_argument("n", help="Num of pts", type=int)
parser.add_argument("o", help="Outfile")
parser.add_argument("-f", help="Facility cost", type=float, default=1)
parser.add_argument("-s", help="Random seed", type=int)
parser.add_argument("-d", help="Distribution", choices=["uar","gaussian","mixture"], default="uar")
parser.add_argument("-mixtures", help="Number of mixtures", type=int, default=10)
parser.add_argument("-fmt", help="File format", choices=["simple","pts"], default="pts")
parser.add_argument("--with-extra-seed", help="Generate also 1+seed instance", action='store_true')
parser.add_argument("--show-plot", help="Show plotted", action='store_true')
parser.add_argument("--save-plot", help="Save plotted", action='store_true')
args = parser.parse_args()

def gen(args, params=None):
    np.random.seed(args.s)
    if args.d == "uar":
        pts = np.random.rand(args.n, 2)
    elif args.d == "gaussian":
        pts = np.random.normal((0.5,0.5),(0.25,0.25),(args.n, 2))
    elif args.d == "mixture":
        min_std = 0.1 / args.mixtures**0.5
        max_std = 0.1 / args.mixtures**0.5
        if params is None:
            centers = np.random.rand(args.mixtures, 2)
            stds = [(np.random.uniform(min_std,max_std),np.random.uniform(min_std,max_std)) for _ in centers]
        else:
            centers, stds = params
        params = centers, stds
        pts = np.concatenate([
            np.random.normal(center,
                            std,
                            (args.n // args.mixtures, 2))
            for center,std in zip(centers,stds)])

    plot_2d(pts, args)
    if args.fmt == "pts":
        fl = FLFPts.from_pts(pts, pts, MetricNames.LP_2.value, args.f)
    elif args.fmt == "simple":
        fl = FLFUncap.from_pts(pts, pts, Metric.lp(2), args.f)
    fl.dump(args.o, '\t'.join(sys.argv))
    
    return params

params = gen(args)

if args.with_extra_seed:
    new_seed = 100+args.s
    args.o = args.o.replace("-s%s" % args.s, "-s%s" % new_seed)
    args.s = new_seed
    
    gen(args, params)
