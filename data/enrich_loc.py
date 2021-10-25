from common import *
# from collections import defaultdict
from scipy import spatial
from tqdm import tqdm
import sys

loc = sys.argv[1]
city_file = "raw/cities500.txt"

cities = []
data = []
with open(city_file) as fin:
    for l in fin:
        if not l or l.startswith("#"): continue
        row = City(*l.split("\t"))
        cities.append(row)
        data.append((row.latitude, row.longitude))
tree = spatial.KDTree(data)

def closest_city(lat,lon):
    _,k = tree.query([(lat,lon)], k=25)
    best = (1e31, None)
    for i in k[0]:
        d = Metric.geo((lat,lon),(cities[i].latitude, cities[i].longitude))
        if d < best[0]:
            best = (d, cities[i])
    return best

def keep(row):
    if row.loc_id == "00000000000000000000000000000000": return False
    if not (-90 <= row.latitude <= 90): return False
    if not (-180 <= row.longitude <= 180): return False
    return True

with open(loc) as fin, open(loc+".enriched","w") as fout:
    for l in tqdm(fin):
        if not l or l.startswith("#"): continue
        if l.split("\t")[1] == '': continue
        row = Chekin(*l.split("\t"))
        if not keep(row): continue
        cc = closest_city(row.latitude,row.longitude)
        fout.write("%s\t%s\t%s\n" % (l.strip(), cc[1].geonameid, cc[1].country_code))

