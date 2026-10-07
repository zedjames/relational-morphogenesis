"""Local I/O and fixed constants; no parent, annotation or evaluation imports."""
import csv,gzip,hashlib,json,datetime
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
B=1024;REPS=('E7.5-R1','E7.5-R2','E7.5-R3');PAIRS=('P12','P13','P23');ENDS={'P12':(0,1),'P13':(0,2),'P23':(1,2)}
SCALES={'q10':.1,'q25':.25,'q50':.5};OVERLAPS=(('O1','L12','L13'),('O2','L21','L23'),('O3','L31','L32'))
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for v in iter(lambda:f.read(1048576),b''):h.update(v)
    return h.hexdigest()
def read(p):
    p=Path(p)
    with (gzip.open if p.suffix=='.gz' else open)(p,'rt',newline='') as f:return list(csv.DictReader(f))
def write(rel,rows):
    r=list(rows);p=ROOT/rel;p.parent.mkdir(parents=True,exist_ok=True);fields=list(dict.fromkeys(k for x in r for k in x)) or ['status']
    with (gzip.open if p.suffix=='.gz' else open)(p,'wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(r)
def js(rel,o):
    p=ROOT/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(o,indent=2,sort_keys=True,allow_nan=False)+'\n')
def md(rel,s):
    p=ROOT/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(s.rstrip()+'\n')
