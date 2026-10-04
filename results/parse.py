import glob, re
from collections import defaultdict

data = defaultdict(lambda: defaultdict(list))
key = None

for fname in glob.glob("linreg_*.out"):
    for line in open(fname):
        if line.startswith("###"):
            cc, opt, n, p, solver, _ = line.split()[1:]
            key = (cc, opt, int(n), int(p), solver)
        elif key and ":" in line:
            name, val = line.split(":", 1)
            m = re.search(r"[-\d.]+", val)
            if m and name.strip() in ("XtX time", "Solve time", "Total time", "Max |beta - beta_true|"):
                data[key][name.strip()].append(float(m.group()))

print("compiler,opt,N,p,solver,runs,XtX_s,Solve_s,Total_s,max_err")
for k in sorted(data, key=lambda k: (k[0], ["O0","O2","O3","Ofast"].index(k[1]), k[2], k[4])):
    d = data[k]
    avg = lambda x: sum(d[x]) / len(d[x])
    print(f"{k[0]},{k[1]},{k[2]},{k[3]},{k[4]},{len(d['Total time'])},"
          f"{avg('XtX time'):.6f},{avg('Solve time'):.6f},{avg('Total time'):.6f},{avg('Max |beta - beta_true|'):.6f}")
