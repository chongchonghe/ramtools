import sys, os
from time import time
import numpy as np
import yt
import ramtools as rt
from ramtools import ytfast
from glob import glob

def fastplot(jobpath, out):
    print("\nTesting ytfast.SlicePlot")
    r = rt.RamsesBase(jobpath)
    if not r.success:
        print(f"Failed to load an output in {jobpath}")
        return
    ds = r.load_ds(out)
    t1 = time()
    p = ytfast.ProjectionPlot(ds, 'x', ('gas', 'density'), center=[0.5, 0.5, 0.5],
                              width=0.8, weight_field=('gas', 'density'))
    p.save('figures/f1.png')
    print("\nFirst time:", time() - t1)
    ds = r.load_ds(out)
    t1 = time()
    p = ytfast.ProjectionPlot(ds, 'x', ('gas', 'density'), center=[0.5, 0.5, 0.5],
                              width=0.8, weight_field=('gas', 'density'))
    p.set_cmap(('gas', 'density'), 'prism')
    rt.plotutils.annotate_box(p, (.5, .5, .5), 0.1, 'x')
    p.save('figures/f1_redo.png')
    print("\nSecond time:", time() - t1)

if __name__ == "__main__":
    if not os.path.exists("figures"):
        os.makedirs("figures")
    print(f"---------------- Running {sys.argv[0]} ---------------------")
    # fastplot("/startrek/chongchong/Projects/2017-RAMSES/Job2.01/")
    fastplot("/startrek/chongchong/Projects/2017-RAMSES/Job2.2.2.v2/", 19)
    print("-----------------------------------------------------------------")
    print()
