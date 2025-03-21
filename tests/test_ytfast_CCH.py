import sys, os
from time import time
import numpy as np
import ramtools as rt
from ramtools import ytfast

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
    print("\nTime elapsed:", time() - t1)
    
    # # ds = r.load_ds(out)
    # t1 = time()
    # ds2 = r.load_ds(out)
    # p2 = ytfast.ProjectionPlot(ds2, 'x', ('gas', 'density'), center=[0.5, 0.5, 0.5],
    #                           width=0.8, weight_field=('gas', 'density'))
    # p2.set_cmap(('gas', 'density'), 'viridis')
    # rt.plotutils.annotate_box(p2, (.5, .5, .5), 0.1, 'x')
    # p2.save('figures/f1_redo.png')
    # print("\nSecond time:", time() - t1)

if __name__ == "__main__":
    job_path = "/priv/avatar/cche/data/2017-RAMSES/Job-highB/Job.M-C.3B"
    job_idx = 20
    if len(sys.argv) > 2:
        job_path = sys.argv[1]
        job_idx = int(sys.argv[2])
    # else:
    #     print("Usage: python test_ytfast_CCH.py <job_path> <job_idx>")
    #     sys.exit(1)
    os.makedirs("figures", exist_ok=True)
    print(f"---------------- Running {sys.argv[0]} ---------------------")
    fastplot(job_path, job_idx)
    # fastplot(job_path, job_idx)
    print("-----------------------------------------------------------------")
