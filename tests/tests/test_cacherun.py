import sys, os
from time import time
import numpy as np
import yt
import ramtools as rt
from glob import glob

def cacherun(job):
    from ramtools import cacherun
    def get_max_pos(ds, field, tag=""):
        return ds.all_data().argmax(field)
    # snap = cacherun.RamsesSnapshot("tests/Job2", 40)
    # t1 = time()
    # print(cacherun.CacheRun(get_max_pos).Run(snap, 'density'))
    # print("First run, time: ", time() - t1)
    # t1 = time()
    # print(cacherun.CacheRun(get_max_pos).Run(snap, 'density'))
    # print("Second run, time: ", time() - t1)
    r = rt.RamsesBase(job)
    out = r.get_first_output()
    if len(glob(f"{job}/output_{out:05d}/hydro_{out:05d}*")) == 0:
        print(f"Simulation output file not complete in {job}/output_{out:05d}. Quitting")
        return
    ds = r.load_ds(out)
    t1 = time()
    print(cacherun.CacheRun(get_max_pos).Run(ds, 'density'))
    print("First run, time: ", time() - t1)
    t1 = time()
    print(cacherun.CacheRun(get_max_pos).Run(ds, 'density'))
    print("Second run, time: ", time() - t1)

def cacherun_base():
    from ramtools.cacherunbase import CacheRun
    import time
    print('start')
    def func1(x, y=1):
        time.sleep(2)
        return x + y
    print('first call:', func1(10, y=2))
    print('second call:', CacheRun(func1, os.path.abspath(__file__), flag='in1')(10, y=2))

def cacherun_base2():
    from ramtools.cacherunbase import CacheRun
    import time
    print('start')
    def func1(x, y=1):
        time.sleep(2)
        return x * x + y
    print('first call:', func1(10, y=2))
    print('second call:', CacheRun(func1, os.path.abspath(__file__), flag='in2')(10, y=2))


if __name__ == "__main__":
    print("---------------- Running test_cacherun.py ---------------------")
    JOB = sys.argv[1]
    def msg(task):
        print("Test {} done!\n".format(task))
    # fastplot(JOB)               # not working
    # fastplot_phase(JOB)         # not working
    # plot(JOB)                   # not working. Replace it with ytfast
    # msg("plot")
    cacherun(JOB)
    msg("cacherun")
    # cacherun_base()
    # msg("cacherun_base")
    # cacherun_base2()
    # msg("cacherun_base2")
    # print('Tests done. Successful!')
    print('Test completed successful!')
    print("-----------------------------------------------------------------")
    print()
