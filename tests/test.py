import sys, os
from time import time
import numpy as np
sys.path.insert(0, "..")
import ramtools as rt

def basic():

    print("Running some tests. By the end of the test, you should see "
        "'Tests done. Successful!'")
    print("\nTest 1: claiming Ramses class on tests/job_test/")
    r = rt.Ramses("tests/Job1")
    print(r.get_time(1))
    r = rt.Ramses(jobid="1", ram_dir="tests")
    print(r.get_time(1))
    rt.set_RAM_DIR("tests")
    r = rt.Ramses(jobid="1")
    print(r.get_time(1))
    print("\nTest 2: load tests/job_test/output_00001/info_00001.txt")
    ds = r.load_ds(1)
    print(f"ds unit length: {ds.length_unit}")

    print("\nTest 3: test utilities.get_sink_mass_from_movie1_for_all_sinks")
    outs, times, sink = rt.utilities.get_sink_mass_from_movie1_for_all_sinks("tests/Job1/movie1")
    print("outs =", outs)
    print("times =", times)
    print("sinks.shape =", np.array(sink).shape)

def plot():
    print("\nTest 4: plotting")
    rp = rt.RamPlot("tests/Job2")
    t1 = time()
    p = rp.plot_prj(40, )
    p.save()
    print("\nFirst time: ", time() - t1)
    t1 = time()
    p = rp.plot_prj(40, )
    p.save()
    print("\nSecond time: ", time() - t1)
    # print(rp.get_time(1))

basic()
if len(sys.argv) >= 2:
    if sys.argv[1] == "plot":
        plot()
print('Tests done. Successful!')
