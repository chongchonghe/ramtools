import sys, os
from time import time
import numpy as np
import yt
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

def fastplot():
    from ramtools import ytfast
    ytfast.set_data_dir("./tests/h5_data")
    print("\nTest 5: fast plotting")
    rp = rt.Ramses("tests/Job2")
    t1 = time()
    p1 = ytfast.ProjectionPlot(
        ds = rp.load_ds(40),
        axis = 2,
        fields = ('gas', 'density'),
        center = 'c',
        width = 0.8,
        axes_unit = 'AU',
        weight_field = ('gas', 'density'),
    )
    p1.save("t.png")
    print("\nFirst time: ", time() - t1)

    print("\nTesting ytfast.SlicePlot")
    r = rt.Ramses(jobid=jobid, ram_dir="..")
    ds = r.load_ds(1)
    t1 = time()
    p = ytfast.SlicePlot(ds, 'x', ('gas', 'density'), width=0.2)
    p.save('f1.png')
    print("\ntime:", time() - t1)
    t1 = time()
    p = ytfast.SlicePlot(ds, 'x', ('gas', 'density'), width=0.2)
    p.save('f2.png')
    print("\ntime:", time() - t1)

def fastplot_phase():
    from ramtools import ytfast
    ytfast.set_data_dir("./tests/h5_data")
    print("\nTest 5: fast plotting")
    rp = rt.Ramses("tests/Job2")
    unitB = rp.unit_B
    def define_field(ds):
        def _Bmag(field, data):
            B2 = ((data[("ramses", "x-Bfield-left")] + data[("ramses", "x-Bfield-right")])**2 +
                (data[("ramses", "y-Bfield-left")] + data[("ramses", "y-Bfield-right")])**2 +
                (data[("ramses", "z-Bfield-left")] + data[("ramses", "z-Bfield-right")])**2) / 4
            return np.sqrt(B2) * unitB.value * 1e6      # in cgs (uGauss)
        ds.add_field(('gas', 'Bmag'), function=_Bmag, sampling_type='cell')
    f_x = ('gas', 'density')
    f_y = ('gas', 'Bmag')
    f_z = ('gas', 'cell_mass')
    t1 = time()
    print("------------------")
    f, ax, p, cb = ytfast.PhasePlot(
        rp.load_ds(40).all_data(),
        f_x,
        f_y,
        [f_z],
        weight_field=None,
        define_field=define_field,
        zlims=[1e31, 1e36]
    )
    print("\nLog 1, time: ", time() - t1)
    f.savefig("t_phaseplot.png", dpi=300)
    print("\nTime: ", time() - t1)

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

def cacherun():
    from ramtools import cacherun
    def get_max_pos(ds, field):
        return ds.all_data().argmax(field)
    # snap = cacherun.RamsesSnapshot("tests/Job2", 40)
    # t1 = time()
    # print(cacherun.CacheRun(get_max_pos).Run(snap, 'density'))
    # print("First run, time: ", time() - t1)
    # t1 = time()
    # print(cacherun.CacheRun(get_max_pos).Run(snap, 'density'))
    # print("Second run, time: ", time() - t1)
    ds = cacherun.my_yt_load("tests/Job2", 40)
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
    if len(sys.argv) == 1:
        main()
    else:
        arg = sys.argv[1]
        if arg == 'all':
            pass
        eval(arg + "()")
    print('Tests done. Successful!')
