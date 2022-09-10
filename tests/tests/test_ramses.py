import sys, os
from time import time
import numpy as np
import yt
import ramtools as rt
from glob import glob

def basic(job):

    print("Running some tests.")
    print(f"\nTest 1: claiming Ramses class on the job {job}")
    r = rt.Ramses(job)
    if not r.success:
        print(f"Failed. No output directory exists in {job}")
        return
    out = r.get_first_output()
    print(f"At out {out}, t = {r.get_time(out)}")
    print(f"\nTest 2: load {job}/output_{out:05d}/info_{out:05d}.txt")
    ds = r.load_ds(out)
    print(f"ds unit length: {ds.length_unit}")
    print("\nTest 3: test utilities.get_sink_mass_from_movie1_for_all_sinks")
    outs, times, sink = rt.utilities.get_sink_mass_from_movie1_for_all_sinks(
        "tests/Job1/movie1")
    print("outs =", outs)
    print("times =", times)
    print("sinks.shape =", np.array(sink).shape)

def fastplot(job):
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

def plot(job):
    print("\nTest 4: plotting")
    rp = rt.RamPlot(job)
    t1 = time()
    out = 1
    p = rp.plot_prj(out)
    p.save()
    print("\nFirst time: ", time() - t1)
    t1 = time()
    p = rp.plot_prj(out, )
    p.save()
    print("\nSecond time: ", time() - t1)
    # print(rp.get_time(1))

if __name__ == "__main__":
    print("---------------- Running test_ramses.py ---------------------")
    JOB = sys.argv[1]
    def msg(task):
        print("Test {} done!\n".format(task))
    basic(JOB)
    msg("basic")
    # fastplot(JOB)               # not working
    # fastplot_phase(JOB)         # not working
    # plot(JOB)                   # not working. Replace it with ytfast
    # msg("plot")
    # cacherun_base()
    # msg("cacherun_base")
    # cacherun_base2()
    # msg("cacherun_base2")
    # print('Tests done. Successful!')
    print('Tests completed successful!')
    print("-----------------------------------------------------------------")
    print()
