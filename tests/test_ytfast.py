import sys, os
from time import time
import numpy as np
import yt
import ramtools as rt
from ramtools import ytfast
from glob import glob

def fastplot(job):
    # ytfast._set_data_dir("./tests/h5_data")
    # print("\nTest: ytfast.ProjectionPlot")
    # rp = rt.RamsesBase(job)
    # t1 = time()
    # p1 = ytfast.ProjectionPlot(
    #     ds = rp.load_ds(out),
    #     axis = 2,
    #     fields = ('gas', 'density'),
    #     center = 'c',
    #     width = 0.8,
    #     axes_unit = 'AU',
    #     weight_field = ('gas', 'density'),
    # )
    # p1.save("figures/f1.png")
    # print("\nFirst time: ", time() - t1)

    print("\nTesting ytfast.SlicePlot")
    r = rt.RamsesBase(job)
    if not r.success:
        print(f"Failed to load an output in {job}")
        return
    ds = r.ds1
    out = r.get_first_output()
    t1 = time()
    p = ytfast.SlicePlot(ds, 'x', ('gas', 'density'), width=0.2)
    p.save('figures/f1.png')
    print("\nFirst time:", time() - t1)
    ds = r.load_ds(out)
    t1 = time()
    p = ytfast.SlicePlot(ds, 'x', ('gas', 'density'), width=0.2)
    p.set_cmap(('gas', 'density'), 'prism')
    p.save('figures/f1_redo.png')
    print("\nSecond time:", time() - t1)

def fastplot_phase():
    from ramtools import ytfast
    ytfast._set_data_dir("./tests/h5_data")
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
    rp = rt.RamPlot("Job1")
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
    if not os.path.exists("figures"):
        os.makedirs("figures")
    print(f"---------------- Running {sys.argv[0]} ---------------------")
    #>>>>> old
    # JOB = sys.argv[1]
    # def msg(task):
    #     print("Test {} done!\n".format(task))
    # fastplot(JOB)
    # msg("fastplot")
    # print('Tests completed successful!')
    #<<<<< old
    #>>>>> new
    plot()
    #<<<<< new
    print("-----------------------------------------------------------------")
    print()
