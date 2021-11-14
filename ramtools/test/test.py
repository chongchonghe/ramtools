import yt
from ramses.ramtools import ramses
from ramses.ramtools.plotutils import quick_plot_prj, den_setup

# make slice plots of your own
ramses.RAM_DIR = "../2017-RAMSES"
jobid = "2.01.sp.zf5.3"
out = 40
job = ramses.Ramses(jobid, True)
sinks = job.get_sink_positions(out)
ds = job.load_ds(out)
s = yt.SlicePlot(ds, 'x', 'density')
den_setup(s, zlim=[1e0, 1e5])   # set density limit in cm^-3
s.save("ramses/ramtools/test/plots/test-slc-plot.pdf")
print("ramses/ramtools/test/plots/test-slc-plot.pdf saved")


# make projection plots using quick_plot_prj(). Might be slow...

# jobdir = "../Job2.01.sp.zf5.3"
# out = 40
# job = ramses.Ramses(jobdir)
# sinks = job.get_sink_positions(out)
# ds = job.load_ds(out)
# p = quick_plot_prj(ds, sinks, max_level=10, width=0.05)
# p.save("ramtools/test/plots/test-prj-plot.pdf")
