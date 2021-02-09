from ramtools import ramses
from ramtools.plotutils import quick_plot_prj

jobdir = "../Job2.01.sp.zf5.3"
out = 40
job = ramses.Ramses(jobdir)
sinks = job.get_sink_positions(out)
ds = job.load_ds(out)
p = quick_plot_prj(ds, sinks, max_level=10)
p.save("ramtools/test/test-quick_plot_prj.pdf")
