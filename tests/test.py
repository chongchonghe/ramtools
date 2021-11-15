import sys
import numpy as np
sys.path.insert(0, "..")
import ramtools as rt

print("Running some tests. By the end of the test, you should see "
      "'Tests done. Successful!'")
print("Test 1: claiming Ramses class on tests/job_test/")
r = rt.Ramses("tests/Job1")
print(r.get_time(1))
r = rt.Ramses(jobid="1", ram_dir="tests")
print(r.get_time(1))
rt.set_RAM_DIR("tests")
r = rt.Ramses(jobid="1")
print(r.get_time(1))
print("Test 2: load tests/job_test/output_00001/info_00001.txt")
ds = r.load_ds(1)
print(f"ds unit length: {ds.length_unit}")

print("Test 3: test utilities.get_sink_mass_from_movie1_for_all_sinks")
outs, times, sink = rt.utilities.get_sink_mass_from_movie1_for_all_sinks("tests/Job1/movie1")
print("outs =", outs)
print("times =", times)
print("sinks.shape =", np.array(sink).shape)

print('Tests done. Successful!')
