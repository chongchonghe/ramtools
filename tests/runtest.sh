#!/bin/bash

set -e

python test_ytfast_CCH.py /priv/avatar/cche/data/2017-RAMSES/Job-highB/Job.M-C.3B 20

exit 0

export PYTHONPATH=..
JOB="$1"
if [[ -z $JOB ]]; then
    JOB="./Job1"
fi
if [[ "$JOB" == "local" ]]; then
  JOB="/Users/chongchonghe/Ext_no_timemachine/Sam/Job2.0.v2.full"
fi
for f in tests/*.py; do
     python $f $JOB;
done

# "/priv/avatar/cche/data/2017-RAMSES/Job-highB/Job.M-C.3B/
# fastplot("/startrek/chongchong/Projects/2017-RAMSES/Job2.2.2.v2/", 19)