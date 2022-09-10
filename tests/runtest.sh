#/bin/zsh

set -e

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
