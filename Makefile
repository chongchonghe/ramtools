UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	YT=python
	AMUSE=python
endif
ifeq ($(UNAME_S),Darwin)
	YT=/Users/chongchonghe/anaconda3/envs/yt-git/bin/python
	AMUSE=/Users/chongchonghe/anaconda3/envs/amuse/bin/python
endif

all:

test:
	for f in dynamics/tests/*.py; do \
		echo; \
		echo "Testing $${f}..."; \
		python $$f; \
	done

zoom:
	${YT} dynamics/zoom/main.py

amuse:
	${AMUSE} -m dynamics.cluster.code.main
