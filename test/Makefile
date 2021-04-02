UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	YT=empty
endif
ifeq ($(UNAME_S),Darwin)
	YT=/Users/chongchonghe/anaconda3/envs/yt-git/bin/python
endif

BASE=/Users/chongchonghe/Academics/Projects/2021-ramses
MODULE=ramses/ramtools

test:
	cd $(BASE); $(YT) $(MODULE)/test/test.py
