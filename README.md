# ramtools

`ramtools` is a convenient and fast tool to post-process RAMSES-RT simulations. 

More instructions are available on the [official ramtools documentation site](https://chongchonghe.github.io/ramtools-pages/). 

Author: Chong-Chong He (che1234@umd.edu)

## Install

Clone this repository and do

```
python -m pip install -e .
```
inside the repository directory.

## Usage

```python
import ramtools as rt
r = rt.ramses.Ramses('path/to/a/ramses/job')
...
```



<!-- to the same folder as your scripts. It will look like this: -->

<!-- ``` -->
<!-- ramtools -->
<!-- ├── README.md -->
<!-- ├── plotutils.py -->
<!-- ├── ramses.py -->
<!-- ├── test -->
<!-- │   └── test.py -->
<!-- ├── units.py -->
<!-- ├── utilities.py -->
<!-- └── yt_field_descrs.py -->
<!-- your_script.py -->
<!-- ``` -->

<!-- In your script, import these tools like a Python module. For instance, `from ramtools import ramses` or `from ramtools.plotutils import den_setup`. See examples in test/test.py -->

<!-- Run the test by modifying the path and output id in test.py and do -->

<!-- ``` -->
<!-- python ramtools/test/test.py -->
<!-- ``` -->
