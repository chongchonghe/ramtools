# ramtools

Clone this repository to the same folder as your scripts. It will look like this:

```
ramtools
├── README.md
├── plotutils.py
├── ramses.py
├── test
│   └── test.py
├── units.py
├── utilities.py
└── yt_field_descrs.py
your_script.py
```

In your script, import these tools like a Python module. For instance, `from ramtools import ramses` or `from ramtools.plotutils import den_setup`. See examples in test/test.py

Run the test by modifying the path and output id in test.py and do

```
python ramtools/test/test.py
```
