import pathlib
from setuptools import setup, find_packages

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name='ramtools',
    version='0.3',
    description='ramtools: fast analysis of RAMSES data',
    long_description=README,
    long_description_content_type="text/markdown",
    url='https://github.com/chongchonghe/ramtools',
    author='Chong-Chong He',
    author_email='che1234@umd.edu',
    license="MIT",
    # packages=find_packages(),
    # packages=["calc"],
    # entry_points={'console_scripts': ['calc=calc:main']},
    # entry_points={'console_scripts': ['calc=calc.__init__:main']},
    install_requires=['yt>=3.6', 'astropy>=4.0', 'f90nml==1.1.1', ],
)
