from glob import glob
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

def readme() :
  with open('README.rst') as f :
    return f.read()

sources = glob("*.pyx")

setup (
  name="pycontacts",
  version="1.0",
  description = "Contact analysis on state lists",
  long_description=readme(),
  platforms=['linux'],
  author="Samuel Genheden",
  author_email="samuel.genheden@gmail.com",
  url="https://github.com/sgenheden/Scripts/pycontacts",
  download_url="https://github.com/sgenheden/Scripts/pycontacts",
  setup_requires=['setuptools_cython','Cython >= 0.10'],
  cmdclass = {'build_ext': build_ext},
  ext_modules = [
    Extension("pycontacts",sources,include_dirs=[numpy.get_include()])
  ]
)
