"""
This imports all of the filters within the FilterClasses file. This is very
important as the filters will not work if this doesn't exist.

This is a dynamic and modular way of doing these imports.

Code is taken from:
https://stackoverflow.com/questions/1057431/how-to-load-all-modules-in-a-folder
"""


from os.path import dirname, basename, isfile
import glob

modules = glob.glob(dirname(__file__) + "/*.py")
__all__ = [basename(f)[:-3] for f in modules if isfile(f) and not f.endswith("__init__.py")]

from . import *
