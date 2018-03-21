from elina_auxiliary_imports import *
from ctypes import util
import random

libc = CDLL(util.find_library('c'))
libmpfr = CDLL(util.find_library('mpfr'))
libgmp = CDLL(util.find_library('gmp'))

printf = libc.printf

cstdout = c_void_p.in_dll(libc, 'stdout')

def to_str(str):
    return bytes(str, 'utf-8')


def print_c(str):
    printf(to_str(str))
