from elina_auxiliary_imports import *
from ctypes import util
from elina_scalar import *
import random


libc = CDLL(util.find_library('c'))
libmpfr = CDLL(util.find_library('mpfr'))
libgmp = CDLL(util.find_library('gmp'))

printf = libc.printf

cstdout = c_void_p.in_dll(libc, '__stdoutp')

def to_str(str):
    return bytes(str, 'utf-8')

def print_c(str):
    printf(to_str(str))

def test_set_init():

    num = c_long(random.randint(0, 99))
    elina_scalar_set_int(scalar1, num)
    print_c('set int : {} scalar: '.format(num))
    elina_scalar_fprint(cstdout, scalar1)
    print_c('\n')


def test_cmp_int():

    r = c_int(random.randint(0, 99))
    print_c('cmp int scalar: ')
    elina_scalar_fprint(cstdout, scalar1)
    print_c(' int: {} scalar<=int: {} scalar==int: {}\n'.format(r, elina_scalar_cmp_int(scalar1, r), elina_scalar_equal_int(scalar1, r)))


def test_set_mpq():

    mpq_t = Mpq_t()
    libgmp.__gmpq_init(mpq_t)
    p = c_long(random.randint(0, 999))
    q = c_ulong(random.randint(0, 19))
    libgmp.__gmpq_set_si(mpq_t, p, q)
    elina_scalar_set_mpq(scalar1, mpq_t)
    print_c('set mpq: ')
    libgmp.__gmpq_out_str(cstdout, 10, mpq_t)
    print_c(' scalar: ')
    elina_scalar_fprint(cstdout, scalar1)
    print_c('\n')
    libgmp.__gmpq_clear(mpq_t)


def test_set_frac():

    p = c_long(random.randint(0, 99))
    q = c_ulong(random.randint(0, 19))
    elina_scalar_set_frac(scalar1, p, q)
    print_c('set frac p : {} q: {} scalar: '.format(p, q))
    elina_scalar_fprint(cstdout, scalar1)
    print_c('\n')


def test_set_double():

    d = c_double(random.random())
    elina_scalar_set_double(scalar1, d)
    print_c('set double: {} scalar: '.format(d))
    elina_scalar_fprint(cstdout, scalar1)
    print_c('\n')


def test_set_mpfr():

    mpfr_t = Mpfr_t()
    libmpfr.mpfr_init(mpfr_t)
    d = c_double(random.random())
    libmpfr.mpfr_set_d(mpfr_t, d, MpfrRnd.MPFR_RNDU)
    elina_scalar_set_mpfr(scalar1, mpfr_t)
    print_c('set mpfr: {}'.format(d))
    libmpfr.__gmpfr_out_str(cstdout, 10, elina_scalar_print_prec, mpfr_t, MpfrRnd.MPFR_RNDU)
    print_c(' scalar: ')
    elina_scalar_fprint(cstdout, scalar1)
    print_c('\n')
    libmpfr.mpfr_clear(mpfr_t)


def test_mpq_set():

    mpq_t = Mpq_t()
    libgmp.__gmpq_init(mpq_t)
    elina_mpq_set_scalar(mpq_t, scalar1, MpfrRnd.MPFR_RNDU)
    print_c('mpq_set scalar: ')
    elina_scalar_fprint(cstdout, scalar1)
    print_c(' mpq: ')
    libgmp.__gmpq_out_str(cstdout, 10, mpq_t)
    libgmp.__gmpq_clear(mpq_t)
    print_c('\n')


def test_double_set():

    elina_scalar_inv(scalar2, scalar1)
    d = c_double()
    elina_double_set_scalar(byref(d), scalar2, MpfrRnd.MPFR_RNDU)
    print_c('double set scalar: ')
    elina_scalar_fprint(cstdout, scalar1)
    print_c(' double: {}\n'.format(d))


def test_mpfr_set():

    elina_scalar_inv(scalar2, scalar1)
    mpfr_t = Mpfr_t()
    libmpfr.mpfr_init(mpfr_t)
    elina_mpfr_set_scalar(mpfr_t, scalar1, MpfrRnd.MPFR_RNDU)
    print_c('mpfr set scalar: ')
    elina_scalar_fprint(cstdout, scalar1)
    print_c(' mpfr: ')
    libmpfr.__gmpfr_out_str(cstdout, 10, elina_scalar_print_prec, mpfr_t, MpfrRnd.MPFR_RNDU)
    print_c('\n')
    libmpfr.mpfr_clear(mpfr_t)


def test_set_scalar():

    elina_scalar_set(scalar2, scalar1)
    print_c('set scalar1: ')
    elina_scalar_fprint(cstdout, scalar1)
    print_c(' scalar2: ')
    elina_scalar_fprint(cstdout, scalar1)
    print_c(' scalar1==scalar2: {}\n'.format(elina_scalar_equal(scalar1, scalar2)))


def test_inv():
    scalar2 = elina_scalar_alloc_set(scalar1)
    elina_scalar_inv(scalar2, scalar1)
    print_c('inversion scalar1: ')
    elina_scalar_fprint(cstdout, scalar1)
    print_c(' scalar2: ')
    elina_scalar_fprint(cstdout, scalar2)
    print_c(' scalar1 <= scalar2: {}\n'.format(elina_scalar_cmp(scalar1, scalar2)))


def test_equality():

    elina_scalar_neg(scalar1, scalar1)
    print_c('equality scalar1: ')
    elina_scalar_fprint(cstdout, scalar1)
    print_c(' scalar2: ')
    elina_scalar_fprint(cstdout, scalar2)
    print_c(' scalar1==scalar2: {}\n'.format(elina_scalar_equal(scalar1, scalar2)))


def test_sgn():

    elina_scalar_neg(scalar1, scalar1)
    print_c('scalar1: ')
    elina_scalar_fprint(cstdout, scalar1)
    print_c(' scalar2: ')
    elina_scalar_fprint(cstdout, scalar2)
    print_c(' sgn(scalar1): {} sgn(scalar2): {}\n'.format(elina_scalar_sgn(scalar1), elina_scalar_sgn(scalar2)))


# def test_swap()
#
#     elina_scalar_swap(scalar1, scalar1)


def test_infty():

    elina_scalar_set_infty(scalar1, 1)
    elina_scalar_set_infty(scalar2, -1)
    print_c('infty scalar1: ')
    elina_scalar_fprint(cstdout, scalar1)
    print_c(' scalar2:')
    elina_scalar_fprint(cstdout, scalar2)
    print_c(' isinfty(scalar1): {} isinfty(scalar2): {}\n'.format(elina_scalar_infty(scalar1), elina_scalar_infty(scalar2)))



scalar1 = elina_scalar_alloc()
scalar2 = elina_scalar_alloc()

# scalar is set to int
test_set_init()

# # compare with an int
test_cmp_int()

# # scalar is set to mpq
test_set_mpq()
test_inv()
test_mpq_set()

# # scalar is set to frac
test_set_frac()

# # scalar is set to double
test_set_double()
test_double_set()

# scalar is set to mpfr
test_set_mpfr()
test_mpfr_set()

# # set to scalar
test_set_scalar()

# # test for equality and checking the sign of scalar
test_equality()
test_sgn()

# # test for swapping
# test_swap()

# # tests with infinities
test_infty()

elina_scalar_free(scalar2)
elina_scalar_free(scalar1)