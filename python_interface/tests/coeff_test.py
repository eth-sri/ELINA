#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
#  This software is distributed under GNU Lesser General Public License Version 3.0.
#  For more information, see the ELINA project website at:
#  http://elina.ethz.ch
#
#  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
#  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
#  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
#  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
#  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY     
#  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
#  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
#  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
#  CONTRACT, TORT OR OTHERWISE).
#
#


import sys
sys.path.insert(0, '../')

from elina_auxiliary_imports import *
from elina_coeff import *
from test_imports import *


def test_set_scalar_int(coeff):
    num = c_long(random.randint(0, 99))
    elina_coeff_set_scalar_int(coeff, num)
    print_c('set scalar int : {} coeff: '.format(num))
    elina_coeff_fprint(cstdout, coeff)
    print_c(' is zero: {}\n'.format(elina_coeff_zero(coeff)))


def test_set_scalar_mpq(coeff):
    mpq_t = Mpq_t()
    libgmp.__gmpq_init(mpq_t)
    p = c_long(random.randint(0, 999))
    q = c_ulong(random.randint(1, 20))
    libgmp.__gmpq_set_si(mpq_t, p, q)
    elina_coeff_set_scalar_mpq(coeff, mpq_t)
    print_c('set scalar mpq: ')
    libgmp.__gmpq_out_str(cstdout, 10, mpq_t)
    print_c(' coeff: ')
    elina_coeff_fprint(cstdout, coeff)
    print_c(' is zero: {}\n'.format(elina_coeff_zero(coeff)))
    libgmp.__gmpq_clear(mpq_t)


def test_set_scalar_frac(coeff):
    p = c_long(random.randint(0, 99))
    q = c_ulong(random.randint(1,20))
    elina_coeff_set_scalar_frac(coeff, p, q)
    print_c('set scalar frac p: {} q: {} coeff '.format(p, q))
    elina_coeff_fprint(cstdout, coeff)
    print_c('\n')


def test_set_scalar_double(coeff):
    d = c_double(random.uniform(-1, 1))
    elina_coeff_set_scalar_double(coeff, d)
    print_c('set scalar double: {} coeff: '.format(d))
    elina_coeff_fprint(cstdout, coeff)
    print_c(' is zero: {}\n'.format(elina_coeff_zero(coeff)))


def test_set_scalar_mpfr(coeff):
    mpfr_t = Mpfr_t()
    libmpfr.mpfr_init(mpfr_t)
    d = c_double(random.uniform(-1, 1))
    libmpfr.mpfr_set_d(mpfr_t, d, MpfrRnd.MPFR_RNDU)
    elina_coeff_set_scalar_mpfr(coeff, mpfr_t)
    print_c('set scalar mpfr: ')
    libmpfr.__gmpfr_out_str(cstdout, 10, elina_scalar_print_prec, mpfr_t, MpfrRnd.MPFR_RNDU)
    print_c(' coeff: ')
    elina_coeff_fprint(cstdout, coeff)
    print_c(' is zero: {}\n'.format(elina_coeff_zero(coeff)))
    libmpfr.mpfr_clear(mpfr_t)


def test_set_interval_int(coeff):
    inf = c_long(random.randint(0, 99))
    sup = c_long(random.randint(0, 99))
    elina_coeff_set_interval_int(coeff, inf, sup)
    print_c('set interval int inf: {} sup: {} coeff: '.format(inf, sup))
    elina_coeff_fprint(cstdout, coeff)
    print_c(' is zero: {}\n'.format(elina_coeff_zero(coeff)))


def test_set_interval_mpq(coeff):
    inf1 = Mpq_t()
    sup1 = Mpq_t()
    libgmp.__gmpq_init(inf1)
    libgmp.__gmpq_init(sup1)
    p = c_long(random.randint(0, 999))
    q = c_ulong(random.randint(1, 20))
    libgmp.__gmpq_set_si(inf1, p, q)
    p = c_long(p.value + random.randint(0, 999))
    libgmp.__gmpq_set_si(sup1, p, q)
    elina_coeff_set_interval_mpq(coeff, inf1, sup1)
    print_c('set interval mpq inf: ')
    libgmp.__gmpq_out_str(cstdout, 10, inf1)
    print_c(' sup: ')
    libgmp.__gmpq_out_str(cstdout, 10, sup1)
    print_c(' coeff: ')
    elina_coeff_fprint(cstdout, coeff)
    print_c(' is zero: {}\n'.format(elina_coeff_zero(coeff)))
    libgmp.__gmpq_clear(inf1)
    libgmp.__gmpq_clear(sup1)


def test_set_interval_frac(coeff):
    p1 = c_long(random.randint(0, 99))
    q1 = c_ulong(random.randint(1, 20))
    p2 = c_long(random.randint(0, 99))
    q2 = c_ulong(random.randint(1, 20))
    elina_coeff_set_interval_frac(coeff, p1, q1, p2, q2)
    print_c('set interval frac p1: {} q1: {} p2: {} q2: {} coeff: '.format(p1, q1, p2, q2))
    elina_coeff_fprint(cstdout, coeff)
    print_c(' is zero: {}\n'.format(elina_coeff_zero(coeff)))


def test_set_interval_double(coeff):
    inf = c_double(random.uniform(-1, 1))
    sup = c_double(random.uniform(-1, 1))
    elina_coeff_set_interval_double(coeff, inf, sup)
    print_c('set interval double inf: {} sup: {} coeff: '.format(inf, sup))
    elina_coeff_fprint(cstdout, coeff)
    print_c(' is zero: {}\n'.format(elina_coeff_zero(coeff)))


def test_set_interval_mpfr(coeff):
    inf = Mpfr_t()
    sup = Mpfr_t()
    libmpfr.mpfr_init(inf)
    libmpfr.mpfr_init(sup)
    d = c_double(random.uniform(-1, 1))
    libmpfr.mpfr_set_d(inf, d, MpfrRnd.MPFR_RNDU)
    d = c_double(random.uniform(-1, 1))
    libmpfr.mpfr_set_d(sup, d, MpfrRnd.MPFR_RNDU)

    elina_coeff_set_interval_mpfr(coeff, inf, sup)
    print_c('set interval mpfr inf: ')
    libmpfr.__gmpfr_out_str(cstdout, 10, elina_scalar_print_prec, inf, MpfrRnd.MPFR_RNDU)
    print_c(' sup: ')
    libmpfr.__gmpfr_out_str(cstdout, 10, elina_scalar_print_prec, sup, MpfrRnd.MPFR_RNDU)
    print_c(' coeff: ')
    elina_coeff_fprint(cstdout, coeff)
    print_c(' is zero: {}\n'.format(elina_coeff_zero(coeff)))
    libmpfr.mpfr_clear(inf)
    libmpfr.mpfr_clear(sup)


def test_coeff_cmp(coeff1, coeff2):
    print_c('Test coeff compare:\n')
    coeff3 = elina_coeff_alloc(ElinaCoeffDiscr.ELINA_COEFF_SCALAR)
    coeff4 = elina_coeff_alloc(ElinaCoeffDiscr.ELINA_COEFF_INTERVAL)

    c = c_int(random.randint(0, 4))
    if c == 0:
        test_set_scalar_int(coeff3)
    elif c == 1:
        test_set_scalar_mpq(coeff3)
    elif c == 2:
        test_set_scalar_frac(coeff3)
    elif c == 3:
        test_set_scalar_double(coeff3)
    else:
        test_set_scalar_mpfr(coeff3)

    c = c_int(random.randint(0, 4))
    if c == 0:
        test_set_interval_int(coeff4)
    elif c == 1:
        test_set_interval_mpq(coeff4)
    elif c == 2:
        test_set_interval_frac(coeff4)
    elif c == 3:
        test_set_interval_double(coeff4)
    else:
        test_set_interval_mpfr(coeff4)

    print_c('cmp scalar vs interval coeff: ')
    elina_coeff_fprint(cstdout, coeff1)
    print_c(' ')
    elina_coeff_fprint(cstdout, coeff2)
    print_c(' coeff1 <= coeff2: {} coeff2 <= coeff1: {}\n'.format(elina_coeff_cmp(coeff1, coeff2),
                                                                  elina_coeff_cmp(coeff2, coeff1)))

    print_c('cmp scalar coeff: ')
    elina_coeff_fprint(cstdout, coeff1)
    print_c(' ')
    elina_coeff_fprint(cstdout, coeff3)
    print_c(' coeff1 <= coeff3: {} coeff3 <= coeff1: {}\n'.format(elina_coeff_cmp(coeff1, coeff3),
                                                                  elina_coeff_cmp(coeff3, coeff1)))

    print_c('cmp interval coeff: ')
    elina_coeff_fprint(cstdout, coeff2)
    print_c(' ')
    elina_coeff_fprint(cstdout, coeff4)
    print_c(' coeff2 <= coeff4: {} coeff4 <= coeff2: {}\n'.format(elina_coeff_cmp(coeff2, coeff4),
                                                                  elina_coeff_cmp(coeff4, coeff2)))


def test_coeff_equality(coeff1, coeff2):
    print_c('Test coeff equality:\n')

    coeff3 = elina_coeff_alloc(ElinaCoeffDiscr.ELINA_COEFF_SCALAR)
    coeff4 = elina_coeff_alloc(ElinaCoeffDiscr.ELINA_COEFF_INTERVAL)

    c = c_int(random.randint(0, 4))
    if c == 0:
        test_set_scalar_int(coeff3)
    elif c == 1:
        test_set_scalar_mpq(coeff3)
    elif c == 2:
        test_set_scalar_frac(coeff3)
    elif c == 3:
        test_set_scalar_double(coeff3)
    else:
        test_set_scalar_mpfr(coeff3)

    c = c_int(random.randint(0, 4))
    if c == 0:
        test_set_interval_int(coeff4)
    elif c == 1:
        test_set_interval_mpq(coeff4)
    elif c == 2:
        test_set_interval_frac(coeff4)
    elif c == 3:
        test_set_interval_double(coeff4)
    else:
        test_set_interval_mpfr(coeff4)

    print_c('equal scalar vs interval coeff: ')
    elina_coeff_fprint(cstdout, coeff1)
    print_c(' ')
    elina_coeff_fprint(cstdout, coeff2)
    print_c(' coeff1 == coeff2: {}\n'.format(elina_coeff_equal(coeff1, coeff2)))

    print_c('equal scalar coeff: ')
    elina_coeff_fprint(cstdout, coeff1)
    print_c(' ')
    elina_coeff_fprint(cstdout, coeff3)
    print_c(' coeff1 == coeff3: {}\n'.format(elina_coeff_equal(coeff1, coeff3)))

    print_c('equal interval coeff: ')
    elina_coeff_fprint(cstdout, coeff2)
    print_c(' ')
    elina_coeff_fprint(cstdout, coeff4)
    print_c(' coeff2 == coeff4: {}\n'.format(elina_coeff_equal(coeff2, coeff4)))


def test_coeff_neg(coeff1):
    print_c('Test coeff neg:\n')
    coeff3 = elina_coeff_alloc(ElinaCoeffDiscr.ELINA_COEFF_SCALAR)
    elina_coeff_neg(coeff3, coeff1)
    print_c('scalar coeff: ')
    elina_coeff_fprint(cstdout, coeff1)
    print_c(' neg coeff: ')
    elina_coeff_fprint(cstdout, coeff3)
    print_c('\n')


def test_coeff_reduce():
    coeff = elina_coeff_alloc(ElinaCoeffDiscr.ELINA_COEFF_INTERVAL)
    num = c_long(random.randint(0,99))
    elina_coeff_set_interval_int(coeff, num, num)
    print_c('Before reduce: ')
    elina_coeff_fprint(cstdout, coeff)
    print_c('\n')
    elina_coeff_reduce(coeff)
    print_c('After reduce: ')
    elina_coeff_fprint(cstdout, coeff)
    print_c('\n')


coeff1 = elina_coeff_alloc(ElinaCoeffDiscr.ELINA_COEFF_SCALAR)
coeff2 = elina_coeff_alloc(ElinaCoeffDiscr.ELINA_COEFF_INTERVAL)

test_set_scalar_int(coeff1)
test_set_scalar_mpq(coeff1)
test_set_scalar_frac(coeff1)
test_set_scalar_double(coeff1)
test_set_scalar_mpfr(coeff1)

test_set_interval_int(coeff2)
test_set_interval_mpq(coeff2)
test_set_interval_frac(coeff2)
test_set_interval_double(coeff2)
test_set_interval_mpfr(coeff2)

test_coeff_cmp(coeff1, coeff2)
test_coeff_equality(coeff1, coeff2)
test_coeff_neg(coeff1)
test_coeff_reduce()
