#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
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


from elina_auxiliary_imports import *
from elina_interval import *
from tests.test_imports import *

def test_set_int(inteval1):

    inf = c_long(random.randint(0, 99))
    sup = c_long(random.randint(0, 99))
    elina_interval_set_int(interval1, inf, sup)

    print_c('set int inf: {} sup: {} interval: '.format(inf, sup))
    elina_interval_fprint(cstdout, interval1)
    print_c(' is bottom: {} is top: {}\n'.format(elina_interval_is_bottom(interval1), elina_interval_is_top(interval1)))


def test_set_mpq(interval1):

    inf1 = Mpq_t()
    sup1 = Mpq_t()
    libgmp.__gmpq_init(inf1)
    libgmp.__gmpq_init(sup1)
    p = c_long(random.randint(0, 999))
    q = c_ulong(random.randint(0, 19))
    libgmp.__gmpq_set_si(inf1, p, q)
    p = c_long(p.value + random.randint(0, 999))
    libgmp.__gmpq_set_si(sup1, p, q)
    elina_interval_set_mpq(interval1, inf1, sup1)

    print_c('set mpq inf: ')
    libgmp.__gmpq_out_str(cstdout, 10, inf1)
    print_c(' sup: ')
    libgmp.__gmpq_out_str(cstdout, 10, sup1)
    print_c(' interval: ')
    elina_interval_fprint(cstdout, interval1)
    print_c(' is bottom: {} is top: {}\n'.format(elina_interval_is_bottom(interval1), elina_interval_is_top(interval1)))
    libgmp.__gmpq_clear(inf1)
    libgmp.__gmpq_clear(sup1)


def test_set_frac(interval1):

    p1 = c_long(random.randint(0, 99))
    q1 = c_ulong(random.randint(1, 19))
    p2 = c_long(random.randint(0, 99))
    q2 = c_ulong(random.randint(1, 19))
    elina_interval_set_frac(interval1, p1, q1, p2, q2)

    print_c('set frac p1: {} q1: {} p2: {} q2: {} interval: '.format(p1, q1, p2, q2))
    elina_interval_fprint(cstdout, interval1)
    print_c(' is bottom: {} is top: {}\n'.format(elina_interval_is_bottom(interval1), elina_interval_is_top(interval1)))


def test_set_double(interval1):

    inf = c_double(random.uniform(-1, 1))
    sup = c_double(random.uniform(-1, 1))
    elina_interval_set_double(interval1, inf, sup)

    print_c('set double inf: {} sup: {} interval: '.format(inf, sup))
    elina_interval_fprint(cstdout, interval1)
    print_c(' is bottom: {} is top: {}\n'.format(elina_interval_is_bottom(interval1), elina_interval_is_top(interval1)))


def test_set_mpfr(interval1):

    inf = Mpfr_t()
    sup = Mpfr_t()
    libmpfr.mpfr_init(inf)
    libmpfr.mpfr_init(sup)
    d = c_double(random.uniform(-1, 1))
    libmpfr.mpfr_set_d(inf, d, MpfrRnd.MPFR_RNDU)
    d = c_double(random.uniform(-1, 1))
    libmpfr.mpfr_set_d(sup, d, MpfrRnd.MPFR_RNDU)
    elina_interval_set_mpfr(interval1, inf, sup)

    print_c('set mpfr inf: ')
    libmpfr.__gmpfr_out_str(cstdout, 10, elina_scalar_print_prec, inf, MpfrRnd.MPFR_RNDU)
    print_c(' sup: ')
    libmpfr.__gmpfr_out_str(cstdout, 10, elina_scalar_print_prec, sup, MpfrRnd.MPFR_RNDU)
    print_c(' interval: ')
    elina_interval_fprint(cstdout, interval1)
    print_c(' is bottom: {} is top: {}\n'.format(elina_interval_is_bottom(interval1), elina_interval_is_top(interval1)))


def test_set_interval(interval1, interval2):
    elina_interval_set(interval2, interval1)
    print_c('set interval1: ')
    elina_interval_fprint(cstdout, interval1)
    print_c(' interval2: ')
    elina_interval_fprint(cstdout, interval2)
    print_c(' interval1 == interval2: {}\n'.format(elina_interval_equal(interval1, interval2)))


def test_cmp(interval1, interval2):
    inf2 = c_long(random.randint(0, 99))
    sup2 = c_long(inf2.value + random.randint(0, 99))
    elina_interval_set_int(interval2, inf2, sup2)

    print_c('cmp interval1: ')
    elina_interval_fprint(cstdout, interval1)
    print_c(' interval2: ')
    elina_interval_fprint(cstdout, interval2)
    print_c(' interval1 <= interval2: {} interva1 == interval2: {}\n'.format(elina_interval_cmp(interval1, interval2),
                                                                           elina_interval_equal(interval1, interval2)))

def test_equality(interval1, interval2):
    elina_interval_set_bottom(interval1)
    print_c('equality interval1: ')
    elina_interval_fprint(cstdout, interval1)
    elina_interval_set_top(interval2)
    print_c(' interval2: ')
    elina_interval_fprint(cstdout, interval2)
    print_c(' interval1 == interval2: {}\n'.format(elina_interval_equal(interval1, interval2)))


def test_neg(interval1, interval2):
    elina_interval_neg(interval1, interval2)
    print_c('neg interval1: ')
    elina_interval_fprint(cstdout, interval1)
    print_c(' interval2: ')
    elina_interval_fprint(cstdout, interval2)
    print_c('\ninterval1 is bottom: {} is top: {}'.format(elina_interval_is_bottom(interval1), elina_interval_is_top(interval1)))
    print_c('\ninterval2 is bottom: {} is top: {}\n'.format(elina_interval_is_bottom(interval2), elina_interval_is_top(interval2)))


def test_interval_array():
    print_c('Testing arrays\n')
    size = c_size_t(random.randint(1, 10))
    interval_array = elina_interval_array_alloc(size)
    for i in range(size.value):
        option = random.randint(0, 4)

        if option == 0:
            test_set_int(interval_array[i])
        elif option == 1:
            test_set_mpq(interval_array[i])
        elif option == 2:
            test_set_frac(interval_array[i])
        elif option == 3:
            test_set_double(interval_array[i])
        else:
            test_set_mpfr(interval_array[i])


interval1 = elina_interval_alloc()
interval2 = elina_interval_alloc()

test_set_int(interval1)
test_set_mpq(interval1)
test_set_frac(interval1)
test_set_double(interval1)
test_set_mpfr(interval1)
test_set_interval(interval1, interval2)
test_cmp(interval1, interval2)
test_equality(interval1, interval2)
test_neg(interval1, interval2)
test_interval_array()

elina_interval_free(interval1)
elina_interval_free(interval2)
