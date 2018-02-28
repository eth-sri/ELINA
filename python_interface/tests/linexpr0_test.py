from elina_auxiliary_imports import *
from elina_linexpr0 import *
from tests.test_imports import *


def test_set_linexpr_scalar_int(linexpr):

    size = elina_linexpr0_size(linexpr)
    new_size = size // 2

    if not new_size:
        new_size = 1

    new_size = c_size_t(new_size)
    elina_linexpr0_realloc(linexpr, new_size)

    num = c_int(random.randint(0, 99))
    elina_linexpr0_set_cst_scalar_int(linexpr, num)
    for i in range(new_size.value):
        num = c_int(random.randint(0, 99))
        elina_linexpr0_set_coeff_scalar_int(linexpr, ElinaDim(i), num)

    print_c('set scalar int size : {} linexpr: '.format(new_size))
    elina_linexpr0_print(linexpr, None)
    print_c(' is linear: {} is quasilinear: {}\n'.format(elina_linexpr0_is_linear(linexpr),
                                                         elina_linexpr0_is_quasilinear(linexpr)))
    intdim = c_size_t(random.randint(0, new_size.value-1))
    realdim = c_size_t(random.randint(0, new_size.value-1))
    print_c('dim: {} is integer: {} dim: {} is real: {}\n'.format(intdim, elina_linexpr0_is_integer(linexpr, intdim),
                                                                  realdim, elina_linexpr0_is_real(linexpr, realdim)))


def test_set_linexpr_scalar_frac(linexpr):

    size = elina_linexpr0_size(linexpr)
    new_size = size // 2

    if not new_size:
        new_size = c_size_t(1)

    new_size = c_size_t(new_size)
    elina_linexpr0_realloc(linexpr, new_size)

    p = c_long(random.randint(0, 999))
    q = c_ulong(random.randint(1, 20))
    elina_linexpr0_set_cst_scalar_frac(linexpr, p, q)

    for i in range(new_size.value):
        p = c_long(random.randint(0, 999))
        q = c_ulong(random.randint(1, 20))
        elina_linexpr0_set_coeff_scalar_frac(linexpr, ElinaDim(i), p, q)

    print_c('set scalar frac size: {} linexpr: '.format(new_size))
    elina_linexpr0_print(linexpr, None)
    print_c(' is linear: {} is quasilinear: {}\n'.format(elina_linexpr0_is_linear(linexpr),
                                                         elina_linexpr0_is_quasilinear(linexpr)))
    intdim = c_size_t(random.randint(0, new_size.value-1))
    realdim = c_size_t(random.randint(0, new_size.value-1))
    print_c('dim: {} is integer: {} dim: {} is real: {}\n'.format(intdim, elina_linexpr0_is_integer(linexpr, intdim),
                                                                  realdim, elina_linexpr0_is_real(linexpr, realdim)))


def test_set_linexpr_scalar_double(linexpr):

    size = elina_linexpr0_size(linexpr)
    new_size = size // 2

    if not new_size:
        new_size = 1

    new_size = c_size_t(new_size)
    elina_linexpr0_realloc(linexpr, new_size)

    d = c_double(random.uniform(-1, 1))
    elina_linexpr0_set_cst_scalar_double(linexpr, d)
    for i in range(new_size.value):
        d = c_double(random.uniform(-1, 1))
        elina_linexpr0_set_coeff_scalar_double(linexpr, ElinaDim(i), d)

    print_c('set scalar double size : {} linexpr: '.format(new_size))
    elina_linexpr0_print(linexpr, None)
    print_c(' is linear: {} is quasilinear: {}\n'.format(elina_linexpr0_is_linear(linexpr),
                                                         elina_linexpr0_is_quasilinear(linexpr)))
    intdim = c_size_t(random.randint(0, new_size.value-1))
    realdim = c_size_t(random.randint(0, new_size.value-1))
    print_c('dim: {} is integer: {} dim: {} is real: {}\n'.format(intdim, elina_linexpr0_is_integer(linexpr, intdim),
                                                                  realdim, elina_linexpr0_is_real(linexpr, realdim)))


def test_set_linexpr_interval_int(linexpr):

    size = elina_linexpr0_size(linexpr)
    new_size = size // 2

    if not new_size:
        new_size = 1

    new_size = c_size_t(new_size)
    elina_linexpr0_realloc(linexpr, new_size)

    inf = c_int(random.randint(0, 99))
    sup = c_int(random.randint(inf.value, 99 + inf.value))
    elina_linexpr0_set_cst_interval_int(linexpr, inf, sup)
    for i in range(new_size.value):
        inf = c_int(random.randint(0, 99))
        sup = c_int(random.randint(inf.value, 99 + inf.value))
        elina_linexpr0_set_coeff_interval_int(linexpr, ElinaDim(i), inf, sup)

    print_c('set interval int size : {} linexpr: '.format(new_size))
    elina_linexpr0_print(linexpr, None)
    print_c(' is linear: {} is quasilinear: {}\n'.format(elina_linexpr0_is_linear(linexpr),
                                                         elina_linexpr0_is_quasilinear(linexpr)))
    intdim = c_size_t(random.randint(0, new_size.value-1))
    realdim = c_size_t(random.randint(0, new_size.value-1))
    print_c('dim: {} is integer: {} dim: {} is real: {}\n'.format(intdim, elina_linexpr0_is_integer(linexpr, intdim),
                                                                  realdim, elina_linexpr0_is_real(linexpr, realdim)))


def test_set_linexpr_interval_frac(linexpr):

    size = elina_linexpr0_size(linexpr)
    new_size = size // 2

    if not new_size:
        new_size = 1

    new_size = c_size_t(new_size)
    elina_linexpr0_realloc(linexpr, new_size)

    numinf = c_long(random.randint(0, 99))
    deninf = c_ulong(random.randint(1, 20))
    numsup = c_long(numinf.value + random.randint(0, 99))
    densup = deninf
    elina_linexpr0_set_cst_interval_frac(linexpr, numinf, deninf, numsup, densup)
    for i in range(new_size.value):
        numinf = c_long(random.randint(0, 99))
        deninf = c_ulong(random.randint(1, 20))
        numsup = c_long(numinf.value + random.randint(0, 99))
        densup = deninf
        elina_linexpr0_set_coeff_interval_frac(linexpr, ElinaDim(i), numinf, deninf, numsup, densup)

    print_c('set interval frac size : {} linexpr: '.format(new_size))
    elina_linexpr0_print(linexpr, None)
    print_c(' is linear: {} is quasilinear: {}\n'.format(elina_linexpr0_is_linear(linexpr),
                                                         elina_linexpr0_is_quasilinear(linexpr)))
    intdim = c_size_t(random.randint(0, new_size.value-1))
    realdim = c_size_t(random.randint(0, new_size.value-1))
    print_c('dim: {} is integer: {} dim: {} is real: {}\n'.format(intdim, elina_linexpr0_is_integer(linexpr, intdim),
                                                                  realdim, elina_linexpr0_is_real(linexpr, realdim)))


def test_set_linexpr_interval_double(linexpr):

    size = elina_linexpr0_size(linexpr)
    new_size = size // 2

    if not new_size:
        new_size = 1

    new_size = c_size_t(new_size)
    elina_linexpr0_realloc(linexpr, new_size)

    inf = c_double(random.uniform(-1, 1))
    sup = c_double(inf.value + random.uniform(-1, 1))
    elina_linexpr0_set_cst_interval_double(linexpr, inf, sup)
    for i in range(new_size.value):
        inf = c_double(random.uniform(-1, 1))
        sup = c_double(inf.value + random.uniform(-1, 1))
        elina_linexpr0_set_coeff_interval_double(linexpr, ElinaDim(i), inf, sup)

    print_c('set interval double size : {} linexpr: '.format(new_size))
    elina_linexpr0_print(linexpr, None)
    print_c(' is linear: {} is quasilinear: {}\n'.format(elina_linexpr0_is_linear(linexpr),
                                                         elina_linexpr0_is_quasilinear(linexpr)))
    intdim = c_size_t(random.randint(0, new_size.value-1))
    realdim = c_size_t(random.randint(0, new_size.value-1))
    print_c('dim: {} is integer: {} dim: {} is real: {}\n'.format(intdim, elina_linexpr0_is_integer(linexpr, intdim),
                                                                  realdim, elina_linexpr0_is_real(linexpr, realdim)))


def test_linexpr_compare(linexpr1, linexpr2, linexpr3, linexpr4):

    print_c('linexpr1: ')
    elina_linexpr0_print(linexpr1, None)

    print_c('\nlinexpr2: ')
    elina_linexpr0_print(linexpr2, None)

    print_c('\nlinexpr3: ')
    elina_linexpr0_print(linexpr3, None)

    print_c('\nlinexpr4: ')
    elina_linexpr0_print(linexpr4, None)

    print_c('\nlinexpr1 <= linexpr2: {} linexpr2 <= linexpr1: {}'.format(elina_linexpr0_compare(linexpr1, linexpr2),
                                                                         elina_linexpr0_compare(linexpr2, linexpr1)))

    print_c('\nlinexpr1 <= linexpr3: {} linexpr3 <= linexpr1: {}'.format(elina_linexpr0_compare(linexpr1, linexpr3),
                                                                         elina_linexpr0_compare(linexpr3, linexpr1)))

    print_c('\nlinexpr1 <= linexpr4: {} linexpr4 <= linexpr1: {}'.format(elina_linexpr0_compare(linexpr1, linexpr4),
                                                                         elina_linexpr0_compare(linexpr4, linexpr1)))

    print_c('\nlinexpr2 <= linexpr3: {} linexpr3 <= linexpr2: {}'.format(elina_linexpr0_compare(linexpr2, linexpr3),
                                                                         elina_linexpr0_compare(linexpr3, linexpr2)))

    print_c('\nlinexpr2 <= linexpr4: {} linexpr4 <= linexpr2: {}'.format(elina_linexpr0_compare(linexpr2, linexpr4),
                                                                         elina_linexpr0_compare(linexpr4, linexpr2)))

    print_c('\nlinexpr3 <= linexpr4: {} linexpr4 <= linexpr3: {}\n'.format(elina_linexpr0_compare(linexpr3, linexpr4),
                                                                         elina_linexpr0_compare(linexpr4, linexpr3)))


def test_linexpr_equality(linexpr1, linexpr2, linexpr3, linexpr4):

    print_c('linexpr1: ')
    elina_linexpr0_print(linexpr1, None)

    print_c('\nlinexpr2: ')
    elina_linexpr0_print(linexpr2, None)

    print_c('\nlinexpr3: ')
    elina_linexpr0_print(linexpr3, None)

    print_c('\nlinexpr4: ')
    elina_linexpr0_print(linexpr4, None)

    print_c('\nlinexpr1 == linexpr2: {}'.format(elina_linexpr0_equal(linexpr1, linexpr2)))
    print_c('\nlinexpr1 == linexpr3: {}'.format(elina_linexpr0_equal(linexpr1, linexpr3)))
    print_c('\nlinexpr1 == linexpr4: {}'.format(elina_linexpr0_equal(linexpr1, linexpr4)))
    print_c('\nlinexpr2 == linexpr3: {}'.format(elina_linexpr0_equal(linexpr2, linexpr3)))
    print_c('\nlinexpr2 == linexpr4: {}'.format(elina_linexpr0_equal(linexpr2, linexpr4)))
    print_c('\nlinexpr3 == linexpr4: {}\n'.format(elina_linexpr0_equal(linexpr3, linexpr4)))

size = c_size_t(random.randint(3, 22))
linexpr1 = elina_linexpr0_alloc(ElinaLinexprDiscr.ELINA_LINEXPR_SPARSE, size)
linexpr2 = elina_linexpr0_alloc(ElinaLinexprDiscr.ELINA_LINEXPR_DENSE, size)
linexpr3 = elina_linexpr0_alloc(ElinaLinexprDiscr.ELINA_LINEXPR_SPARSE, size)
linexpr4 = elina_linexpr0_alloc(ElinaLinexprDiscr.ELINA_LINEXPR_SPARSE, size)

test_set_linexpr_scalar_int(linexpr1)
test_set_linexpr_scalar_int(linexpr2)

test_set_linexpr_scalar_frac(linexpr1)
test_set_linexpr_scalar_frac(linexpr2)

test_set_linexpr_scalar_double(linexpr1)
test_set_linexpr_scalar_double(linexpr2)

test_set_linexpr_interval_int(linexpr3)
test_set_linexpr_interval_int(linexpr4)

test_set_linexpr_interval_frac(linexpr3)
test_set_linexpr_interval_frac(linexpr4)

test_set_linexpr_interval_double(linexpr3)
test_set_linexpr_interval_double(linexpr4)

test_linexpr_compare(linexpr1, linexpr2, linexpr3, linexpr4)
test_linexpr_equality(linexpr1, linexpr2, linexpr3, linexpr4)
