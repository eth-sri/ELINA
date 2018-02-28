from tests.test_imports import *
from elina_lincons0 import *
from elina_scalar import *
from elina_dimension import *


def generate_random_lincons0_array(dim, nbcons):

    lincons0_array = elina_lincons0_array_make(nbcons)
    # elina_lincons0_array_fprint(cstdout, lincons0_array, None)

    for i in range(nbcons.value // 3):
        option = random.randint(0, 1)
        if option:
            lincons0_array.p[i].constyp = c_uint(ElinaConstyp.ELINA_CONS_SUPEQ)
        else:
            lincons0_array.p[i].constyp = c_uint(ElinaConstyp.ELINA_CONS_EQ)

        d = c_double(random.uniform(-1, 1))
        linexpr0 = elina_linexpr0_alloc(ElinaLinexprDiscr.ELINA_LINEXPR_SPARSE, dim)
        cst = pointer(linexpr0.contents.cst)
        elina_scalar_set_double(cst.contents.val.scalar, d)

        k = 0
        for j in range(dim.value // 3 + 1):
            linterm = pointer(linexpr0.contents.p.linterm[k])
            linterm.contents.dim = ElinaDim(j)
            coeff = pointer(linterm.contents.coeff)
            d = c_double(random.uniform(-1, 1))
            elina_scalar_set_double(coeff.contents.val.scalar, d)
            k += 1

        elina_linexpr0_realloc(linexpr0, c_size_t(k))
        lincons0_array.p[i].linexpr0 = linexpr0
        # elina_lincons0_fprint(cstdout, lincons0_array.p[i], None)
        # print_c('\n')

    for i in range(nbcons.value // 3, 2 * nbcons.value // 3):

        lincons0_array.p[i].constyp = ElinaConstyp.ELINA_CONS_SUPEQ
        d = c_double(random.uniform(-1, 1))
        linexpr0 = elina_linexpr0_alloc(ElinaLinexprDiscr.ELINA_LINEXPR_SPARSE, dim)
        cst = pointer(linexpr0.contents.cst)
        elina_scalar_set_double(cst.contents.val.scalar, d)

        k = 0
        for j in range(dim.value // 3 + 1, 2 * dim.value // 3):
            linterm = pointer(linexpr0.contents.p.linterm[k])
            linterm.contents.dim = ElinaDim(j)
            coeff = pointer(linterm.contents.coeff)
            d = c_double(random.uniform(-1, 1))
            elina_scalar_set_double(coeff.contents.val.scalar, d)
            k += 1

        elina_linexpr0_realloc(linexpr0, c_size_t(k))
        lincons0_array.p[i].linexpr0 = linexpr0
        # elina_lincons0_fprint(cstdout, lincons0_array.p[i], None)
        # print_c('\n')

    for i in range(2 * nbcons.value // 3, nbcons.value):

        lincons0_array.p[i].constyp = ElinaConstyp.ELINA_CONS_SUPEQ
        d = c_double(random.uniform(-1, 1))
        linexpr0 = elina_linexpr0_alloc(ElinaLinexprDiscr.ELINA_LINEXPR_SPARSE, dim)
        cst = pointer(linexpr0.contents.cst)
        elina_scalar_set_double(cst.contents.val.scalar, d)

        k = 0
        for j in range(2 * dim.value // 3, dim.value):
            linterm = pointer(linexpr0.contents.p.linterm[k])
            linterm.contents.dim = ElinaDim(j)
            coeff = pointer(linterm.contents.coeff)
            d = c_double(random.uniform(-1, 1))
            elina_scalar_set_double(coeff.contents.val.scalar, d)
            k += 1

        elina_linexpr0_realloc(linexpr0, c_size_t(k))
        lincons0_array.p[i].linexpr0 = linexpr0
        # elina_lincons0_fprint(cstdout, lincons0_array.p[i], None)
        # print_c('\n')

    return lincons0_array


dim = c_size_t(random.randint(3, 52))
nbcons = c_size_t(random.randint(3, 52))
lincons0_array = generate_random_lincons0_array(dim, nbcons)

print_c('Lincons array\n')
elina_lincons0_array_fprint(cstdout, byref(lincons0_array), None)
print_c('is linear: {} is quasilinear: {} type of array: {}\n'.format(
    elina_lincons0_array_is_linear(byref(lincons0_array)),
    elina_lincons0_array_is_quasilinear(byref(lincons0_array)),
    elina_lincons0_array_type(byref(lincons0_array))))

perm = elina_dimperm_alloc(dim)
for i in range(dim.value):
    perm.contents.dim[i] = (i+4) % dim.value

perm_array = elina_lincons0_array_permute_dimensions(byref(lincons0_array), perm)
print_c('permutation\n')
elina_dimperm_fprint(cstdout, perm)
print_c('permuted array\n')
elina_lincons0_array_fprint(cstdout, byref(perm_array), None)

intdim = c_size_t(random.randint(1, 5))
realdim = c_size_t(random.randint(1, 5))
dimchange = elina_dimchange_alloc(intdim, realdim)
for i in range(0, intdim.value + realdim.value):
    dimchange.contents.dim[i] = random.randint(0, dim.value)

add_array = elina_lincons0_array_add_dimensions(byref(lincons0_array), dimchange)
print_c('dimension add array\n')
elina_dimchange_fprint(cstdout, dimchange)
print_c('array after adding dimension\n')
elina_lincons0_array_fprint(cstdout, byref(add_array), None)

elina_lincons0_array_clear(byref(lincons0_array))
elina_lincons0_array_clear(byref(perm_array))
elina_lincons0_array_clear(byref(add_array))
