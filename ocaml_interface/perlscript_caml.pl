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

while (<>) {
    s/Manager.ap_manager_ptr/'a t Apron.Manager.t/g; 
    s/internal_ptr/internal/g;  
    s/external elina_poly_/external /g;  
    s/external manager_alloc_loose : unit -> 'a t Apron.Manager.t/external manager_alloc_loose : unit -> loose t Apron.Manager.t/g;
    s/external manager_alloc_strict : unit -> 'a t Apron.Manager.t/external manager_alloc_strict : unit -> strict t Apron.Manager.t/g;
    s/external manager_alloc_equalities : unit -> 'a t Apron.Manager.t/external manager_alloc_equalities : unit -> equalities t Apron.Manager.t/g;
    print;
}
