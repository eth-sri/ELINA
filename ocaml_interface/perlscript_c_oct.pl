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

while (<>) {
    s/extern void camlidl_apron_lincons0_ml2c\(value, ap_lincons0_t \*\)/extern void camlidl_apron_lincons0_ml2c\(value, ap_lincons0_t \*, camlidl_ctx\)/g;
    s/\#define camlidl_ml2c_lincons0_ap_lincons0_t\(v,c,ctx\) camlidl_apron_lincons0_ml2c\(v,c\)/\#define camlidl_ml2c_lincons0_ap_lincons0_t\(v,c,ctx\) camlidl_apron_lincons0_ml2c\(v,c,ctx\)/g;
    s/extern void camlidl_apron_tcons0_ml2c\(value, ap_tcons0_t \*\)/extern void camlidl_apron_tcons0_ml2c\(value, ap_tcons0_t \*, camlidl_ctx\)/g;
    s/\#define camlidl_ml2c_tcons0_ap_tcons0_t\(v,c,ctx\) camlidl_apron_tcons0_ml2c\(v,c\)/\#define camlidl_ml2c_tcons0_ap_tcons0_t\(v,c,ctx\) camlidl_apron_tcons0_ml2c\(v,c,ctx\)/g;
    s/struct ap_texpr_op_t/ap_texpr_op_t/g;
    s/struct ap_texpr_rtype_t/ap_texpr_rtype_t/g;
    s/struct ap_texpr_rdir_t/ap_texpr_rdir_t/g;
    print;
}

