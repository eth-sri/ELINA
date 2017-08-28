(* File generated from elina_poly.idl *)

type internal_ptr

(** Convex Polyhedra and Linear Equalities abstract domains *)

type loose
type strict
  (** Two flavors for convex polyhedra: loose or strict.

      Loose polyhedra cannot have strict inequality constraints like [x>0].
      They are algorithmically more efficient
      (less generators, simpler normalization).

      Convex polyhedra are defined by the conjunction of a set of linear
      constraints of the form [a_0*x_0 + ... + a_n*x_n + b >= 0] or
      [a_0*x_0 + ... + a_n*x_n + b > 0]
      where [a_0, ..., a_n, b, c] are constants and [x_0, ..., x_n] variables.
  *)
type 'a t
(** Type of convex polyhedra/linear equalities, where ['a] is [loose], [strict] or [equalities].

    Abstract values which are convex polyhedra have the type
    [(loose t) Apron.Abstract0.t] or [(loose t) Apron.Abstract1.t] or
    [(strict t) Apron.Abstract0.t] or [(strict t) Apron.Abstract1.t].

    Abstract values which are conjunction of linear equalities have the type
    [(equalities t) Apron.Abstract0.t] or [(equalities t) Apron.Abstract1.t].

    Managers allocated by ELINA have the type ['a t Apron.Manager.t].
*)

(** Create an ELINA manager for loose convex polyhedra. *)
external elina_poly_manager_alloc_loose : unit -> Manager.ap_manager_ptr
	= "camlidl_elina_poly_elina_poly_manager_alloc_loose"

(** Create aa ELINA manager for strict convex polyhedra. *)
external elina_poly_manager_alloc_strict : unit -> Manager.ap_manager_ptr
	= "camlidl_elina_poly_elina_poly_manager_alloc_strict"

(** {2 Type conversions} *)

val manager_is_elina_poly : 'a Apron.Manager.t -> bool
val manager_is_elina_poly_loose : 'a Apron.Manager.t -> bool
val manager_is_elina_poly_strict : 'a Apron.Manager.t -> bool


val manager_of_elina_poly : 'a t Apron.Manager.t -> 'b Apron.Manager.t
val manager_of_elina_poly_loose : loose t Apron.Manager.t -> 'a Apron.Manager.t
val manager_of_elina_poly_strict : strict t Apron.Manager.t -> 'a Apron.Manager.t


val manager_to_elina_poly : 'a Apron.Manager.t -> 'b t Apron.Manager.t
val manager_to_elina_poly_loose : 'a Apron.Manager.t -> loose t Apron.Manager.t
val manager_to_elina_poly_strict : 'a Apron.Manager.t -> strict t Apron.Manager.t


module Abstract0 : sig
  val is_elina_poly : 'a Apron.Abstract0.t -> bool
  val is_elina_poly_loose : 'a Apron.Abstract0.t -> bool
  val is_elina_poly_strict : 'a Apron.Abstract0.t -> bool
  
  val of_elina_poly : 'a t Apron.Abstract0.t -> 'b Apron.Abstract0.t
  val of_elina_poly_loose : loose t Apron.Abstract0.t -> 'a Apron.Abstract0.t
  val of_elina_poly_strict : strict t Apron.Abstract0.t -> 'a Apron.Abstract0.t
  
  val to_elina_poly : 'a Apron.Abstract0.t -> 'b t Apron.Abstract0.t
  val to_elina_poly_loose : 'a Apron.Abstract0.t -> loose t Apron.Abstract0.t
  val to_elina_poly_strict : 'a Apron.Abstract0.t -> strict t Apron.Abstract0.t
end

module Abstract1 : sig
  val is_elina_poly : 'a Apron.Abstract1.t -> bool
  val is_elina_poly_loose : 'a Apron.Abstract1.t -> bool
  val is_elina_poly_strict : 'a Apron.Abstract1.t -> bool

  val of_elina_poly : 'a t Apron.Abstract1.t -> 'b Apron.Abstract1.t
  val of_elina_poly_loose : loose t Apron.Abstract1.t -> 'a Apron.Abstract1.t
  val of_elina_poly_strict : strict t Apron.Abstract1.t -> 'a Apron.Abstract1.t

  val to_elina_poly : 'a Apron.Abstract1.t -> 'b t Apron.Abstract1.t
  val to_elina_poly_loose : 'a Apron.Abstract1.t -> loose t Apron.Abstract1.t
  val to_elina_poly_strict : 'a Apron.Abstract1.t -> strict t Apron.Abstract1.t
end


(**
{2 Compilation information}

See {!Introduction.compilation} for complete explanations. 
We just show examples with the file [mlexample.ml].

{3 Bytecode compilation}

{[ocamlc -I $MLGMPIDL_PREFIX/lib -I $APRON_PREFIX/lib -o mlexample.byte \
  bigarray.cma gmp.cma apron.cma elina_poly.cma mlexample.ml]}

{[ocamlc -I $MLGMPIDL_PREFIX/lib -I $APRON_PREFIX/lib -make-runtime -o myrun \
  bigarray.cma gmp.cma apron.cma elina_poly.cma

ocamlc -I $MLGMPIDL_PREFIX/lib -I $APRON_PREFIX/lib -use-runtime myrun -o mlexample.byte \
  bigarray.cma gmp.cma apron.cma elina_poly.cma mlexample.ml ]}

{3 Native-code compilation}

{[ocamlopt -I $MLGMPIDL_PREFIX/lib -I $APRON_PREFIX/lib -o mlexample.opt \
  bigarray.cmxa gmp.cmxa apron.cmxa elina_poly.cmxa mlexample.ml ]}

{3 Without auto-linking feature}

{[ocamlopt -I $MLGMPIDL_PREFIX/lib -I $APRON_PREFIX/lib -noautolink -o mlexample.opt \
  bigarray.cmxa gmp.cmxa apron.cmxa elina_poly.cmxa mlexample.ml \
  -cclib "-L$MLGMPIDL_PREFIX/lib -L$APRON_PREFIX/lib \
	  -loptpoly_caml_debug -loptpoly \
	  -lapron_caml_debug -lapron_debug \
	  -lgmp_caml -L$MPFR_PREFIX/lib -lmpfr -L$GMP_PREFIX/lib -lgmp \
	  -L$CAMLIDL_PREFIX/lib/ocaml -lcamlidl \
	  -lbigarray" ]}

*)
