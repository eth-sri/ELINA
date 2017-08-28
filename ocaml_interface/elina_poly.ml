(* File generated from elina_poly.idl *)

type internal

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

external manager_alloc_loose : unit -> loose t Apron.Manager.t
	= "camlidl_elina_poly_elina_poly_manager_alloc_loose"

external manager_alloc_strict : unit -> strict t Apron.Manager.t
	= "camlidl_elina_poly_elina_poly_manager_alloc_strict"


let manager_is_elina_poly man =
  let str = Apron.Manager.get_library man in
  let str =
    try String.sub str 0 5
    with Invalid_argument _ -> ""
  in
  (String.compare str "elina_poly")==0
let manager_of_elina_poly (man:'a t Apron.Manager.t) : 'b Apron.Manager.t = Obj.magic man
let manager_to_elina_poly (man:'a Apron.Manager.t) : 'b t Apron.Manager.t =
  if manager_is_elina_poly man then
    Obj.magic man
  else
    failwith "Elina_poly.to_elina_poly: the argument manager is not an Elina_poly manager"

let manager_is_elina_poly_loose man =
  let str = Apron.Manager.get_library man in
  (String.compare str "elina_poly, loose mode")==0
let manager_of_elina_poly_loose (man:loose t Apron.Manager.t) : 'a Apron.Manager.t = Obj.magic man
let manager_to_elina_poly_loose (man:'a Apron.Manager.t) : loose t Apron.Manager.t =
  if manager_is_elina_poly_loose man then
    Obj.magic man
  else
    failwith "Elina_poly.to_elina_poly_loose: the argument manager is not a loose Elina_poly manager"

let manager_is_elina_poly_strict man =
  let str = Apron.Manager.get_library man in
  (String.compare str "elina_poly, strict mode")==0
let manager_of_elina_poly_strict (man:strict t Apron.Manager.t) : 'a Apron.Manager.t = Obj.magic man
let manager_to_elina_poly_strict (man:'a Apron.Manager.t) : strict t Apron.Manager.t =
  if manager_is_elina_poly_strict man then
    Obj.magic man
  else
    failwith "Elina_poly.to_elina_poly_strict: the argument manager is not a strict Elina_poly manager"


module Abstract0 = struct
  let is_elina_poly abs =
    manager_is_elina_poly (Apron.Abstract0.manager abs)
  let is_elina_poly_loose abs =
    manager_is_elina_poly_loose (Apron.Abstract0.manager abs)
  let is_elina_poly_strict abs =
    manager_is_elina_poly (Apron.Abstract0.manager abs)

  let of_elina_poly (abs: 'a t Apron.Abstract0.t) : 'b Apron.Abstract0.t = Obj.magic abs
  let of_elina_poly_loose (abs: loose t Apron.Abstract0.t) : 'a Apron.Abstract0.t = Obj.magic abs
  let of_elina_poly_strict (abs: strict t Apron.Abstract0.t) : 'a Apron.Abstract0.t = Obj.magic abs

  let to_elina_poly (abs:'a Apron.Abstract0.t) : 'b t Apron.Abstract0.t =
    if is_elina_poly abs then
      Obj.magic abs
    else
      failwith "Elina_poly.Abstract0.to_elina_poly: the argument value is not a elina_poly value"
  let to_elina_poly_loose (abs:'a Apron.Abstract0.t) : loose t Apron.Abstract0.t =
    if is_elina_poly_loose abs then
      Obj.magic abs
    else
      failwith "Elina_poly.Abstract0.to_elina_poly_loose: the argument value is not a loose elina_poly value"
  let to_elina_poly_strict (abs:'a Apron.Abstract0.t) : strict t Apron.Abstract0.t =
    if is_elina_poly_strict abs then
      Obj.magic abs
    else
      failwith "Elina_poly.Abstract0.to_elina_poly_strict: the argument value is not a strict elina_poly value"
end

module Abstract1 = struct
  let is_elina_poly abs =
    manager_is_elina_poly (Apron.Abstract1.manager abs)
  let is_elina_poly_loose abs =
    manager_is_elina_poly_loose (Apron.Abstract1.manager abs)
  let is_elina_poly_strict abs =
    manager_is_elina_poly (Apron.Abstract1.manager abs)

  let of_elina_poly (abs: 'a t Apron.Abstract1.t) : 'b Apron.Abstract1.t = Obj.magic abs
  let of_elina_poly_loose (abs: loose t Apron.Abstract1.t) : 'a Apron.Abstract1.t = Obj.magic abs
  let of_elina_poly_strict (abs: strict t Apron.Abstract1.t) : 'a Apron.Abstract1.t = Obj.magic abs

  let to_elina_poly (abs:'a Apron.Abstract1.t) : 'b t Apron.Abstract1.t =
    if is_elina_poly abs then
      Obj.magic abs
    else
      failwith "Elina_poly.Abstract1.to_elina_poly: the argument value is not a elina_poly value"
  let to_elina_poly_loose (abs:'a Apron.Abstract1.t) : loose t Apron.Abstract1.t =
    if is_elina_poly_loose abs then
      Obj.magic abs
    else
      failwith "Elina_poly.Abstract1.to_elina_poly_loose: the argument value is not a loose elina_poly value"
  let to_elina_poly_strict (abs:'a Apron.Abstract1.t) : strict t Apron.Abstract1.t =
    if is_elina_poly_strict abs then
      Obj.magic abs
    else
      failwith "Elina_poly.Abstract1.to_elina_poly_strict: the argument value is not a strict elina_poly value"
end

