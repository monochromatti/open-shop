import OpenSHOP
import Lean.Util.CollectAxioms

open Lean Elab Command in
run_cmd do
  let env ← getEnv
  let mut checked := 0
  for (name, info) in env.constants.toList do
    if name.getRoot == `OpenSHOP then
      let axioms ← collectAxioms name
      let unexpected := axioms.filter fun ax =>
        ax != ``propext && ax != ``Classical.choice && ax != ``Quot.sound
      unless unexpected.isEmpty do
        throwError "{name} depends on unapproved axioms: {unexpected}"
      if info matches .thmInfo _ then
        checked := checked + 1
  if checked == 0 then
    throwError "No model theorems were checked"
  logInfo m!"Checked {checked} theorem declarations, including generated auxiliaries; only standard Lean logical axioms are used."
