using Aqua
using JET

Aqua.test_all(Chaoslab)
# Report only problems raised inside Chaoslab's own modules (Base and KernelAbstractions internals are not ours).
const OWN_MODULES = (Chaoslab, (m for m in (getfield(Chaoslab, n) for n in names(Chaoslab; all = true))
                                if m isa Module && parentmodule(m) === Chaoslab && m !== Chaoslab)...)
JET.test_package(Chaoslab; target_modules = OWN_MODULES)
