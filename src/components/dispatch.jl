# Type-stable dispatch over the ported model types.
#
# An assembled system holds its components in a Vector{AbstractComponent} (case order
# matters: later components read outputs the earlier ones wrote). Calling the kernels
# through that abstract vector costs a dynamic dispatch and boxed views per call, which
# dominated the residual. These helpers branch on the concrete types first, so each call
# is compiled for its model; any other type falls back to dynamic dispatch. Add a new
# model type to both functions.

@inline function _outputs_any!(y, c::AbstractComponent, x, u, rec::ModeRecorder)
    c isa GENROU_PHTRUE && return _outputs!(y, c, x, u, c.p, rec)
    c isa GENSAL_PHTRUE && return _outputs!(y, c, x, u, c.p, rec)
    c isa IEEET1_PHTRUE && return _outputs!(y, c, x, u, c.p, rec)
    c isa IEEEG1_PHTRUE && return _outputs!(y, c, x, u, c.p, rec)
    c isa IEEEG3_PHTRUE && return _outputs!(y, c, x, u, c.p, rec)
    c isa COMPLEXLOAD && return _outputs!(y, c, x, u, c.p, rec)
    return _outputs!(y, c, x, u, params(c), rec)
end

@inline function _step_any!(dx, y, c::AbstractComponent, x, u, rec::ModeRecorder)
    c isa GENROU_PHTRUE && return _step!(dx, y, c, x, u, c.p, rec)
    c isa GENSAL_PHTRUE && return _step!(dx, y, c, x, u, c.p, rec)
    c isa IEEET1_PHTRUE && return _step!(dx, y, c, x, u, c.p, rec)
    c isa IEEEG1_PHTRUE && return _step!(dx, y, c, x, u, c.p, rec)
    c isa IEEEG3_PHTRUE && return _step!(dx, y, c, x, u, c.p, rec)
    c isa COMPLEXLOAD && return _step!(dx, y, c, x, u, c.p, rec)
    return _step!(dx, y, c, x, u, params(c), rec)
end
