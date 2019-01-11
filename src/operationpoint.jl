# (C) 2018 authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

using PowerDynBase
using NLsolve

"""
    function operationpoint(start::AbstractState)
    function operationpoint(grid::GridDynamics, start)

Find the operation point (fixed point of the right-hand-side) of a power grid (either `grid` or inferred)
by using the initial value `start`.
"""
function operationpoint(start::AbstractState{G, V, T}) where {G, V, T}
    if SlackAlgebraic ∉ start |> Nodes .|> parametersof .|> typeof
        throw(OperationPointError("currently not making any checks concerning assumptions of whether its possible to find the fixed point"))
    end
    if SwingEq ∈ start |> Nodes .|> parametersof .|> typeof
        throw(OperationPointError("found SwingEq node but these should be SwingEqLVS (just SwingEq is not yet supported for operation point search)"))
    end
    grid = GridDynamics(start)
    rootfct = RootFunction(grid)
    initp = convert(AbstractVector{V}, start)
    res = nlsolve(rootfct, initp)
    if ~res.f_converged
        throw(OperationPointError("solver did not converge"))
    end
    State(grid, res.zero)
end
function operationpoint(grid::GridDynamics, start::AbstractState)
    @assert grid === GridDynamics(start)
    operationpoint(start)
end
function operationpoint(grid::GridDynamics, start::AbstractVector)
    operationpoint(State(grid, start))
end
