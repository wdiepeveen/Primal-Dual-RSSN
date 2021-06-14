using Manopt

abstract type tProblem end

function getcost(p::P,x) where {P <: tProblem}
    return p.costFunction(x)
end
