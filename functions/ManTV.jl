include(joinpath("..","functions","ManTVdualbasepoint.jl"))
using Manopt

function ManTV(M::Power,f::PowPoint,α::Real,p::PowPoint,anisotropic=true)
    d = M.powerSize
    if anisotropic || length(d)==1
        tv = costL2TV(M,f,α,p)/α
    else
        datafid = 1/(2*α)*distance(M,f,p)^2

        N = Power(M.manifold,(length(d),))
        R = CartesianIndices(d)
        Λ = forwardLogs(M,p)
        P = ManTVdualbasepoint(p)
        dataTV = 0
        for i in R
            dataTV += norm(N,p[i,:],Λ[i,:])
        end
        tv = datafid + dataTV
    end

    return tv
end
