import JuMP
import SDDP

function estimatevf(sp::JuMP.Model, ts)
  vf = []
  for c in aldcuts(sp)
    vf_i = [c.v + dot(c.l, ti - c.xi) - c.rho*norm(ti - c.xi, 1) for ti in ts]
    push!(vf, vf_i)
  end
  maximum(hcat(vf...), 2)
end

function qbar(sp::JuMP.Model)
  ex = SDDP.ext(sp)
  if SDDP.hasnoises(sp)
    vs = []
    for i in 1:length(ex.noiseprobability)
      SDDP.setnoise!(sp, ex.noises[i])
      JuMP.solve(sp)
      push!(vs, JuMP.getobjectivevalue(sp))
    end
    return mean(vs)
  else
    JuMP.solve(sp)
    return JuMP.getobjectivevalue(sp)
  end
end


