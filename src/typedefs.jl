import JuMP
import SDDP

# Struct to store state information for ALD
# xin, xout    are the incoming/outgoing SDDP.@state variables
# av           is the extra variable for calculating absolute values of deviation
# cc           is the constraint xin_{t+1} = xout_t
# av_c1, av_c2 are the constraints to model the absolute value |z - xin| during ALD, when we # relax cc
# lb, ub       are the lower/upper bounds for xout variable, needed for the next stage
struct ALDState
  xin   :: JuMP.Variable
  xout  :: JuMP.Variable
  av    :: JuMP.Variable
  cc    :: SDDP.LinearConstraint
  av_c1 :: SDDP.LinearConstraint
  av_c2 :: SDDP.LinearConstraint
  lb    :: Float64
  ub    :: Float64
end

# Struct to store cut information for ALD
# xi       basepoint for cut
# v,l,rho  cut parameters
struct ALDCut
  xi  :: Vector{Float64}
  v   :: Float64
  l   :: Vector{Float64}
  rho :: Float64
end

# Struct to track ALD
struct ALDExtension
  # The extended state variable information needed for ALD cuts
  states :: Vector{ALDState}

  # Dynamic information (reset at each backwards pass)
  # xin information (lb/ub might depend on Markov state, current value always),
  xin_lb :: Vector{Float64}
  xin_v  :: Vector{Float64}
  xin_ub :: Vector{Float64}
  # (v,l) pairs,  push!-ed for each noise
  vstore :: Vector{Float64}
  lstore :: Vector{Vector{Float64}}

  # "Evolution information": current \rho, list of ALD cuts added
  rho    :: Vector{Float64}
  cuts   :: Vector{ALDCut}
end

aldstates(m) = m.ext[:ALD].states
aldcuts(m)   = m.ext[:ALD].cuts

lagrangian(m) = m.ext[:Lagrangian]


immutable MixedSolvers
    LP::JuMP.MathProgBase.AbstractMathProgSolver
    MIP::JuMP.MathProgBase.AbstractMathProgSolver
end
