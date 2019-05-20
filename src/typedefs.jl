# Copyright Bernardo Freitas Paulo da Costa 2019
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
  xhat:: Vector{Float64}
  v   :: Float64
  l   :: Vector{Float64}
  rho :: Float64
end

# Struct for ALD parameters
# Change (slowly) over iterations
mutable struct ALDparams
  Lip     :: Float64
  rho_line:: Tuple{Float64,Float64}
  useSB   :: Bool
  tents   :: Bool
  maxcuts :: Int64
  dropcut :: Int64
  # savecut :: Union{Void,String}
end

# Struct to track ALD
struct ALDExtension
  # The extended state variable information needed for ALD cuts
  states :: Vector{ALDState}
  params :: ALDparams

  # Dynamic information (reset at each backwards pass)
  # xin information (lb/ub might depend on Markov state, current value always),
  # current \rho
  xin_lb :: Vector{Float64}
  xin_v  :: Vector{Float64}
  xin_ub :: Vector{Float64}
  rho    :: Vector{Float64}
  # (v,l,rho) triple,  push!-ed for each noise
  vstore :: Vector{Float64}
  lstore :: Vector{Vector{Float64}}
  rhostore :: Vector{Float64}

  # "Evolution information": list of ALD cuts added
  cuts   :: Vector{ALDCut}
end

aldstates(m) = m.ext[:ALD].states
aldparams(m) = m.ext[:ALD].params
aldcuts(m)   = m.ext[:ALD].cuts

lagrangian(m) = m.ext[:Lagrangian]


immutable MixedSolvers
    LP::JuMP.MathProgBase.AbstractMathProgSolver
    MIP::JuMP.MathProgBase.AbstractMathProgSolver
end
