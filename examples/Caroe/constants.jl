module Const
# ================
# simulation config
run_names  = ["SB", "ALD parallel"]
ramp_modes = [:None, :parallel]
iters      = [100, 200]

export run_names, ramp_modes, iters
end
