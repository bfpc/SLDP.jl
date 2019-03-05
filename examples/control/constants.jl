module Const
# =============
# default noise
srand(11111)
A_noise   = 0.4
num_noise = 5
noise = randn(num_noise)
noise = A_noise * [noise; -noise]

# =======
# discount
discount = 0.9

# ================
# simulation config
run_names  = ["SB", "ALD simple", "ALD parallel", "ALD parallel2"]
ramp_modes = [:None, :simple, :parallel, :parallel2]
iters      = [100, 200, 200, 200]

export noise, discount, run_names, ramp_modes, iters
end
