import SLDP

include("../examples/Caroe/caroe.jl")

m = caroe_model(num_noise=21)
i = 0
n = 21
ts = linspace(0,5,n)
fcf2 = []
for x in ts
    @time for y in ts
        i += 1
        print(i, ",")
        v = SLDP.Qtilde(m,1,1,x,y)
        push!(fcf2, v)
    end
end

