# Stochastic Lipschitz Dynamic Programming
## and Stochastic augmented Dual Dynamic Programming for Lipschitz problems

This package contains an extension to SDDP.jl-v0 for using
[SLDP methods](https://arxiv.org/abs/1905.02290)
to solve Multistage Stochatic MIPs.
It is largely inspired in [SDDiP](https://github.com/lkapelevich/SDDiP.jl),
and our Lagrangian submodule is a very reduced version of the original one.

It requires the user to provide an upper bound $\rho_n$
for the Lipschitz constant of each node.
This has different meanings according to the cutting method.

More documentation is forthcoming.

## Cuts
At present, we use either

* reverse 1-norm cuts
* augmented Lagrangian dual cuts

### Reverse 1-norm cuts

Here, the user must provide a valid upper bound
for the Lipschitz constant (relative to the 1-norm at the domain).
This will be used to construct the cut

$$ \theta \geq Q_n^k(x^k) - \rho_n \|x - x^k \|. $$

### Augmented Lagrangian dual cuts

Here, the user must provide, besides a Lipschitz bound $\rho_n$,
a policy for increasing $\rho_n$.
This is given by the coefficients $a_n$ and $b_n$ for a linear function,
which is then saturated between $0$ and $\rho_n$.
The effective augmenting parameter for node $n$ at iteration $k$ will then be
$p_n = \text{clip}(a_n k + b_n, 0, \rho_n)$.

Then, we obtain a Benders multiplier $\pi$ by solving the LP relaxation,
and solve the augmented Lagrangian problem corresponding to $(\pi, p_n)$.
The cut is therefore

$$ \theta \geq Q_n^{k, \text{AL}}(x^k) + \pi^\top(x - x^k) - p_n \|x - x^k \|. $$

This cut, which we call **Strenghtened augmented Benders cut**,
is valid by construction, so the user-provided $\rho_n$
is only an upper bound to what the _algorithm_ will try.
If it is too small, it will not close the duality gap.

Ideally, one would solve for the optimal Lagrange multiplier,
given the augmenting term $p_n \|x - x^k\|$.
But another method for producing tight cuts
is finding the optimal (lowest) augmenting parameter $\rho$.
For the moment, a simple bisection algorithm for $\rho \in (0, \rho_n)$ is provided.
Tolerances are still hardcoded.
