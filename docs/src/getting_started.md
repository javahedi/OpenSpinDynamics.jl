# Getting Started

OpenSpinDynamics.jl provides a compact interface for simulating closed and open quantum spin dynamics.

A typical workflow has four steps:

1. define the coupling,
2. construct a `SpinModel`,
3. prepare an initial state and observables,
4. call `evolve`.

This page introduces the core public API with small examples.

## Installation

OpenSpinDynamics.jl can be added from its repository while the package is under active development.

```julia
using Pkg

Pkg.add(
    url="https://github.com/javahedi/OpenSpinDynamics.jl",
)
```

Then load the package with

```@example getting_started
using OpenSpinDynamics
using SparseArrays
```

## Build a spin model

As a first example, consider a four-spin nearest-neighbor XXZ chain.

```@example getting_started
N = 4

coupling = NearestNeighborCoupling(
    0.0,
    N,
)

model = SpinModel(
    N;
    Jxy=1.0,
    Jz=0.5,
    coupling=coupling,
)

size(model.hamiltonian)
```

For `N` spin-1/2 degrees of freedom, the Hilbert-space dimension is ``2^N``.

The main model keywords are:

- `Jxy`: transverse exchange strength,
- `Jz`: longitudinal interaction strength,
- `hx`: local ``x`` fields,
- `hz`: local ``z`` fields,
- `coupling`: a coupling object,
- `Jmn`: an explicit coupling matrix for custom interactions.

For example, a model with local fields can be written as

```@example getting_started
field_model = SpinModel(
    N;
    Jxy=1.0,
    Jz=0.5,
    hx=fill(0.2, N),
    hz=zeros(N),
    coupling=coupling,
)

nothing # hide
```

## Couplings

OpenSpinDynamics.jl provides several coupling constructors.

### Nearest-neighbor coupling

```@example getting_started
nearest = NearestNeighborCoupling(
    0.2,
    N,
)

nearest.matrix
```

The parameter `δ` controls alternating nearest-neighbor bond strengths.

### Clean long-range coupling

```@example getting_started
long_range = LongRangeCouplingClean(
    2.0,
    N,
)

long_range.matrix
```

For this coupling,

```math
J_{ij}
\propto
\frac{1}{|i-j|^\alpha}.
```

### Disordered long-range coupling

```@example getting_started
disordered = LongRangeCouplingDisorder(
    2.0,
    20,
    N;
    reordered=false,
)

size(disordered.matrix)
```

Here `N` spins are sampled from a lattice containing `L=20` available sites.

## Initial states

Two common product states are available directly.

### Néel state

```@example getting_started
ψ_neel = neel_state(
    N;
    direction=:z,
)

size(ψ_neel)
```

This produces the alternating state

```math
|\uparrow\downarrow\uparrow\downarrow\cdots\rangle.
```

### Polarized state

```@example getting_started
ψ_polarized = polarized_state(
    N;
    direction=:z,
)

size(ψ_polarized)
```

This produces

```math
|\uparrow\uparrow\uparrow\uparrow\cdots\rangle.
```

Both state constructors also support `direction=:x`.

## Spin operators

Use `spin_operators` to construct single-site observables on the full Hilbert space.

```@example getting_started
ops = spin_operators(N)

size(ops.z[1])
```

The returned `SpinOperators` object contains:

- `ops.x[i]`,
- `ops.z[i]`,
- `ops.plus[i]`,
- `ops.minus[i]`.

For example, the average longitudinal magnetization is

```math
M_z =
\frac{1}{N}
\sum_{i=1}^{N}\sigma_i^z.
```

```@example getting_started
Mz = sum(
    ops.z[i]
    for i in 1:N
) / N

nothing # hide
```

## Closed-system evolution

For pure-state dynamics, pass a wavefunction to `evolve`.

```@example getting_started
times = collect(
    range(0.0, 2.0; length=51)
)

result = evolve(
    model,
    ψ_neel,
    times,
    [Mz];
    method=:krylov,
)

size(result)
```

The rows correspond to time points and the columns correspond to observables.

The closed-system methods are:

- `method=:krylov`,
- `method=:arnoldi`.

For example,

```@example getting_started
result_arnoldi = evolve(
    model,
    ψ_neel,
    times,
    [Mz];
    method=:arnoldi,
)

maximum(
    abs.(
        result[:, 1] .-
        result_arnoldi[:, 1]
    )
)
```

## Lindblad evolution

For Markovian open-system dynamics, provide a density matrix and collapse operators.

Consider local spontaneous emission,

```math
L_i = \sqrt{\gamma}\,\sigma_i^-.
```

```@example getting_started
γ = 0.2

lindblad_ops = [
    sqrt(γ) * ops.minus[i]
    for i in 1:N
]

ρ0 = sparse(
    ψ_neel * ψ_neel'
)

nothing # hide
```

The density matrix can then be evolved with

```@example getting_started
open_result = evolve(
    model,
    ρ0,
    times,
    [Mz];
    method=:expm,
    lindblad_ops=lindblad_ops,
)

size(open_result)
```

The supported density-matrix methods are:

- `method=:expm`,
- `method=:ode`.

For small Hilbert spaces, `:expm` is often a convenient reference method.

## Quantum trajectories

The same dissipative problem can be treated with stochastic quantum trajectories.

```@example getting_started
trajectory_mean, trajectory_std = evolve(
    model,
    ψ_neel,
    times,
    [Mz];
    method=:trajectories,
    lindblad_ops=lindblad_ops,
    num_samples=20,
)

size(trajectory_mean)
```

The trajectory method returns both the sample mean and sample standard deviation.

For accurate Monte Carlo estimates, increase `num_samples`.

## Updating a model

An existing model can be updated in place.

```@example getting_started
update_model!(
    model;
    hz=fill(0.1, N),
)

model.hz
```

Couplings can also be replaced directly:

```@example getting_started
update_model!(
    model;
    coupling=LongRangeCouplingClean(2.0, N),
)

size(model.hamiltonian)
```

Unspecified parameters retain their current values.

## Where to go next

The worked examples explore the main workflows in more detail:

- [Closed XXZ dynamics](examples/closed_xxz.md)
- [Amplitude damping](examples/amplitude_damping.md)
- [Dissipative XXZ chain](examples/dissipative_xxz.md)
- [Quantum trajectories vs Lindblad dynamics](examples/trajectories_vs_lindblad.md)

For all exported functions and types, see the [API Reference](api.md).
