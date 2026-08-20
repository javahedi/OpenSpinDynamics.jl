# OpenSpinDynamics.jl

**OpenSpinDynamics.jl** is a Julia package for simulating quantum spin dynamics in closed and open systems.

It provides:

- sparse quantum spin-model construction
- Krylov and Arnoldi time evolution
- Lindblad master-equation dynamics
- stochastic quantum trajectories
- nearest-neighbor, long-range, and disordered couplings
- common product-state initial conditions

## Quick example

```@example index
using OpenSpinDynamics

N = 4

coupling = NearestNeighborCoupling(0.0, N)

model = SpinModel(
    N;
    Jxy=1.0,
    Jz=1.0,
    coupling=coupling,
)

ψ0 = neel_state(N)

ops = spin_operators(N)

times = collect(
    range(0.0, 2.0; length=51)
)

result = evolve(
    model,
    ψ0,
    times,
    [ops.z[1]];
    method=:krylov,
)

size(result)
```

This evolves a four-spin Néel state under a nearest-neighbor XXZ Hamiltonian and measures the local magnetization
``\langle \sigma_1^z(t) \rangle``.

## Open-system dynamics

The same `SpinModel` interface can be used for dissipative dynamics with Lindblad master equations or stochastic quantum trajectories.

For example, local spontaneous emission can be introduced through collapse operators

```math
L_i = \sqrt{\gamma}\,\sigma_i^-.
```

See the worked examples for complete simulations and plots.

## Documentation

Start with the [Getting Started](getting_started.md) guide for an introduction to the package.

### Examples

- [Closed XXZ dynamics](examples/closed_xxz.md)
- [Amplitude damping](examples/amplitude_damping.md)
- [Dissipative XXZ dynamics](examples/dissipative_xxz.md)
- [Quantum trajectories vs Lindblad dynamics](examples/trajectories_vs_lindblad.md)

See the [API Reference](api.md) for the exported types and functions.
