# Sensitivity Approximation in Julia

The motivation for this reference implementation is an approximation
of the Fisher Information FI matrix for optimization and sampling
methods. For these purposes we require a fast FI that does not need to
be accurate to machine precision.

The collection of julia scripts in this repository is supplemental
material to the publication:
```
Eriksson, O., Kramer, A., Milinanni, F., and Nyquist, P., “Sensitivity Approximation by the Peano-Baker Series”, arXiv e-prints, 2021
```
published as pre-print [arXiv:2109.00067 [math.NA]](https://arxiv.org/abs/2109.00067)

For citation, you can also use the BibTeX entry in [SensApprox.bib](./SensApprox.bib) (or directly from arXiv, linked above).

## Brief Summary

We use an approximative method with an error that scales with `h²`,
where `h` is the maximum step-size of the integration method (steps
are adaptive).

For initial value problems that converge to a stable steady state
(node or focus and not e.g. stable limit cycles), the method switches
from Peano Baker Series approximation to near-steady-state
approximation, once the system is sufficiently close to the steady
state.

## Abstract

In this paper we develop a new method for numerically approximating
sensitivities in parameter-dependent ordinary differential equations
(ODEs). Our approach, intended for situations where the standard
forward and adjoint sensitivity analysis become too computationally
costly for practical purposes, is based on the Peano-Baker series from
control theory. We give a representation, using this series, for the
sensitivity matrix `S` of an ODE system and use the representation to
construct a numerical method for approximating `S`. We prove that,
under standard regularity assumptions, the error of our method scales
as `O(h²)`, where `h` is the largest time step used when numerically
solving the ODE. We illustrate the performance of the method in
several numerical experiments, taken from both the systems biology
setting and more classical dynamical systems. The experiments show the
sought-after improvement in running time of our method compared to the
forward sensitivity approach. For example, in experiments involving a
random linear system, the forward approach requires roughly `sqrt(d)`
longer computational time, where `d` is the dimension of the parameter
space, than our proposed method.
