---
sidebar_position: 1
---

# Boltzmann Transport Theory

The Boltzmann Transport Equation (BTE) is given by
$$
\partial_t n_i + \boldsymbol{v}_i\cdot\nabla_{\boldsymbol{r}} + \boldsymbol{F}\cdot\nabla_{\boldsymbol{k}} n_i + \Gamma_i = 0
$$
where $n_i(\boldsymbol{r})$ is the distribution function of carriers in combined state index $i$ (encompassing the wavevector $\boldsymbol{k}$ and band $\nu$), $\boldsymbol{v}_i$ is the group velocity of state $i$, $\boldsymbol{F}$ is an external driving force, and $\Gamma_i[\{n_j(\boldsymbol{r})\}]$ is the collision operator, whichspecifiestherateat which carriers scatter into and out of state $i$ as a function of the full carrier distribution $n_j$ at position $\boldsymbol{r}$.

SpaRTaNS solves the BTE with a few simplifying assumptions:

1. Steady state $\partial_t n_i = 0$.
2. Linearized collisions: write $\delta n_i = n_i - n_{0,i}$, where $n_{0,i}$ is an equilibrium distribution, then expand $\Gamma_i$ to linear order in $\delta n_i$.

These allow us to write the BTE as
$$
\left(\frac{\partial \Gamma_i}{\partial n_j} + \delta_{ij} \boldsymbol{v}_i\cdot\nabla_{\boldsymbol{r}}\right)\delta n_j = S_i,
$$
where summation is implied over repeated indices and we have written the forcing term as a *source* of carriers $S_i = \boldsymbol{F}\cdot\nabla_{\boldsymbol{k}} n_i$.

Next, we separate the collision operator into diagonal terms, representing decay with lifetime $\tau_i$, and off-diagonal ‘mixing’ terms:
$$
\frac{\partial \Gamma_i}{\partial n_j} = \tau_i^{-1} \delta_{ij} - M_{ij}.
$$
Using this decomposition, we write the BTE as
$$
\left(1 + \tau_i\boldsymbol{v}_i\cdot\nabla_{\boldsymbol{r}}\right)\delta n_j = \tau_i S_i + G_{ij} \delta n_j,
$$
where summation is again implied over repeated indices. Here 
$$
G_{ij} \equiv \tau_i M_{ij}
$$
with no summation implied.

SpaRTaNS solves this equation iteratively, by expressing $\delta n_j$ as a power series in $G$:
$$
\delta n_j = \delta n_j^0 + \delta n_j^1 + ...
$$
 where
$$
\begin{aligned}
(1+\tau_i \boldsymbol{v}_i\cdot\nabla_{\boldsymbol{r}})\delta n_i^0 &= \tau_i S_i \\
(1+\tau_i \boldsymbol{v}_i\cdot\nabla_{\boldsymbol{r}})\delta n_i^k &= G_{ij} \delta n_j^{k-1}
\end{aligned}
$$
This approach converges so long as the spectral radius of $G_{ij}$ is less than unity. Otherwise more sophisticated approaches like [Jacobi weighting](https://arxiv.org/abs/1811.01059) or alternate decompositions (e.g. choosing artificially smaller $\tau_{ij}^{-1}$, so that $M_{ij}$ has non-zero diagonal entries) must be used.

### References

[Varnavides+2018](https://arxiv.org/abs/1811.01059)

[Jermyn+2017](https://arxiv.org/abs/1707.07060)
