---
sidebar_position: 1
---

# Boltzmann Transport Equation

:::caution

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
\left(\frac{\partial \Gamma_i}{\partial n_j} + \delta_{ij} \boldsymbol{v}_i\cdot\nabla_{\boldsymbol{r}}\right)\delta n_j = S_i
$$
where we have written the forcing term as a *source* of carriers $S_i = \boldsymbol{F}\cdot\nabla_{\boldsymbol{k}} n_i$.

Under Construction.

:::

