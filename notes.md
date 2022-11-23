# Complementary Dynamics Notes

$$
\newcommand{\ds}{\displaystyle}
\newcommand{\x}{\mathbf x}
\newcommand{\u}{\mathbf u}
\newcommand{\p}{\mathbf p}
\newcommand{\M}{\mathbf M}
\newcommand{\J}{\mathbf J}
\newcommand{\f}{\mathbf f}
\newcommand{\norm}[1]{\lVert #1 \rVert}
\DeclareMathOperator*{\argmin}{argmin}
$$

## Notation

- $n$ is the num of vertices in the mesh
- $d$ is dimension
- $\u^r : \R^m \to \R^{dn}$ is rig displacement
  - $m$ is number of parameters
    - parameters can be written as $\p \in \R^m$
  - rig maps parameters (e.g. bone positions) to mesh vertex displacements
  - $\u^r_t = \u^r(\p_t)$, i.e. *rig displacement* given parameters for time $t$
  - $\ds \J_t = \frac{\partial \u^r}{\partial \p}\Bigg\vert_{\p_t} \in \R^{dn \times m}$ is the rig jacobian, i.e. how the rig is changing with respect to the parameters at time $t$
- $\u^c \in \R^{dn}$ is complementary displacement
  - what we want to solve for
- $\u_t = \u^r_t + \u^c_t \in \R^{dn}$ is the final mesh displacement
  - i.e. final mesh *position* is $\mathbf V + \u_t$ where $\mathbf V \in \R^{dn}$ is the rest vertex positions of the mesh

- $h$ is step size
- $E_t$ is energy to minimize at each time step
- $\Phi$ is potential energy
- $\M \in \R^{dn \times dn}$ is mass matrix
- $\f : \R^{dn} \to \R^{dn}$ is external force at each time
- dots denote finite differences over time steps
  - $\ds \dot\u = \frac{\u_t - \u_{t-h}}{h}$
  - $\ds \ddot\u = \frac{\dot\u_t - \dot\u_{t-h}}{h}$

- $\norm{\x}^2_M = \x^T \M \x$ for any $x \in \R^{dn}$

## Goal

Given $\p_t$ for every $t$, we want to solve for all complementary displacements $\u^c_t$  so that

- $\u^c$ reacts to all forces, both internal and external
  - (physics based animation)
- $\u^c$ does not *undo* rig displacements $\u^r$

## Energy to minimize

Let's set up our goal as a constrained optimization problem for every time $t$:
$$
\u_t = \argmin_{\u_t} E_t(\u_t)
$$
We can break this into the different physically significant energies we want to consider:

- $\Phi(\u_t)$ is the potential energy
  - energy function not time dependent because it's a physical calculation
- $\ds \frac{h^2}{2} \ddot\u_t^T \M \ddot\u_t$ is momentum term
- $-\u_t^T \f(\u_t)$ is external work

So what we really want to minimize is
$$
\u_t = \argmin_{\u_t} E_t(\u_t) =  \argmin_{\u_t}\  \Phi(\u_t) + \frac{h^2}{2}\ddot\u_t^T \M \ddot\u_t - \u_t^T \f(\u_t)
$$
However, we only want to find complementary displacements, we already know what our rig displacements are, so really we want to solve
$$
\u^c_t = \argmin_{\u^c_t} E_t(\u^r + \u^c)
$$
treating $\u^r$ as a constant

We need constraints, otherwise the solution would be $\u_c = -\u_r$ which is not useful

We don't want $\u^c$ to interfere with the rig, or equivalently, we don't want to be able to reach a displacement using $\u^c$ that we might reach just by tweaking the parameters. Thus, the closest we should be able to get to $\u_t$ by tweaking the parameters is by using our already set $\p_t$.

Writing this as an equation, what we want is
$$
\begin{align*}
\p_t &= \argmin_{\p} \frac{1}{2} \norm{\u_t - \u^r(\p)}^2_M \\
\p_t &= \argmin_{\p} \frac{1}{2} \norm{(\u^r(\p_t) + \u^c) - \u^r(\p)}^2_M
\end{align*}
$$
Let's take the derivative of both sides:
$$
\begin{align*}
0 &= \frac{\partial}{\partial \p} \left( \frac{1}{2} \norm{(\u^r(\p_t) + \u^c) - \u^r(\p)}^2_M \right) \\
0 &= \left( \frac{\partial}{\partial \p} (\u^r(\p_t) + \u^c - \u^r(\p)) \right)^T \M (\u^r(\p_t) - \u^r(\p) + \u^c) \\
0 &= \left( \frac{\partial \u^r}{\partial \p} \right)^T \M (\u^r(\p_t) - \u^r(\p) + \u^c)
\end{align*}
$$
This is $0$ at our critical point $\p_t$, so let's evaluate there:
$$
\begin{align*}
0 &= \left( \frac{\partial \u^r}{\partial \p}\Bigg\vert_{\p_t} \right)^T \M (\u^r(\p_t) - \u^r(\p_t) + \u^c) \\
0 &= \left( \frac{\partial \u^r}{\partial \p}\Bigg\vert_{\p_t} \right)^T \M \u^c
\end{align*}
$$
If we define our *rig jacobian* $\ds \J_t = \frac{\partial \u^r}{\partial \p}\Bigg\vert_{\p_t}$, then we can rewrite this as
$$
0 = \J_t^T \M \u^c
$$
This means that $\u^c$ must be "orthogonal" to the way the rig is changing with respect to the parameters at time $t$.

So in the end, our constrained optimization is
$$
\u^c_t = \argmin_{\u^c_t} E_t(\u^r + \u^c) \ \text{ subject to } \ \J_t^T \M \u^c = 0
$$
