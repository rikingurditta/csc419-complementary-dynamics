# Complementary Dynamics

## Introduction
Animating models can be very difficult task that that can require a lot of work, especially when trying to simulate real life physics in an artistic manner. When using programs to add physical effects to animations one issue is that the physics may undo the hand animation and the intent of the artist. Complementary dynamics is a system that is designed to enhance the animation while not interfering with the artists intentions. It achieves this by finding complementary displacements such that they are algebraically orthogonal to the animated displacements.
## Computing Complementary Displacement

### Overview

In short, we want to generate complementary displacements that are orthogonal to the rig displacements. To do so from a high level what would be used is newtons method to minimize the energy with constraints, where the energy is the Neohookean elasticity. These would affect the rig Jacobian which has calculations that would prevent allowing $\mathbf u^c_t$ from being able to undo $\mathbf u^r_t$.

### Minimizing Energy

In order to create interesting animations that would not undo the artists hand animation we minimize the potential energy of the system described, in our case, as the Neohookean elasticity. Given $\mathbf u_t^r, \mathbf J_t$ we want to find a a direction to move the function in such that it satisfies the orthogonal constraint so as to not undo the input displacements. We solve a system involving a block matrix

$$
\begin{bmatrix} \mathbf Q & \mathbf C^T\\ \mathbf C & 0 \end{bmatrix}\begin{bmatrix} \mathbf x \\ \lambda \end{bmatrix} = \begin{bmatrix} \mathbf l \\ 0 \end{bmatrix}
$$

where the values of $\mathbf Q, \mathbf l$ are determined due to solving the Lagrange multipliers at each iteration of the newtons method. The requirement that $\mathbf J_t^T \mathbf M \mathbf x = \mathbf 0$ is a constraint to ensure that the complementary effects would not undo the input rig animation. Once we solve this system, we can use $\mathbf x$ to calculate our search direction to find a suitable $\mathbf u_{t+h}^c$:

```c++
Eigen::VectorXd d = x - Uc;  // search in direction towards x
double alpha = 1.;  // search step size
while (f(Uc + alpha * d + Ur) > f(Uc + Ur) + 1e-8 * alpha * tmp_g.dot(d) and alpha > 1e-8) {
    alpha = alpha * 0.5;  // scaling factor is 0.5
}
Uc = Uc + alpha * d;  // update Uc
```


### Neohookean Elasticity

The algorithm is essentially minimizing physical energy subject to constraints that uphold the artist's intent. In our implementation, our energy function is Neohookean elasticity. The formula for this energy for a single tetrahedron is
$$
\Phi_{tet} = vol_{tet} \cdot \left( C \cdot (\det(\mathbf F)^{-2/3}\cdot \text{tr}(\mathbf F^T \mathbf F) - 3) + D \cdot (\det(\mathbf F) - 1)^2 \right)
$$
where $C$ and $D$ are physical constants that depend on the material being modelled. In our implementation, we use values of 107143 and 428571, based on the values given in a CSC417 assignment. $\mathbf F$ is the deformation gradient of the tetrahedron, a matrix which encodes in what way the tetrahedron has been deformed. To calculate the global energy, just sum up the energies for each tetrahedron.

To calculate the gradient $d\Phi/d\mathbf u^c$, we use the chain rule and combine $d\Phi/d\mathbf F$ and $d\mathbf F/d\mathbf u^c$, which are easier to compute (though both very large formulas). Here is an outline of how we calculate $d\Phi/d\mathbf u^c$ for a single tetrahedron:

```c++
// get deformation gradient
Eigen::Matrix3d F;
deformation_gradient(V, tet, U, F);
// calculate dF/dU, i.e. derivative of deformation gradient with respect to displacements
Eigen::Matrix912d B;
dF_dU_flattened(V, tet, B);
// calculate derivative of neohookean potential energy with respect to deformation gradient
Eigen::Vector9d dphi;
dphi_neo_hookean_dF(dphi, F, neohookean_C, neohookean_D);
// calculate derivative of neohookean potential energy with respect to displacements for current tet (chain rule)
Eigen::Vector12d g_tet = tet_volume * B.transpose() * dphi;
```

Just like with energy, we sum up the gradients for each tetrahedron to get the global gradient. We do this by distributing:

```c++
g.segment(tet(0) * 3, 3) += g_tet.segment(0, 3);
g.segment(tet(1) * 3, 3) += g_tet.segment(3, 3);
g.segment(tet(2) * 3, 3) += g_tet.segment(6, 3);
g.segment(tet(3) * 3, 3) += g_tet.segment(9, 3);
```

The hessian is calculated per-tetrahedron and distributed in a similar way, except it is distributed to a global sparse square matrix rather than a global dense vector.

### Rig Jacobian and Complementarity

The crucial part of the algorithm is the rig Jacobian $\mathbf J$. This is where the artist's animation comes into play.

Artists create animations using "rigs" that consist a (relatively small) number of parameters to control the mesh. In our case, we use linear blend skinning, which is where the artist places a number of "bones" around the mesh. Each bone influences the vertices around it by some weight (chosen by the artist), and for each vertex, the sum of bone weights adds up to 1. The artist can transform the bones, inducing a weighted average transformation on each vertex. In this setup, the parameters of the rig are affine transformations of the bones.

Given a matrix `W` of weights `W(i, j)` of bone `j` on vertex `i`, we can use `igl::lbs_matrix` to calculate the skinning matrix, so that if `T` is a matrix of affine transformations, then `W * T` gives us the deformed vertex positions.

```c++
Eigen::MatrixXd LBS;
igl::lbs_matrix(V, weights, LBS);
TU = LBS * T;
viewer.data().set_vertices(TU);
```

`T` parametrizes this deformation, so we define $\mathbf p$ to be a vector encoding of `T`. We can find the derivative of the skinning transformation, giving us $d\mathbf u^r/d\mathbf p = \mathbf J$, the rig Jacobian. The calculations for this step are outlined in equations 13 and 14 of Complementary Dynamics, but essentially boil down to computations of
$$
\mathbf J_{ij} = w_{ij} \mathbf I \otimes \begin{pmatrix}\mathbf v_i^T & 1\end{pmatrix}
$$
where $w_{ij}$ is the weight of the $j^{\text{th}}$ bone on the $i^{\text{th}}$ vertex, $\mathbf I$ is a 3 by 3 identity matrix, $\otimes$ is the [Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product), and $\mathbf v_i$ is the rest position of the $i^{\text{th}}$ vertex. This equation is implemented as

```c++
// v1 = [v^T 1]
Eigen::MatrixXd v1(1, 4);
v1.block<1, 3>(0, 0) = V.row(i);
v1(3) = 1;
// calculate kronecker product Id (*) v1
Eigen::MatrixXd kronecker = Eigen::MatrixXd::Zero(3, 12);
kronecker.block<1, 4>(0, 0) = v1;
kronecker.block<1, 4>(1, 4) = v1;
kronecker.block<1, 4>(2, 8) = v1;
// each block of J is w_ij * Id (*) [v^T 1]
J.block<3, 12>(i * 3, j * 12) = W(i, j) * kronecker;
```

Our derivations above are looking to optimize $\mathbf u^c$ so that $\mathbf J^T \mathbf M \mathbf u^c = \mathbf 0$. Intuitively, this ensures that $-\mathbf u^c$ does not correspond to any rig transformation, so $\mathbf u^c$ does not undo any of the artist's work.

## Conclusion
conclusion(1 paragraph at most, going over intuition)

