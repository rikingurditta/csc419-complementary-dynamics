# Complementary Dynamics

## Introduction

- motivation
  - animating is a lot of work
  - physics based animation is hard to control artistically
  - usually hard to combine them because physics might undo hand animation
- complementary dynamics is a system for physics that doesn't interfere with animation
  - does this by finding "complementary" displacements based on physics
  - these are displacements that are algebraically orthogonal to animated displacements

## Computing Complementary Displacement

### Overview?

general overview of the algorithm(1-2 paragraphs maybe? Still should be fairly short)

### Minimizing energy

- newton's method stuff
- what exactly are we minimizing and how? this is the core of the algorithm
- newtons method(motivate but mostly go over the derivatives then explain the steps for newtons method, maybe 2-3 paragraphs. Should have a section on line search part i guess. I think important to talk about what this is analagous to and why it is doing what it does)
- talk about the derivatives and solving newtons method(1-2 paragraphs I don't think this needs to be that in depth)


### Neohookean Elasticity

The algorithm is essentially minimizing physical energy subject to constraints that uphold the artist's intent. In our implementation, our energy function is Neohookean elasticity. The formula for this energy for a single tetrahedron is
$$
\Phi_{tet} = vol_{tet} \cdot \left( C \cdot (\det(\mathbf F)^{-2/3}\cdot \text{tr}(\mathbf F^T \mathbf F) - 3) + D \cdot (\det(\mathbf F) - 1)^2 \right)
$$
where $C$ and $D$ are physical constants that depend on the material being modelled. In our implementation, we use values of 107143 and 428571, based on the values given in a CSC417 assignment.

- above, we discussed minimizing the potential energy of the system
- in practice, potential energy defined by the type of material being modelled
- neohookean potential is for continuum solids made of elastic material
- gradient is force, hessian is "stiffness"

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
// calculate kronecker product of identity (*) v1
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

