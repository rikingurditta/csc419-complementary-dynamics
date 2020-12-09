# Complementary Dynamics

## Introduction
Animating models can be very difficult task that that can require a lot of work, especially when trying to simulate real life physics in an artistic manner. When using programs to add physical effects to animations one issue is that the physics may undo the hand animation and the intent of the artist. Complementary dynamics is a system that is designed to enhance the animation while not interfering with the artists intentions. It achieves this by finding complementary displacements such that they are algebraically orthogonal to the animated dispalcements.
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

- above, we discussed minimizing the potential energy of the system
- in practice, potential energy defined by the type of material being modelled
- neohookean potential is for continuum solids made of elastic material
- gradient is force, hessian is "stiffness"

### Rig Jacobian and Complementarity

- rig works with linear blend skinning
- parameters of rig are affine transformations of "bones"
- rig function is weighted average of bone transformations
  - maybe have a picture of this
- rig jacobian is derivative of mapping from bone transformations to vertex displacements
- physical interpretation of $\mathbf J^T \mathbf M \mathbf u^c = \mathbf 0$:
  - modifying rig parameters, aka coming up with affine transformations of bones, could not produce any part of $\mathbf u^c$ (or $\mathbf u^c$)
  - thus $\mathbf u^c$ cannot be undoing any transformation that the artist produced

## Conclusion

conclusion(1 paragraph at most, going over intuition)

