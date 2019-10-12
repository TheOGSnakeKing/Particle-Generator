# Physcis based simulation of particles.

Simulates upto 30000 particles at a stable frame rate. The particles have initial attributes including position,velocity, lifespan and color being drawn from uniform and guassian probability distribution.

Included collisions. For demonstration purposes, the particles collide with one plane in the program and not with the other. Uses barycentric co-ordinates for collision.

Includeds a circular vortex force to the upper right side of the particle simulator for interesting particle choreography.

Optimized memory allocation by storing random variables by guassian and normal distribution and storing them in a data structure. The program refers to these whenever a random number is required and performs operations to get the number needed.

Uses parallel arrays to dynamically store and delete the location and other parameters of particles.

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Needs freeGlut and openGl to run.

Copyright reserved © Nagaraj Raparthi.

