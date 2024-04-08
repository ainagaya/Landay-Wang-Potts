4. Landau-Wang simulation of the 2D Potts model.
References:
1.F. Wang and D.P.Landau, PHYSICAL REVIEW E, VOLUME 64, 056101 (2001).
2. D.P. Landau, S-H. Tsai, M. Exler, Am. J. Phys. 72 (2004) for pedagogical discussion.
First submission: points 1 and 2 below.
1. Implement the Landau-Wang algorithm for the 2D q-state Potts model. The algorithm is
described in Project 3 above.
Input: q, L, Emin, Emax, finitial, ff inal, niter (for example Landau and Wang use finitial = e,
niter = 104
). Output: energy histogram and final n(E).
Suggestion: set Emin to the ground state energy, and Emax = −Emin (what configuration
gives this energy?) for small L. For large L, Emax can be reduced to Emax = 0.2 in order to
improve statistics, since positive values of the energy do not contribute significantly to the
canonical distribution at reasonable temperatures anyway.
2. To test the program, run it with small L (e.g. L = 20) and q = 2, which corresponds to the
Ising model, and then proceed as in point 2 of Project 3 (note that n(Egroundstate) = q.
3. Perform simulations with q = 10, L = 40 (larger L can be tried, but it will probably will
require to reduce the energy interval). Plot the microcanonical entropy, compute the probability distribution of the energy density E/N and plot it for different β, and determine βc
as the value at which the two peaks of the distribution have the same height.
4. Plot the the internal energy, the entropy, free energy, and specific heat (from hE
2
i − hEi
2
),
as a function of the temperature. You should get results similar to Fig.2 the paper by Wang
and Landau.
5. (Optional) Test the effect of changing the size L.
