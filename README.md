HBT_correlation_
===============

### This code computes particle pair correlation function and HBT radii from the freeze-out hypersurface with Cooper-fyre emission function

The correlation function is computed according to,

C(q, K) = 1 + \frac{|\int d^4 x exp[i q \cdot x] s(x, K)|^2}{|\int d^4 x s(x, K)|^2},

where s(x, K) is the particle emission function based on Cooper-Fyre formula.

Since momentum q has 3 dimensions, this program will calculate the correlation function either on regular lattice points or at randomly sampled points.

It allows the users to compute both azimuthal averaged and azimuthal dependent correlation function (as function of \vec{KT}). Both KT-differential and KT-integrated correlation function.

