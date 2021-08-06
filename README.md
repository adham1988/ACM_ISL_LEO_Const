# ACM_ISL_LEO_Const
To run the simulation, run sim_final.m
It takes about 3568.852995 seconds, 1 hour approximately.

The scripts simulate Adaptive Coding and Modulation between two LowEarth  Orbit  satellites  in  the following inter-satellite  linkcommunication:
-Walker-star intra-plan link.
-Walker-star inter-plan.
-Walker-star cross-seam link.
-Walker-delta inter-plan.
Implementation   of Extended Hamming for Forward Error Correc-tion,  BPSK,  QPSK  and  QAM  for  modulation,the  satellite  constellation  simulator,  a  Signal-Noise  Ratio  estimator,  as  well  as  an  AdaptiveCoding and Modulation controller.

Satelitte orbits and constellation are implemented through kep2cart.m where the positions of the satellites at given time are saved. This, to be later used to estimate the distance which is used in path loss calculation.  
Algorithm implemented in Matlab R2020b that calculates Cartesian state vectors from a set of Keplerian orbit elements.
The input to the algorithm is a set of Keplerian orbit elements:
The semi-major axis a,
eccentricity e,
argument of periapsis w,
longitude of ascending node omega,
inclination i,
and the mean anomaly M0.
Additionally the time t has to be specified.
The outputs of the algorithm are the Cartesian position and velocity vectors r and v.

some notations:
-the simulation is made of 291 step, from 0 upt to 290 with 10s increment each. At each step, 50,000 bits are transmitted from one satellite another one.
-At first iteration, a conservative choice of modulation and coding with Hamming(64,57) and BPSK. 
