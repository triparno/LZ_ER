# The low-energy edge of the WS RoI acceptance

The working behind the figures Sec. 4 quotes: the floor is 0.62 keVee, the
width 19.40 keVee, and the three treatments span 2.2%, which is contained in
the 20% assigned to dN/dE.

## Why the yield parametrisation cannot give the floor

Charge yield is bounded above by the case in which no ion recombines,

    Qy_max = 1 / [ W (1 + Nex/Ni) ] = 62.6 e-/keV     (W = 13.532 eV, Nex/Ni = 0.1812)

The tuned ER model tends to C*m2 = 85.3 e-/keV as E -> 0, i.e. m2 = 84.9
before the density correction C = 1.0051. That exceeds the ceiling by 36%.
The two cross at E = 0.97 keVee: below about 1 keVee the parametrisation
returns more electrons than there are ions, and is being used outside the
range over which it was tuned.

    E [keV]   Qy(E)   Ne = Qy*E   Ni = E/[W(1+r)]   Ne > Ni?
      0.193   84.60       16.32            12.07      YES
      0.26    82.9        21.6             16.3       YES
      0.50    78.0        39.0             31.3       YES
      1.0     70.6        70.6             62.6       YES (marginal)
      5.0     44.1       220             313          no
     20.0     20.5       409             1251         no
     70.0     17.7      1237            4379          no

## The floor from all three cuts

The WS RoI is three cuts, not one:

    S1c > 3 phd      -> Nph > 27.2
    S2c > 10^2.75    -> Ne  > 16.30
    S2  > 645 phd    -> Ne  > 18.70     (binding: stronger than the S2c floor)

Jointly, Nq > 27.2 + 18.7 = 45.9, so E > 45.9 * W = 0.621 keVee.

At that Nq: Ni = 38.9, Nex = 7.0, and the two requirements can be met only
for r in a knife-edge interval around 0.52. The capped model gives full
ionisation there (r = 0), hence Nph = Nex = 7.0 and S1c = 0.78 phd, which
fails the S1c cut. So sub-keV electron recoils enter the RoI only on an
upward excursion in recombination, and the 0.62 keVee figure is a boundary
of the acceptance rather than a threshold that typical events cross.

## The three treatments

    treatment                        floor [keVee]   width [keVee]
    S2c cut alone, uncapped yield        0.193           19.83
    S2c cut alone, full ionisation       0.261           19.76
    all three cuts jointly               0.621           19.40

Ceiling in all three cases: Qy(E)E = 409.5 gives E = 20.02 keVee, Qy = 20.46.

Spread in dN/dE across the three: 2.2%. The first is listed only to show
what the parametrisation returns if its sub-keV extrapolation is taken at
face value; discarding it leaves 1.9%.
