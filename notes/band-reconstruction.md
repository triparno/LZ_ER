# Reconstructing the published ER band: an attempt that did not conclude

Fig. 2 of the Letter draws the ER band as "the calibrated NEST model, taking
flat energy spectra and evaluating along lines of fixed S1c". We implemented
that construction (band.py) to ask how much the distribution families of the
paper deviate from the published band, rather than from the response at the
event-conditioned energy. The attempt did not reach a usable answer and the
paper uses the fixed-energy quantity instead.

## What was verified against NEST.cpp at 541b643

 - the variance:   Variance = lambdaChen*r*(1-r)*Ni + omega^2*Ni^2   (line 359)
 - the skew-normal placement:
       widthCorrection = sqrt(1 - (2/pi) delta^2)
       scale    = sqrt(Variance)/widthCorrection
       location = (1-r)Ni - scale*delta*sqrt(2/pi)                   (lines 401-412)
   which is exactly Eqs. (5a)-(5b) of the paper.
 - omega extraction from QuantaResult::Variance.

So the response model as implemented in the paper matches the source.

## What did not reproduce

Our reconstructed band is wider than the published one:

    S1c    E[keVee]   published 10-90 [dex]   ours    ratio
    100       17            0.134            0.148    0.91
    300       48            0.172            0.214    0.81
    540       90            0.246            0.301    0.82

The median agrees well (within 0.011 dex everywhere), and two independent
digitisations of Fig. 2 agree with each other to ~0.005 dex, so the
disagreement is in the width, not the position or the trace.

## Why we did not report it

At least three candidate causes sit inside our own assumptions and we cannot
distinguish them from public information:

 1. The range of the "flat energy spectra" is not stated. Truncating at
    ~100 keVee reproduces the published width at S1c = 540 almost exactly
    (0.250 against 0.246) but changes nothing below S1c = 400.
 2. A residual ~10% disagreement persists at low S1c, which the truncation
    does not touch.
 3. The published band is symmetric about its median; ours is skewed, with
    the longer tail above, as a positive gamma_ER requires. At S1c = 540 the
    published (90-50) and (50-10) are 0.125 and 0.121; ours are 0.158 and
    0.138. This is a difference in the skewness treatment, not a width scale.

Effects we omit (Fano smearing of Nq, the binomial draw of Ni) would widen
our band further, so they cannot account for the gap -- but (1) narrows it
and is unknown to us, so no claim can be made in either direction.

## Files

    band.py                 the construction (collapsed, exact for flat spectra)
    band_check.py           compares against band_published.csv
    trace_band.py           our digitisation of Fig. 2
    band_published.csv      digitised 10/50/90 percentile curves
