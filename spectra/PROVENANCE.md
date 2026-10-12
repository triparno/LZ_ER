# Provenance of the beta spectra

These CSV files are DIGITISATIONS of published figures from

  S. J. Haselschwardt, J. Kostensalo, X. Mougeot and J. Suhonen,
  "Improved calculations of beta decay backgrounds to new physics in
  liquid xenon detectors", Phys. Rev. C 102, 065501 (2020), arXiv:2007.13686.

They are NOT the original authors' tabulated data and must not be cited as
such. Anyone needing the underlying calculation should contact those authors.

Curve digitised: the g_A = 0.85 + Exchange solid red curve, and its shaded
uncertainty band where resolvable. The dashed "No Exchange" curve, which lies
inside the band, is treated as band interior and never as an edge.

Extraction is by colour segmentation; trace.py, trace212.py and trace85.py are
the scripts. Each takes the source image and output CSV as arguments:

    python trace.py [source.png] [output.csv]

The source images are crops of the published figures and are NOT redistributed
here. Make them by cropping the relevant panel of the figure in the cited paper
to the plot frame; the pixel calibration constants at the top of each script
assume that framing and must be re-measured for a different crop.

The overlay images trace_check.png, trace212_check.png and trace85_check.png
in this directory were produced during the original digitisation and show the
traced curves on the source figures. They are the check on the extraction; the
scripts themselves do not regenerate them.

    pb214_gA085_exchange.csv   2 keV grid, 0-1000 keV   (Q = 1019 keV)
    pb212_gA085_exchange.csv   1 keV grid, 0-559 keV    (Q = 574 keV)
    kr85_exchange_1stforb.csv  1 keV grid, 0-685 keV    (Q = 687 keV)

Columns: E_keV, intensity, and for the lead isotopes intensity_upper,
intensity_lower, band_resolved. Where band_resolved = 0 the band is narrower
than the drawn line and the bounds are the line's own extent -- an upper limit
on the width, not a measurement.

Intensities are arbitrary units. Only ratios enter the paper, in which the
normalisation cancels. The lead curves stop short of their endpoints, so
extrapolate to zero if an absolute normalisation is needed.
