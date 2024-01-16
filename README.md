# Synthetic Data Generation Using the Robinson/Freyer/Hindriks Model

## Overview

This is a library written to implement the neural activity model described
in Robinson 2002, Freyer 2011, and Hindriks 2023. The model is intended
to simulate the cortico-thalamic loop and alpha-band oscillations that this
loop supports, as well as additional resonances within and between the
cortex and thalamus.

Key references are:

* __Robinson 2002__
(https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.041924) -
Describes a model of four interacting regions, finds steady-state operating
points, and analyzes perturbations around those points to find oscillation
modes.

* __Freyer 2011__ (https://www.jneurosci.org/content/31/17/6353.short) -
Describes an extension to the model where noise is modulated by the
network's output, providing a closer match to the distribution of
oscillation modes found in biological recordings.

* __Hindriks 2023__ (https://www.nature.com/articles/s42003-023-04648-x) -
Describes an extension to the model that couples together several instances
of the Robinson/Freyer model to simulate co-oscillation of differenet brain
regions.

## Folders

* `examples` -- Sample code that uses this library.
* `library` -- Matlab library function code.
* `manual` -- Project documentation PDFs.
* `manual-src` -- LaTeX files for rebuilding project documentation.
* `reference-code` -- Original code implemented by Hindriks and Tewarie for
the 2023 paper. Duplicated with permission.


_(This is the end of the file.)_
