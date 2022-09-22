# cyclostationarity-uac

This repository contains the Matlab code to accompany the following publication:
F.-X. Socheleau, *Cyclostationarity of Communication Signals in Underwater Acoustic Channels*,  accepted for publication in IEEE Journal of Oceanic Engineering. Note that the real data presented in the paper are not provided (data owned by third parties).

This code is published under the MIT License (see full license in the 'LICENSE' file)

## Prerequisites and installation
Required software: Matlab (tested on R2018b) with the *signal processing* and *statistics* toolboxes.

## Content
The __code__ folder contains the main scripts to run:
- ``demo1_MSML_Doppler.m`` is an illustration of the time-varying Doppler scale estimator described in Sec. III-C of the paper. It is applied to a QPSK signal filtered by a simulated multiscale-multilag channel in additive white Gaussian noise.
- ``demo2_MSML_CS_signature_detection.m`` is an illustration of the cyclostationary-signature detectors described in Sec. III-D. It is applied to a CP-OFDM signal (with a RRC window) filtered by a simulated multiscale-multilag channel in additive white alpha-stable noise.
- ``demo3_DISP_DSSS_symb.m`` is an illustration of the symbol-rate estimator described in Sec. IV-C. It is applied to a DSSS signal filtered by a low-frequency dispersive channel in additive white Gaussian noise.
- ``run_all_demos.m`` runs the three demo scripts in sequence

Some parameters can be modified such as SNR or the number of hydrophones in the third script. However, note that this is a demonstration code only designed to facilitate the understanding of the paper. Application to any other dataset may require modifications (optimization strategies, thresholds, etc.) for which the author will not provide any support.

The __subroutines__ folder contains all the functions called by the main scripts and the __data__ folder contains .mat files with signals and simulation parameters. Figures are stored in the results __folder__.
