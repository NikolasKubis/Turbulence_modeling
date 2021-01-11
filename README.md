# Turbulence_modeling

This project provides algorithms for the problem of turbulence modeling in adaptive optics.

# Introduction

High resolution imaging in the visible spectrum from ground-based telescopes is seriously hampered byatmospheric turbulence. Adaptive optics (AO) is the system that corrects in real-time atmosphericaberrations. Incoming light is split in two beams and directed toward awavefront sensor (WFS) and the instrument camera. The information provided by the sensors isprocessed to estimate a future wavefront and converted into voltage for the actuators located undera deformable mirror (DM). The DM  flattens the wavefront in order to retrieve the original imagequality as if there was no atmosphere.

# Methods:

 - Static wavefront reconstruction
 - Vector Auto-Regressive model of Order 1
 - Subspace Identification

# Estimation

After a linear state-space model is established a Kalman filter is utilized to estimate the system's state and to further control the deformable mirror.