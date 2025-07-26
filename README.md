# ISB2025: Workshop "Biophysical muscle models for musculoskeletal simulation"
Authors: Lena Ting, Surabhi Simha, Hansol Ryu, Tim van der Zee, Friedl De Groote.

**Synopsis**

To date, there are extremely few musculoskeletal movement simulations that use biophysics-based muscle models. Most simulations use phenomenological models of muscle force generation, i.e. “Hill-type models” that provide black-box input-output relationships between muscle activation, length, and force, offering no insight into biological processes of muscle. Because data underlying Hill-type models is collected in controlled isometric or isotonic experiments, Hill-type models have poor fidelity in behaviourally-relevant movements in which muscles’ operating conditions are very different from those in the experiments that underlie the model. For example, we previously needed to extend Hill-type models with a phenomenological description of short-range stiffness to capture the response of the muscle to stretch in clinical tests and during standing balance. Recently, we have developed biophysical muscle models based on muscle crossbridge dynamics from which both force-velocity and short-range stiffness properties emerge. In this workshop, we present (1) these biophysical muscle models as well as (2) musculoskeletal simulations using these biophysical models using open-source software based on OpenSim. Hands-on tutorials focus on simulating a clinical test of joint hyper-resistance and perturbed standing balance.


**Dependencies**

This code requires the following:
- MATLAB (license required)- can be downloaded from: www.mathworks.com
- OpenSim 4.5 (open/free) - can be downloaded from: https://simtk.org/
- MATLAB – OpenSim API (open/free) - follow the steps described here: https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53089380/Scripting+with+Matlab
- Casadi (open/free) – can be downloaded from: https://web.casadi.org/get/. Make sure that the Casadi is added to your MATLAB path. In MATLAB go to Home->Set Path->Add with Subfolders ... and navigate to and select the casadi folder (e.g. casadi-windows-matlabR2016a-v3.5.5)

**Getting started**

- To get started, download or clone this repository.
- Next, click on a folder (e.g. Part 1, Part 2 - OpenSim) and follow the instructions shown there

**Contact**

tim.vanderzee@kuleuven.be
