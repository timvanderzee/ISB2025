# ISB2025: Workshop "Biophysical muscle models for musculoskeletal simulation"
Authors: Lena Ting, Surabhi Simha, Hansol Ryu, Tim van der Zee, Friedl De Groote.

**Synopsis**

To date, there are extremely few musculoskeletal movement simulations that use biophysics-based muscle models. Most muscluloskeletal simulations use phenomenological models of muscle force generation, i.e. “Hill-type models” that compute muscle force from  muscle activation, length, and velocity based on empirical measurements and do not incorporate biophysical mechanisms of muscle force generation. Data used to generate Hill-type models are collected in controlled isometric or isotonic or isokinetic experiments and have proven useful to predict muscle force in steady-state conditions, particularly during muscle shortening. However, Hill-type models have poor fidelity in unsteady conditions where transient properties of muscle force are important, particularly when muscles are stretched and geenrate force eccentrically. Muscle forces in response to stretch can have a transient increase in force referred to as short-range stiffness, caused when attached muscle cross-bridges are stretched, and then rapidly detach. Such transient forces are critical to maintaining the length of the muscle in response to environmental perturbations. Previously, we augemnted Hill-type models with a phenomenological description of muscle short-range stiffness to enable more realistic simulations of perturbation responses to standing balance, and in clinical test of joint hyper-resistance, i.e. spasticity, where the lower limb is dropped under the force of gravity. However, these phenomenological short-range stiffness model is not generalizable to artibrary movement conditions as the short range stiffness is highly dependent on the prior movement and activation history of the muscle when the muscle is stretched. Therefore, we have developed biophysical muscle models based on muscle cross-bridge dynamics from which both force-velocity and short-range stiffness properties emerge. In this workshop, we present (1) a biophysical muscle model as well as (2) musculoskeletal simulations using the biophysical model using open-source software based on OpenSim. Hands-on tutorials focus on (1) exploring the force generating properties of the cross-bridge model parameters, stretch protocols and muscle activation in comparison to Hill-type models derived from the biophysical model, and (2) simulating a clinical test of joint hyper-resistance and perturbed standing balance.


**Dependencies**

Part 1 and Part 2 - Custom requires the following:
- MATLAB (license required)- can be downloaded from: www.mathworks.com

Part 2 - OpenSim code requires the following:
- MATLAB (license required)- can be downloaded from: www.mathworks.com
- OpenSim 4.5 (open/free) - can be downloaded from: https://simtk.org/projects/opensim/ </br>
  On Mac you may need to go to Settings > Privacy & Security > allow apps from unidentified developers </br>
  Be sure to run OpenSim to create your OpenSim working folder before moving on to the next step.
- MATLAB – OpenSim API (open/free) - follow the steps described here: https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53089380/Scripting+with+Matlab
- Casadi (open/free) – can be downloaded from: https://web.casadi.org/get/. Make sure that the Casadi is added to your MATLAB path. In MATLAB go to Home->Set Path->Add with Subfolders ... and navigate to and select the casadi folder (e.g. casadi-windows-matlabR2016a-v3.5.5) </br>
Note that casadi for Apple Silicon running R2023b version does not work, please use the 2018b version.

**Getting started**

- To get started, download or clone this repository.
- Next, click on a folder (e.g. Part 1, Part 2 - OpenSim) and follow the instructions shown in README.md

**Contact**

tim.vanderzee@kuleuven.be

Part 1: hansol.ryu@gatech.edu, surabhi.n.simha@emory.edu
Part 2 OpenSim: tim.vanderzee@kuleuven.be
Part 2 Custom: hansol.ryu@gatech.edu


