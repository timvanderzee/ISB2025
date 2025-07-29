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

_______________________________________
**Acknowledgements**

This work was funded in part by an NIH R01 R01 HD90642 Multi-scale models of proprioceptive encoding to reveal mechanisms of impaired sensorimotor controland an NIH Supplement to Support Enhancement of Software Tools for Open Science (NIH R01 HD90642-S1)

**Citing this work**

If you find this work helpful in your research please cite this repository and the following papers (can be found in /literature) as relevant: 

The structure of a two- four-state cross-bridge model to reproduce history dependence in biophysical muscle spindle models was presented in:

(2-state) Blum, K.P., Horslen, B., Campbell K.S., Nardelli P., Cope T.C., Ting L.H.* (2020) Diverse and complex muscle spindle firing properties emerge from multiscale muscle mechanics.  Elife, Dec 28;9:e55177. doi: 10.7554/eLife.55177. ( https://doi.org/10.7554/eLife.55177)

(4-state) Simha, S.N., Ting, L.H. (2024) Intrafusal cross-bridge dynamics shape history-dependent muscle spindle responses to stretch. Journal of Experimental Physiology, Jan;109(1):112-124. doi: 10.1113/EP090767. (https://doi.org/10.1113/EP090767)
 
The parameters of the cross-bridge model presented here have been optimized to fit history-dependence muscle fiber data described in:
Horslen, B.C., Milburn, G.N., Simha, S.N., Blum K.P., Campbell, K.S., Ting, L.H. (2023) History-dependent muscle resistance to stretch remains high after small, posturally-relevant pre-movements. Journal of Experimental Biology, Sept 15;226(18):jeb245456. doi: 10.1242/jeb.245456. (https://doi.org/10.1242/jeb.245456)

We also present a version of the four-state model using a Gaussian approximation of the cross-bridge distribution as described in : 
Zahalak, G. I.  (1981).   A   distribution-moment   approximation   for   kinetic   theories   of   muscular
contraction. Math. Biosci. 55, 89–114. (https://doi.org/10.1016/0025-5564(81)90014-6)

The perturbed standing balance models are based on:
De Groote, F., Allen, J. L., Ting, L. H.* (2017) Contribution of muscle short-range stiffness to initial changes in joint kinetics and kinematics during perturbations to standing balance: A simulation study. Journal of Biomechanics, Apr 11;55, 71-77. (https://doi.org/10.1016/j.jbiomech.2017.02.008)

The pendulum test models of joint hyper-resistance (spasiticity) in cerebral palsy are based on:
De Groote, F., Blum, K.P., Horslen, B.C., Ting, L.H.* (2018) Interaction between muscle tone, short-range stiffness and increased sensory feedback gains explains key kinematic features of the pendulum test in spastic cerebral palsy: A simulation study. PLoS ONE, Oct 18;13(10):e0205763. (https://doi.org/10.1371/journal.pone.0205763)

Willaert, J., Desloovere, K., Van Campenhout, A.C.., Ting, L.H., De Groote, F. (2020) Movement history influences pendulum test kinematics in children with spastic cerebral palsy. Frontiers in Bioengineering and Biotechnology, Aug 7;8:920. (https://doi.org/10.3389/fbioe.2020.00920)

Willaert, J., Desloovere, K., Van Campenhout, A., Ting L.H., De Groote, F. (2024) Identification of neural and non-neural origins of joint hyper-resistance based on a novel neuromechanical model. IEEE Trans Neural Syst Rehabil Eng. 2024:32:1435-1444. (http://doi.org/10.1109/TNSRE.2024.3381739)



