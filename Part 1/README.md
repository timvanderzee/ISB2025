# Part 1 - Muscle models #
In this part of the workshop, you will run a biophysical muscle model based on cross-bridge dynamics for a simple muscle fiber model. 
You will see the effect of cross-bridge model parameters, length changes, and activation on force development. You will then compare the force response from the biophysical muscle model to the force response from a phenomenological (Hill-type) muscle model that is derived from the biophysical model.  

The model consists of:

**Learning Objectives**
At the end of this tutorial, you should be able to:
- Simulate force output from biophysical muscle model under varying length, velocity, and activation conditions
- Understand the effect of cross-bridge parameters on the steady-state and transient force generation properties of a muscle
- Compare and contrast force-generating properties of a biophysical and Hill-type muscle model

**Assignment 1.0: Preparation**
  -	Go to ISB2025\Part 1\tutorial and open the script called ‘GUI_XB_n_Hill.m’
  -	Run the code (F5 on Windows) to open the GUI windon

<img width="1718" height="892" alt="GUIscreenshot1" src="https://github.com/user-attachments/assets/ce7c5514-7f83-4317-a790-bb5dd40c4078" />

The panel in the upper right allows you to set the values for parameters that change the rate functions for crossbridge attachment (f,w) and detachment (k11, k12, k21, k22). Although the model also contains states for cooperativity, these parameters cannot be modified within the GUI.

If you would like to simulate the biophysical muscle model outside of the GUI, please see the Alternate Preparation section at the end of the this document.


**Assignment 1.1: Effect of cross-bridge rate functions on cross-bridge distribution and force development (at fixed velocity and activation)**
   While running the GUI (‘GUI_XB_n_Hill.m’), 
  - Click on the force time series (bottom right) to visualize the crossbridge distribution (bottom left) underlying the force at that moment in time.
  - Change the attachment and detachment rate parameters and calcium concentration (activation) by clicking on the upper left plot. Re-run the simulation to see how crossbridge distribution and force development are affected.

**Assignment 1.2: Effect of muscle length change on cross-bridge distributions and force generation (at fixed velocity and activation)**
   While running the GUI (‘GUI_XB_n_Hill.m’), 
  - Select the desired attachment and detachment rate parameters and calcium concentration.
  - Click on the desired muscle stretch pattern (ramp-and-hold, stretch-shorten cycles).
  - Run the protocol to observe movement history effects on crossbridge distribution and force generation. 

**Assignment 1.3: Simulations characterizing steady-state properties of the muscle fiber: force-velocity**
   While running the GUI (‘GUI_XB_n_Hill.m’), 
  - Run the protocol to simulate stretches to generate a force-length and force-velocity curve based on your XB model parameters. (click "Generate Hill") 
  - Compare cross-bridge and Hill-type forces.
  - Try different protocols and compare forces. 

**Alternative 1.0: Preparation**
If you like, you can simulate these muscle models without GUI and edit the code directly. 
  - *XB_n_Hill.m* allows you to set up attachment/detachment parameters, change activation and simulation protocols. 
  - *FL_FV_FpCa.m* allows you to set up attachment/detachment parameters and simulate stretches to generate a force-length, force-velocity and force-pCa curves. You can save these curves and load into *XB_n_Hill.m* to simulate muscle force based on Hill-type muscle model you generated.  
  
