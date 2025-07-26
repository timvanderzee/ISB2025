# Part 1 - Muscle models #
In this part of the workshop, you will run a biophysical muscle model based on cross-bridge dynamics for a simple muscle fiber model. 
You will see the effect of cross-bridge model parameters, protocols, activation on force development, and will compare it to phenomenological (Hill-type) muscle model.  

**Learning Objectives**
At the end of this tutorial, you should be able to:
- Simulate force output from biophysical muscle model under varying length, velocity, and activation conditions
- Understand the effect of cross-bridge properties on the steady-state and transient force generation properties of a muscle
- Compare and contrast force-generating properties of a biophysical and Hill-type muscle model

**Assignment 1.0: Preparation**
  -	Go to ISB2025\Part 1\tutorial and open the script called ‘GUI_XB_n_Hill.m’
  -	Run the code (F5 on Windows)

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

**Alternative way**
If you like, you can simulate these muscle models without GUI and edit the code directly. 
  - *XB_n_Hill.m* allows you to set up attachment/detachment parameters, change activation and simulation protocols. 
  - *FL_FV_FpCa.m* allows you to set up attachment/detachment parameters and simulate stretches to generate a force-length, force-velocity and force-pCa curves. You can save these curves and load into *XB_n_Hill.m* to simulate muscle force based on Hill-type muscle model you generated.  
  
