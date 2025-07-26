# Part 1 - Muscle models #
In this part of the workshop, you will run a biophysical muscle model that is based on cross-bridge dynamics to simulate force from a muscle fiber.

## **Learning Objectives**  
At the end of this tutorial, you should be able to:  
- Simulate force output from biophysical muscle model under varying length, velocity, and activation conditions  
- Understand the effect of cross-bridge parameters on the steady-state and transient force generation properties of a muscle  
- Compare and contrast force-generating properties of a biophysical and Hill-type muscle model

## Muscle Model
The model consists of:  
1) a contractile element that simulates muscle active force from crossbridge dynamics  
2) elastic elements, in series and in parallel with the contractile element, that simulate muscle passive force

### **Contractile element:**  
![Alt text](images/xbridgeModel.png)

1. Concentration of calcium ions (pCa) activate actin sites and make them available for crossbridge binding. Attached corssbirdges can further increase the available actin sites through the process of coopertativity (Campbell et al. 2014).  
2. Crossbridges cycle between detached and attached states governed by attachment rate f(∆x) and detachement rate g(∆x), where ∆x is the length of crossbridge relative to its resting length.  
3. Myosin heads (equvivalent to detached crossbridges) can also enter and leave a 'super-relaxed' state where they cannot attach to actin to form crossbridges. This is govered by entering rate function f<sub>SRX</sub> (muscle force) and leaving rate function g <sub> SRX </sub> (muscle force), where high muscle force can recruit more myosins from super relaxed state to be available for crossbridge attachment.

## **Assignment 1.0: Start muscle model**
1. Go to `ISB2025\Part 1\tutorial` and open the script called `GUI_XB_n_Hill.m`
2. Run the code (F5 on Windows, F5 or Fn+F5 on mac). This will open a GUI as shown below
3. Click at different time points on the force time series (bottom right) to see what the crossbridge distribution (bottom left) is at any given time.
4. Can you intuit how the crossbridge distribution at a given time gives the force (blue line, bottom right) at that time?

![Alt text](images/first_gui_view.png)

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
  
