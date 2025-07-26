# Part 1 - Muscle models #
In this part of the workshop, you will run a biophysical muscle model based on cross-bridge dynamics for a simple muscle fiber model. 
You will see the effect of cross-bridge model parameters, protocols, activation on force development, and will compare it to phenomenological (Hill-type) muscle model.  

**Learning Objectives**
At the end of this tutorial, you should be able to:
- Simulate force output from biophysical muscle model under varying length, velocity, and activation conditions
- Understand the effect of cross-bridge properties on the steady-state and transient force generation properties of a muscle
- Compare and contrast force-generating properties of a biophysical and Hill-type muscle model

**Assignment 2.1: Knee pendulum test simulation demonstrates history-dependence of the biophysical model**
- Code: Part 2-Custom/kneeTest.m
- Run the code (F5 on Windows)
- Click on the force time series (middle row) to visualize the crossbridge distribution (bottom left) underlying the force at that moment in time. Observe the difference between isometric vs. pre-movement conditions. 
- Bonus: Change baseline activation and run the code again. Try pCa = 6.7, 6.9, 7.2, 8

**Assignment 2.2: Standing balance simulation with agonist-antagonist muscle pair**
- Code: Part 2-Custom/standingBalance.m
- Run the code (F5 on Windows)
- Bonus: Change baseline activation and run the code again.
