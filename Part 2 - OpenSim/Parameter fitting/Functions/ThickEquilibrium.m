function [error, dX] = ThickEquilibrium(Q0, F, DRX, dDRXdt, k1, k2, kF, act)

[J1, J2] = ThickFilament_Dynamics(F, DRX, k1, k2, kF, act);

dX = J1 - J2;

error = dDRXdt - dX; 

end