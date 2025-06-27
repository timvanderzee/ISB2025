function [error, dX] = ThickEquilibrium(F, DRX, dDRXdt, k1, k2, kF, Ntot)

[J1, J2] = ThickFilament_Dynamics(F, DRX, k1, k2, kF, Ntot);

dX = J1 - J2;

error = dDRXdt - dX; 

end