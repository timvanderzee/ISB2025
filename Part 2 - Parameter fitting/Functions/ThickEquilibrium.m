function [error, dX] = ThickEquilibrium(Q0, F, DRX, dDRXdt, k1, k2, kF, act)

SRX = 1 - DRX;

J1 = k1 * (1 + kF * max(F,0)/max(act, 1e-5)) .* SRX;
J2 = k2 .* DRX;

dX = J1 - J2;

error = dDRXdt - dX; 

end