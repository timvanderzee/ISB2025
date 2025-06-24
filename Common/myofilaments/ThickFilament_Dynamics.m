function[J1, J2] = ThickFilament_Dynamics(F, DRX, k1, k2, kF, Ntot)

SRX = 1 - DRX;

% J1 = k1 * (1 + kF * max(F,0)/max(act, 1e-5)) .* SRX;
J1 = k1 * (1 + kF * F/Ntot) .* SRX;
J2 = k2 .* DRX;
    
end