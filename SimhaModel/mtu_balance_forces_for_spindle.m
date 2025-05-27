function new_length = mtu_balance_forces_for_spindle(obj)

x_adj = fzero(@(x) tempForce(x,obj), 0, optimset('display','off','TolFun',1e-10));

new_length = obj.hs_length + x_adj;
end


function zeroF = tempForce(x,obj)
% Adjust for filament compliance
delta_x = x * obj.compliance_factor;
% Shift populations by interpolation
interp_positions = obj.x_bins - delta_x;
% THIS IS WHERE YOU SOLVE FOR BIN POPULATIONS %
temp_bin_pops = interp1(obj.x_bins,obj.bin_pops,interp_positions, ...
    'linear',0)';
cbF = obj.cb_number_density * obj.k_cb * 1e-9 * ...
    sum((obj.x_bins + obj.power_stroke).* temp_bin_pops');

pF = obj.k_passive * (obj.hs_length - obj.hsl_slack);

% Calculate what would be the new force based on the muscle
hsF = cbF + pF;

% Calculate what would be the new force based on the tendon
newTendonLength = obj.cmd_length - (obj.hs_length + delta_x);
if newTendonLength<0, newTendonLength=0;end
newTendonForce = obj.tendon_stiffness * newTendonLength;


% Compare the muscle and tendon estimate of force (goal is to minimize)
zeroF = hsF - newTendonForce;

end