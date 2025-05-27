function [hsE, mtData] = musTenDriverForVdz4State(t,delta_cdl,sarcE)

hsE = halfSarc_4stateFromVdz();
hsE.power_stroke = sarcE.power_stroke;
hsE.hsl_slack = sarcE.hsl_slack;
hsE.k_passive = sarcE.k_passive;
hsE.compliance_factor = sarcE.compliance_factor;
% hsE.tendon_stiffness=sarcE.tendon_stiffness;
hsE.passive_force_mode = sarcE.passive_force_mode;
hsE.passive_L = sarcE.passive_L;
hsE.passive_sigma = sarcE.passive_sigma;
hsE.bin_min = sarcE.bin_min;      % min x value for myosin distributions in nm
hsE.bin_max = sarcE.bin_max;       % max x value for myosin distributions in nm
hsE.bin_width = sarcE.bin_width;
hsE.k_cb = sarcE.k_cb ;
hsE.thick_filament_length = sarcE.thick_filament_length;
hsE.thin_filament_length = sarcE.thin_filament_length;
hsE.bare_zone_length = sarcE.bare_zone_length;
hsE.k_falloff = sarcE.k_falloff; 
hsE.max_rate = sarcE.max_rate;
hsE.cb_number_density = sarcE.cb_number_density;

for i = 1:numel(t)
    hsE.pCa_perStep = sarcE.pCa(i);
    
    
    if sarcE.isTendon==1
        hsE.tendon_stiffness = sarcE.tendon_stiffness_vec(i);
        if i > 1
            time_step = t(i) - t(i-1);
        else
            hsE.cmd_length = sarcE.hs_length;
            hsE.hs_length = sarcE.hs_length;
            time_step = t(2) - t(1);
            
            x_new = mtu_balance_forces_for_spindle(hsE);
            x_adj = x_new - hsE.hs_length;
            hsE.forwardStep(0,x_adj,0,0,1);
            
        end
        
        
        hsE.forwardStep(time_step,0,0,1,0)
        
        x_new = mtu_balance_forces_for_spindle(hsE);
        x_adj = x_new - hsE.hs_length;
        hsE.forwardStep(0,x_adj,delta_cdl(i),0,1);
    else
        if i > 1
            time_step = t(i) - t(i-1);
        else
            hsE.cmd_length = sarcE.hs_length;
            hsE.hs_length = sarcE.hs_length;
            time_step = t(2) - t(1);
        end
        delta_hsl = delta_cdl(i);
        hsE.forwardStep(time_step,delta_hsl,delta_cdl(i),1,1);
    end
    
    % if mod(i,100)==0, disp(['done with t' num2str(i)]);end
    
    mtData.t(i) = t(i);
    
    mtData.f_bound(i) = sum(hsE.bin_pops);
    mtData.f_overlap(i) = hsE.f_overlap;
    mtData.cb_force(i) = hsE.cb_force;
    
    
    mtData.passive_force(i) = hsE.passive_force;
    mtData.hs_force(i) = hsE.hs_force;
    mtData.hs_length(i) = hsE.hs_length;
    mtData.cmd_length(i) = hsE.cmd_length;
    mtData.bin_pops(:,i) = hsE.bin_pops;
    mtData.no_detached(i) = hsE.no_detached;
end

end

