function Noverlap = calc_filament_overlap(parms, lce)
% CALC_FILAMENT_OVERLAP calculate filament overlap. by HXRyu.  
% second order approximation, maximum of 1 at lce0, 0 at +-0.5(half-sarc len) 

half_s_len_norm = parms.s/2/parms.h; 
curv_w = 0.5;

len = (lce - 0); % assuming 0 is optimal length
Noverlap = max(0.00001, -(len-half_s_len_norm*curv_w)*(len+half_s_len_norm*curv_w)/(half_s_len_norm*curv_w)^2);
% 0.0001 is there so that we avoid dividing by 0  

end

