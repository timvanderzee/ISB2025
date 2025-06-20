function [dlMdt] = CalculateDLMDT(vMtilda, params_OS)
%Calculate derivative of muscle length
vMtilda_ext  = vMtilda(1,:); 
vMtilda_flex = vMtilda(2,:);

lMo        = params_OS.MT(2,:); 
vMtildamax = params_OS.MT(5,:);
dlMdt_ext  = vMtilda_ext.*  vMtildamax(1)./ lMo(1);  % Extensor - change vMtilda
dlMdt_flex = vMtilda_flex.* vMtildamax(2)./ lMo(2);  % Flexor - change vMtilda
dlMdt      = [dlMdt_ext; dlMdt_flex]; 
end

