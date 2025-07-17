function[lM, lMtilda, cos_alpha, dL] = get_lM_from_fse(fse, lMT, params, kT)

lMo = params(2);
lTs = params(3);
alphao = params(4);

lTtilda = fse/kT + 1;

dL = lMT-lTs.*lTtilda;

lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilda).^2); % [m]
cos_alpha = (lMT-lTs.*lTtilda)./lM;

% ignoring pennation
% lM = (lMT-lTs.*lTtilda);
% cos_alpha = 1;

lMtilda = lM./lMo; % []

end