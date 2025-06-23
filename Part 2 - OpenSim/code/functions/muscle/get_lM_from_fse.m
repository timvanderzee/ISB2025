function[lM, lMtilda, cos_alpha] = get_lM_from_fse(fse, lMT, params, kT)

lMo = params(2);
lTs = params(3);
alphao = params(4);

lTtilda = fse/kT + 1;

lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilda).^2); % [m]

lMtilda = lM./lMo; % []
cos_alpha = (lMT-lTs.*lTtilda)./lM;

end