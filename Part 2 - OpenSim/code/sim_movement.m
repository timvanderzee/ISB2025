function[tFW,sFW] = sim_movement(modelname, t0, tf, s0, input, osimModel, osimState, auxdata)

if strcmp(modelname, 'Hill')
    M = 1;
else
    M = 5;
end

x00 = zeros(1, M);
x_isom = zeros(auxdata.NMuscles,M);
lMtilda_isom = zeros(auxdata.NMuscles,1);
fse_isom = zeros(auxdata.NMuscles,1);

t_input = [0; 10];

% ls = {'-','--'};
% act = 1;
vec_kT = 35;

for m = 1:auxdata.NMuscles
    A = [input.act(m); input.act(m)];
    LMT = [input.lMT(m); input.lMT(m)];
    VMT = [input.vMT(m); input.vMT(m)];

    % simulate
    [t,state] = ode15s(@TendonForceOdeVecSRS,[0,10],x00,[],t_input, A,LMT,VMT, auxdata.params(:,m), auxdata.Fvparam, auxdata.Fpparam, auxdata.Faparam, input.lMT(m), input.vMT(m), vec_kT, modelname);

    fse = state(:,1);
    fse_isom(m) = fse(end);
    x_isom(m,:) = state(end,:);

%     figure(1)
%     nexttile
%     plot(t, fse); hold on

%     FMo = ones(size(fse,1),1)*auxdata.params(1,m);
    lMo = ones(size(fse,1),1)*auxdata.params(2,m);
    lTs = ones(size(fse,1),1)*auxdata.params(3,m);
    alphao = ones(size(fse,1),1)*auxdata.params(4,m);

%     FT = fse .* FMo;
    lTtilda = fse/vec_kT + 1;
    lM = sqrt((lMo.*sin(alphao)).^2+(input.lMT(m)-lTs.*lTtilda).^2);
    lMtilda = lM./lMo;
    lMtilda_isom(m) = lMtilda(end);
  
  
end


%% simulate movement
auxdata.ksrs = 0;
auxdata.kT = vec_kT;

% x0 = reshape([x_isom(:,j)], size(x_isom{1,j},2), M)';

input.lMtilda_isom = lMtilda_isom;
[tFW,sFW] = ode15s(@compute_state_derivatives,[t0 tf],[s0(:); x_isom(:)],[], input, osimModel, osimState, auxdata, modelname);


end
