function[tFW,sFW] = sim_movement(modelname, t0, tf, s0, input, osimModel, osimState, auxdata)

if strcmp(modelname, 'Hill')
    M = 1;
else
    M = 5;
end

x00 = zeros(1, M);
x_isom = zeros(auxdata.NMuscles,M);
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
    [t,state] = ode15s(@TendonForceOdeVecSRS,[0,10],x00,[],t_input, A,LMT,VMT, auxdata.params(:,m), auxdata.Fvparam, auxdata.Fpparam, auxdata.Faparam, vec_kT, modelname);

    fse = state(:,1);
    fse_isom(m) = fse(end);
    x_isom(m,:) = state(end,:); 
  
end

%% simulate pre-movement
% if ~isfield(input, 'platacc') % check whether we're doing the pendulum test
%     for m = 1:auxdata.NMuscles
%         A = [input.act(m); input.act(m)];
%         LMT = [input.lMT(m); input.lMT(m)];
%         VMT = [input.vMT(m); input.vMT(m)];
% 
%         % simulate
%         [t,state] = ode15s(@TendonForceOdeVecSRS,[0,10],x00,[],t_input, A,LMT,VMT, auxdata.params(:,m), auxdata.Fvparam, auxdata.Fpparam, auxdata.Faparam, vec_kT, modelname);
% 
%         fse = state(:,1);
%         fse_isom(m) = fse(end);
%         x_isom(m,:) = state(end,:); 
% 
%     end
% end
%% simulate movement
auxdata.kT = vec_kT;

% simulate
[tFW,sFW] = ode15s(@compute_state_derivatives,[t0 tf],[s0(:); x_isom(:)],[], input, osimModel, osimState, auxdata, modelname);


end
