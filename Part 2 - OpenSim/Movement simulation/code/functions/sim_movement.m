function[tFW,sFW] = sim_movement(modelname, t0, tf, s0, input, osimModel, osimState, auxdata)

if strcmp(modelname, 'Hill')
    M = 1;
else
    M = 5;
end

x00 = zeros(1, M);
x_isom = zeros(auxdata.NMuscles,M);

disp('Simulation isometric ...')
t_input = [0; 5];

% ls = {'-','--'};
% act = 1;
vec_kT = 35;

for m = 1:auxdata.NMuscles
    A = [input.act(m); input.act(m)];
    LMT = [input.lMT(m); input.lMT(m)];
    VMT = [input.vMT(m); input.vMT(m)];

    % simulate
    [t,state] = ode15s(@TendonForceOdeVecSRS,t_input,x00,[],t_input, A,LMT,VMT, auxdata.params(:,m), auxdata.Fvparam, auxdata.Fpparam, auxdata.Faparam, vec_kT, modelname);

    x_isom(m,:) = state(end,:); 
  
end

%% simulate pre-movement
t_pre = linspace(t_input(1), t_input(2), 1000);
f = 5;
AMP = auxdata.AMP;
% T = 1/f;

x0 = x_isom;

if ~isfield(input, 'platacc') % check whether we're doing the pendulum test
    disp('Simulation pre-movement ...')
    for m = 5
        
%         disp(num2str(m))
        
        % activation
        A = input.act(m) * ones(size(t_pre));
        
        % pre movement
        LMT_pre = (1 + AMP) * input.lMT(m) - AMP * input.lMT(m) * cos(2*pi*f*t_pre);
        VMT_pre = 2 * pi * f * AMP * input.lMT(m) * sin(2*pi*f*t_pre);
        
        % simulate
        [t,state] = ode15s(@TendonForceOdeVecSRS,t_input,x_isom(m,:),[],t_pre, A,LMT_pre,VMT_pre, auxdata.params(:,m), auxdata.Fvparam, auxdata.Fpparam, auxdata.Faparam, vec_kT, modelname);

        x0(m,:) = state(end,:); 

    end
end

%% simulate movement

disp('Simulation movement ...')
auxdata.kT = vec_kT;

% simulate
[tFW,sFW] = ode15s(@compute_state_derivatives,[t0 tf],[s0(:); x0(:)],[], input, osimModel, osimState, auxdata, modelname);


end
