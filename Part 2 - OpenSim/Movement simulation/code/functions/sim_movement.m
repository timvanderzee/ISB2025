function[tFW,sFW] = sim_movement(modelname, t0, tf, s0, input, osimModel, osimState, auxdata)

if strcmp(modelname, 'Hill')
    M = 1;
else
    M = 5;
end

x00 = zeros(1, M);

x_isom = zeros(auxdata.NMuscles,M);

disp(modelname)
disp('Simulating isometric ...')
t_input = [0; 5];
vec_kT = 35;

for m = 1:auxdata.NMuscles
    disp(m)
    x00(1) = auxdata.fse0(m);
    
    A = [input.act(m); input.act(m)];
    LMT = [input.lMT(m); input.lMT(m)];
    VMT = [input.vMT(m); input.vMT(m)];

    if ~strcmp(modelname, 'Hill')

        % Length and overlap
        fse = max(auxdata.fse0(m), 0);
        [~, lMtilda, cos_alpha] = get_lM_from_fse(fse, input.lMT(m), auxdata.params(:,m), vec_kT);
        % Parallel force
        [Fpe, ~] = get_parallel_force(lMtilda, auxdata.Fpparam);

        % Contractile element force
        FMce = max(fse./cos_alpha - Fpe, 0);

        delta = 1.9207;

        Q0 = FMce / delta;
        x00(2) = Q0;
       
    end
        
    % simulate
    [t,state] = ode15s(@TendonForceOdeVecSRS,t_input,x00,[],t_input, A,LMT,VMT, auxdata.params(:,m), auxdata.Fvparam, auxdata.Fpparam, auxdata.Faparam, vec_kT, modelname, auxdata.parms);

    x_isom(m,:) = state(end,:); 
  
end

%% simulate pre-movement
t_pre = linspace(t_input(1), t_input(2), 1000);
f = 5;
AMP = auxdata.AMP;
% T = 1/f;

x0 = x_isom;

if ~isfield(input, 'platacc') % check whether we're doing the pendulum test
    disp('Simulating pre-movement ...')
    for m = 4:5
        
%         disp(num2str(m))
        
        % activation
        A = input.act(m) * ones(size(t_pre));
        
        % pre movement
        LMT_pre = (1 + AMP(m)) * input.lMT(m) - AMP(m) * input.lMT(m) * cos(2*pi*f*t_pre);
        VMT_pre = 2 * pi * f * AMP(m) * input.lMT(m) * sin(2*pi*f*t_pre);
        
        % simulate
        [t,state] = ode15s(@TendonForceOdeVecSRS,t_input,x_isom(m,:),[],t_pre, A,LMT_pre,VMT_pre, auxdata.params(:,m), auxdata.Fvparam, auxdata.Fpparam, auxdata.Faparam, vec_kT, modelname, auxdata.parms);

        x0(m,:) = state(end,:); 

    end
end

%% simulate movement

disp('Simulating movement ...')
auxdata.kT = vec_kT;

% simulate
[tFW,sFW] = ode15s(@compute_state_derivatives,[t0 tf],[s0(:); x0(:)],[], input, osimModel, osimState, auxdata, modelname);

end
