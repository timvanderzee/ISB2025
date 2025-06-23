%% S1
addpath(genpath('C:\Users\u0167448\Documents\GitHub\BIOMUS\soft'))

mainfolder = 'C:\From_Friedl_backup\GT_essential';
mainfolder = 'C:\Users\timvd\Documents';
cd([mainfolder, '\dynamic_optimization_full'])

% Add paths for software and dataset
% ----------------------------------
warning off
voeg_paden_toe('SCP_simple');
warning on

load_simulation_parameters

trials = [15 16 28 35 36 37 42];
trials = 15;

% without SRS 50ms, without SRS 100ms, with SRS 50ms, with SRS 100ms

%                  1          2          3          4          5          6          7        8           9        10            11        12
%                red     fade red     green    fade green    blue     fade blue   black     gray       yellow    orange        plum      lime
% color_palette= [255 0 0; 255 150 150; 0 114 54; 153 255 102; 0 0 255; 44 176 207; 0 0 0; 155 155 155; 255 153 0; 255 185 0; 153 51 102; 171 218 77]*(1/255);
% color=color_palette([12 3 11 5 10 1 11 5 10 4 3 5 12 7 8 9 1 10 4 6 11 3 2 5 12],:);

color_palette= [0 114 54; 75 171 80; 153 228 108; 153 51 102; 205 105 210; 221 160 230]*(1/255);
color=color_palette([1 2 3 4 5 6],:);

dofNames = {'ankle' 'knee' 'hip' 'hip AB' 'hip ROT'};

avgQ = [-10 -10 -10; -10 -10 -10];

% filenames = {'_ebase5_0_50ms_ksrs0_web3', '_ebase5_0_50ms_ksrs0_kT70_web30', '_ebase5_0_50ms_ksrs0_web0', '_ebase5_0_50ms_ksrs280_web3', '_ebase5_0_50ms_ksrs280_kT70_web30','_ebase5_0_50ms_ksrs280_web0'};
% filenames = {'_ebase3_0_50ms_ksrs0_web30_rev', '_ebase3_0_50ms_ksrs0kT70_web30_rev', '_ebase3_0_50ms_ksrs0_web0_rev', '_ebase3_0_50ms_ksrs280_web30_rev', '_ebase3_0_50ms_ksrs280kT70_web30_rev','_ebase3_0_50ms_ksrs280_web0_rev'};
filenames = {'_ebase3_min10_50ms_ksrs0_web30_rev', '_ebase3_min10_50ms_ksrs0kT70_web30_rev', '_ebase3_min10_50ms_ksrs0_web0_rev', '_ebase3_min10_50ms_ksrs280_web30_rev', '_ebase3_min10_50ms_ksrs280kT70_web30_rev','_ebase3_min10_50ms_ksrs280_web0_rev'};
vec_ksrs = [0 0 0 280 280 280];
vec_kT = [35 70 35 35 70 35];
n_est = length(filenames);

RMSEall = zeros(length(trials), n_est, 3);
RMSE50all = zeros(length(trials), n_est, 3);

avg_eb = zeros(length(trials), n_est);

figure()
for tr = 1:length(trials)
    RMSE = zeros(n_est, 3);
    RMSE50 = zeros(n_est, 3);
    
    trial = trials(tr)
    
    IKpath = [mainfolder,'\Models\SCP_gait10dof18musc_Trial' num2str(trial) '_ks.mot'];
    IKsolution = importdata(IKpath);
    
    % PARAMETERS
    auxdata.NMuscles = M;
    
    load Fvparam
    Fvparam(1) = 1.475*Fvparam(1);
    Fvparam(2) = 0.25*Fvparam(2);
    Fvparam(3) = Fvparam(3) + 0.75;
    Fvparam(4) = Fvparam(4) - 0.027;
    auxdata.Fvparam = Fvparam;
    
    load Faparam
    auxdata.Faparam = Faparam;
    
    e0 = 0.6;
    kpe = 4;
    t50 = exp(kpe * (0.2 - 0.10e1) / e0);
    pp1 = (t50 - 0.10e1);
    t7 = exp(kpe);
    pp2 = (t7 - 0.10e1);
    Fpparam = [pp1;pp2];
    auxdata.Fpparam = Fpparam;
    
    
    load([mainfolder, '\OSIM_posturePerts_Friedl\matlab_files\PlatOnsets.dat'])
    load_simulation_parameters
    
    eval(['load Trial' num2str(trial) '_InputDataTime'])
    tOnset = PlatOnsets(trial);
    time = InputDataTime;
    
    indeces = find(time > tOnset);
    sampleBegin = indeces(1) - 1 + round(fs * 0.02); % start simulation 100 ms before perturbation
    t0 = InputDataTime(sampleBegin);
    sampleEnd = indeces(1) - 1 + round(fs * 0.12); % end simulation 50 ms after platform onset
    tf = InputDataTime(sampleEnd);
    
    
    % INPUT
    
    pelvis_tilt_exp = IKsolution.data(sampleBegin,2);
    pelvis_tilt_mod = IKsolution.data(sampleBegin,7) + IKsolution.data(sampleBegin,6) +IKsolution.data(sampleBegin,5);
    delta_pelvis_tilt = pelvis_tilt_exp + pelvis_tilt_mod;
    
    s0 = zeros(8,1);
    s0(2:4) = [-IKsolution.data(sampleBegin,7)+delta_pelvis_tilt IKsolution.data(sampleBegin,6) -IKsolution.data(sampleBegin,5)]/180*pi;
    
    
    
    load([mainfolder,'\OSIM_posturePerts_Friedl\processedData\Trial' num2str(trial)], 'data');
    
    platpos = data.analog.platemotion(:,2)/100;
    input.platacc = data.analog.platemotion(:,6);
    input.time = [1:1:length(data.analog.platemotion(:,6))]/data.analog.samplerate;
    
    % Import the OpenSim modeling classes
    import org.opensim.modeling.*
    
    % Read in the osim model
    osimModel = Model([mainfolder,'\Models\SCP_gait4dof9musc.osim']); %
    osimModel3D = Model([mainfolder,'\Models\SCP_gait10dof18musc.osim']); %
    nom = length(namen);
    params = zeros(5, nom);
    
    muscles = osimModel3D.getMuscles();
    
    for i = 1:nom
        muscle = muscles.get(namen{i});
        params(3,i) = muscle.getTendonSlackLength();
        params(2,i) = muscle.getOptimalFiberLength();
        params(1,i) = 2*muscle.getMaxIsometricForce(); % one leg model
        params(4,i) = muscle.getPennationAngleAtOptimalFiberLength();
    end
    params(5,:) = 10*params(2,:);
    auxdata.params = params;
    
    % Initialize the model (this builds the system and initialize the state)
    osimState = osimModel.initSystem();
    
    for k = 1:length(filenames)
        
        load(['musdyn_SRS_dlmt_Trial' num2str(trial) filenames{k}], 'solution')
        input.act =  max(solution.parameter(1:9), 0.01);
        avg_eb(tr,k) = mean(input.act);
        
        % Update model state with current values
        osimState.setTime(0);
        numVar = osimState.getNY();
        for i = 0:numVar-1
            osimState.updY().set(i, s0(i+1,1));
        end
        
        % Compute the isometric muscle length!
        
        % Get the number of states, coordinates, muscles and controls from the model;
        % in this case the number of controls equals the number of muscles
        Nstates       = osimModel.getNumStateVariables();
        Ncontrols     = osimModel.getNumControls();
        Ncoord        = osimModel.getNumCoordinates();
        % model_muscles = osimModel.getMuscles();
        % Nmuscles      = model_muscles.getSize();
        
        % Get the names of the states from the model
        stateNames = cell(Nstates,1);
        for i = 1:Nstates
            stateNames(i,1) = cell(osimModel.getStateVariableNames().getitem(i-1));
        end
        
        auxdata.NStates = Nstates;
        
        % % Get the names of the states from the model
        % states_all = cell(Nstates,1);
        % for i = 1:Nstates
        %    states_all(i,1) = cell(osimModel.getStateVariableNames().getitem(i-1));
        % end
        
        % Get the names of the controls/muscles from the model
        osimModel.computeStateVariableDerivatives(osimState); % To be able to access muscle fiber velocity.
        Actuators = osimModel.getActuators();
        % Muscles = osimModel.getMuscles();
        controlNames = cell(Ncontrols,1);
        controls_all = cell(Ncontrols,1);
        lMT = zeros(1,M);
        vMT = zeros(1,M);
        for i = 1:Ncontrols
            currentActuator = Actuators.get(i-1);
            if strcmp(currentActuator.getConcreteClassName(), 'PathActuator')
                shouldBePathAct = osimModel.getActuators().get(currentActuator.getName());
                currentPathAct = PathActuator.safeDownCast(shouldBePathAct);
                lMT(i) = currentPathAct.getLength(osimState);
                vMT(i) = currentPathAct.getLengtheningSpeed(osimState);
            end
            controls_all(i,1) = cell(currentActuator.getName());
            controlNames(i,1) = cell(currentActuator.getName());
        end
        
        t_input = [0;1];
        lMtilda_isom = zeros(M,1);
        fse_isom = zeros(M,1);
        x0 = [0 0 0];
        
        %%
        FMo = params(1,:);
        lMo = params(2,:);
        lTs = params(3,:);
        alphao = params(4,:);
        
        l0 = (sqrt((lMT - lTs).^2 + (lMo .* sin(alphao)).^2)) ./ lMo;
        
        return
        
        %% simulate isometric 
        close all
        clc
        
        t_input = [0; 10];
        
        modelnames = {'TendonForceOdeVecSRS', 'TendonForceOdeVecSRS_BP'};
        x0 = {0, [0 0 0 0 0 0]};
        x_isom = cell(M,2);
        
        for m = 1:M
            figure(1)
            nexttile
            
            for j = 1:2
                A = [input.act(m); input.act(m)];
                LMT = [lMT(m); lMT(m)];
                VMT = [vMT(m); vMT(m)];
                
                % simulate
                [t,state] = ode15s(eval(['@', modelnames{j}]),[0,10],x0{j},[],t_input, A,LMT,VMT,params(:,m), Fvparam, Fpparam, Faparam, lMT(m), 0, vec_kT(k));

                fse = state(:,1);
                fse_isom(m,j) = fse(end);
                x_isom{m,j} = state(end,:);
                
                plot(t, fse); hold on              
                
                FMo = ones(size(fse,1),1)*params(1,m);
                lMo = ones(size(fse,1),1)*params(2,m);
                lTs = ones(size(fse,1),1)*params(3,m);
                alphao = ones(size(fse,1),1)*params(4,m);

                FT = fse .* FMo;
                lTtilda = fse/vec_kT(k) + 1;
                lM = sqrt((lMo.*sin(alphao)).^2+(lMT(m)-lTs.*lTtilda).^2);
                lMtilda = lM./lMo;
                lMtilda_isom(m,j) = lMtilda(end);
            end
        end
        
        return
        
        %% simulate movement
        close all
        
        auxdata.ksrs = vec_ksrs(k);
        auxdata.kT = vec_kT(k);
        modelnames = {'compute_state_derivatives', 'compute_state_derivatives_BP'};
        
        for j = 1:2

            x0 = reshape([x_isom{:,j}], size(x_isom{1,j},2), M)';

%             if j == 2
%                 x0(:,end-1) = x0(:,end-1) * .9;
%             end
            
            input.lMtilda_isom = lMtilda_isom(:,j);
            [tFW,sFW] = ode15s(eval(['@',modelnames{j}]),[t0,tf],[s0; x0(:)],[], input, osimModel, osimState, auxdata);

            fse = sFW(:,auxdata.NStates+1:auxdata.NStates+auxdata.NMuscles);

            figure(1)
            for m = 1:M
                subplot(3,3,m)
                plot(tFW, fse(:,m)); hold on
                box off
            end
        
        end
        
        return
        
        %%
        close all
        clc
        
        for m = 1:M
            figure(1)
            nexttile
            
            LMT = [lMT(m); lMT(m)];
            VMT = [vMT(m); vMT(m)];
            
            VMT = .5 * ones(1,2);
            
            A = [input.act(m); input.act(m)];
            A = [1 1];
            
            for j = 1:2
                
                if j == 1
                    [t,state] = ode15s(@TendonForceOdeVecSRS,[0,.1], fse_isom(m),[],t_input, A, LMT,VMT,params, m, Fvparam, Fpparam, Faparam, lMT(m), 0, vec_kT(k));
                elseif j == 2
                    [t,state] = ode15s(@TendonForceOdeVecSRS_BP,[0,.1], x_isom(m,1:3),[],t_input, A,LMT,VMT,params, m, Fvparam, Fpparam, Faparam, lMT(m), 0, vec_kT(k));
                end
                
                fse = state(:,1);
                
                plot(t, fse); hold on
            end
        end
        
        
        %%
        delta = 2.5;
        
        Q0 = state(:,2);
        Q1 = state(:,1)/delta - Q0;
                Q2 = state(:,3);
        
        % interpolate
        ti = linspace(0, max(t),1000);
        xi = interp1(t, [Q0 Q1 Q2], ti);
        si = linspace(-2,2,1000);
        
        Q0c = max(Q0, eps); % correct because can't be zero
        
        eps = 1e-6;
        p = xi(:,2)./(xi(:,1)+eps); % Eq. 52
                q = sqrt(max(xi(:,3)./(xi(:,1)+eps) - p.^2, 1e-3));  % Eq. 52
        %         q = .1;
        
        % n
        n_func = @(Q0, p, q, xi) Q0 ./ (sqrt(2*pi)*q) * exp(-((xi-p).^2) / (2*q^2));  % Eq. 52, modified
        
        figure(2)
        
        for i = 1:length(xi)
            n = n_func(xi(i,1),p(i),q(i), si);
            
            plot(si, n);
            axis([-2 2 0 1])
            drawnow
%                         pause
        end
        
 
        %%
        
%         dx = compute_state_derivatives_BP(t0, [s0; x_isom(:,1); x_isom(:,2)], input, osimModel, osimState, auxdata);
        dx =  compute_state_derivatives(t0,[s0;fse_isom(:,1)], input, osimModel, osimState, auxdata);
        
        
        %%
        
  
        
        

        
        %%
        clc
        close all
        
        for m = 1:M
            figure(1)
            nexttile
            
            for j = 1:2
                if j == 1
                    
                    states = sFWB(:,1:auxdata.NStates);
                    fse = sFWB(:,auxdata.NStates+1:auxdata.NStates+auxdata.NMuscles);
                    t = tFWB;
                else
                    
                    states = sFWH(:,1:auxdata.NStates);
                    fse = sFWH(:,auxdata.NStates+1:auxdata.NStates+auxdata.NMuscles);
                    t = tFWH;
                    
                end
                %         Q0 = sFW(:,auxdata.NStates+auxdata.NMuscles+1:(auxdata.NStates+auxdata.NMuscles*2));
                %         Q2 = sFW(:,(auxdata.NStates+auxdata.NMuscles*2)+1:(auxdata.NStates+auxdata.NMuscles*3));
                
                
                
                plot(t, fse(:,m)); hold on
                box off
                xlabel('Time (s)')
                ylabel('fse (-)')
            end
            
        end
        
        
        legend('Biophysical','Hill','location','best')
        legend boxoff


        
        %%
 
         
        % Write model states to an OpenSim STO file
        StatesData.name = ['SCP_Trial' num2str(trial) '_' num2str(k) '_FW_States'];
        StatesData.nRows = length(tFW);
        StatesData.nColumns = Nstates + 1;
        StatesData.labels = [{'time'}; stateNames]';
        StatesData.inDegrees = false;
        StatesData.data = [tFW,sFW(:,1:Nstates)];
        writeOpenSimStatesFile(StatesData)
        
        % Create data structure for the controls file
        ControlData.name = ['SCP_Trial' num2str(trial) '_' num2str(k) '_FW_Controls'];
        ControlData.nRows = size(tFW);
        ControlData.nColumns = Ncontrols+1; %All the controls + time
        ControlData.inDegrees = false;
        ControlData.labels= [{'time'}; controlNames]';
        
        mass = 2000;
        acc = interp1(input.time, input.platacc, tFW)*9.81;
        controlFX = mass*acc/10000;
        ControlData.data = [tFW,sFW(:,1+Nstates:end), controlFX];
        writeOpenSimControlFile(ControlData)
        
        %
        
        positionsFW = sFW(:, 1:4);
        positionsFW(:, 2:4) = positionsFW(:, 2:4)*180/pi;
        qFW = [tFW positionsFW];
        positionsFW(:, 1) = 1000*positionsFW(:, 1);
        
        t_IK = IKsolution.data(sampleBegin:sampleEnd,1);
        s_IK = [-IKsolution.data(sampleBegin:sampleEnd,7)+delta_pelvis_tilt IKsolution.data(sampleBegin:sampleEnd,6) -IKsolution.data(sampleBegin:sampleEnd,5)];
        
        posFW = interp1(tFW, positionsFW(:, 2:4), t_IK(t_IK<tFW(end)));
        
        t50 = (t0+tf)/2;
        posFW50 = interp1(tFW, positionsFW(:, 2:4), t_IK(t_IK<t50));
        
        RMSE(k,:) = sqrt(mean((posFW - s_IK(t_IK<tFW(end),:)).^2));
        RMSE50(k,:) = sqrt(mean((posFW50 - s_IK(t_IK<t50,:)).^2));
        
        
        if trial == 15 | trial == 35
            if trial == 15
                extr = 1;
            else
                extr = 2;
            end
            
            for i = 1:3
                subplot(4,7,(extr-1)*7+i)
                plot(1000*(tFW-t0), positionsFW(:,i+1), 'Color', color(k,:), 'LineWidth', 2); hold on;
            end
            
        end
        
    end
    
    %%
    RMSEall(tr,:,:) = RMSE;
    RMSE50all(tr,:,:) = RMSE50;
    
    if trial == 15 | trial == 35
        for i = 1:3
            
            if trial == 15
                extr = 1;
            else
                extr = 2;
            end
            
            subplot(4,7,(extr-1)*7+i)
            
            if extr == 2
                xlabel('t [ms]')
            end
            if i == 1
                plot(1000*(IKsolution.data(sampleBegin:sampleEnd,1)-t0), -IKsolution.data(sampleBegin:sampleEnd,7)+delta_pelvis_tilt, 'k', 'LineWidth', 3); hold on;
                ylabel('q [^o]')
                axis([0 100 avgQ(extr,i)-10 avgQ(extr,i)+10])
                set(gca,'XTick',[0 50 100])
                set(gca,'YTick',[-20 -10 0])
                box off
            elseif i == 2
                plot(1000*(IKsolution.data(sampleBegin:sampleEnd,1)-t0), IKsolution.data(sampleBegin:sampleEnd,6), 'k', 'LineWidth', 3); hold on;
                %         ylabel('q [^o]')
                axis([0 100 avgQ(extr,i)-10 avgQ(extr,i)+10])
                set(gca,'XTick',[0 50 100])
                set(gca,'YTick',[-20 -10 0])
                box off
            else
                plot(1000*(IKsolution.data(sampleBegin:sampleEnd,1)-t0), -IKsolution.data(sampleBegin:sampleEnd,5), 'k', 'LineWidth', 3); hold on;
                %         ylabel('q [^o]')
                axis([0 100 avgQ(extr,i)-10 avgQ(extr,i)+10])
                set(gca,'XTick',[0 50 100])
                set(gca,'YTick',[-20 -10 0])
                box off
                if extr == 1
                    legend('no SRS - w_{opt}', 'no SRS - w_{opt} - k_T high', 'no SRS - w_{0}', 'SRS - w_{opt}', 'SRS - w_{opt} - k_T high', 'SRS  - w_{0}', 'measured')
                end
            end
            if extr == 1
                title(dofNames{i})
            end
            
        end
    end
    
    
end

avg_avg_eb = mean(avg_eb)

avgRMSE50 = squeeze(mean(RMSE50all,1));
stdRMSE50 = squeeze(std(RMSE50all,1,1));

avgRMSE = squeeze(mean(RMSEall,1));
stdRMSE = squeeze(std(RMSEall,1,1));

subplot(4,7,[15:17])
x = 1:3;
bar_h = bar(x, avgRMSE50'); hold on;
set(gca,'XTickLabel',{'ankle' 'knee' 'hip'});
% set(gca,'YTick',[0 1 2 3])
colormap(color(1:6,:))
sb = 0.1325;
X = [x'-2.5*sb x'-1.5*sb x'-0.5*sb x'+0.5*sb x'+1.5*sb x'+2.5*sb];
h=errorbar(X,avgRMSE50',min(avgRMSE50', stdRMSE50') ,stdRMSE50','k'); set(h,'linestyle','none')
ylabel('RMSE [^o]')
axis([0.5 3.5 0 6])
set(gca,'YTick',[0 5])
box off

subplot(4,7,[19:21])
x = 1:3;
bar_h = bar(x,avgRMSE'); hold on;
set(gca,'XTickLabel',{'ankle' 'knee' 'hip'});
% set(gca,'YTick',[0 1 2 3])
colormap(color(1:6,:))
sb = 0.1325;
X = [x'-2.5*sb x'-1.5*sb x'-0.5*sb x'+0.5*sb x'+1.5*sb x'+2.5*sb];
h=errorbar(X,avgRMSE',min(avgRMSE', stdRMSE') ,stdRMSE','k'); set(h,'linestyle','none')
ylabel('RMSE [^o]')
axis([0.5 3.5 0 6])
set(gca,'YTick',[0 5])
box off

legend('no SRS - w_{opt}', 'no SRS - w_{opt} - k_T high', 'no SRS - w_{0}', 'SRS - w_{opt}', 'SRS - w_{opt} - k_T high', 'SRS  - w_{0}')

disp('S1')
avg50 = mean(avgRMSE50,2)
std50 = std(avgRMSE50, 0, 2)


%% S2

% Add paths for software and dataset
% ----------------------------------
warning off
voeg_paden_toe('AMPC04_simple');
warning on

load_simulation_parameters

trials = [39 48 82 151 14 22 49 56 89 162 10 15];

RMSEall = zeros(length(trials), 4, 3);

dofNames = {'ankle' 'knee' 'hip' 'hip AB' 'hip ROT'};

avgQ = [-7 -9 -20; -7 -9 -20];

% filenames = {'_ebase5_0_50ms_ksrs0_web30', '_ebase5_0_50ms_ksrs0_kT70_web30', '_ebase5_0_50ms_ksrs0_web0', '_ebase5_0_50ms_ksrs280_web30', '_ebase5_0_50ms_ksrs280_kT70_web30','_ebase5_0_50ms_ksrs280_web0'};
% filenames = {'_ebase3_0_50ms_ksrs0_web30_rev', '_ebase3_0_50ms_ksrs0kT70_web30_rev', '_ebase3_0_50ms_ksrs0_web0_rev', '_ebase3_0_50ms_ksrs280_web30_rev', '_ebase3_0_50ms_ksrs280kT70_web30_rev','_ebase3_0_50ms_ksrs280_web0_rev'};
filenames = {'_ebase3_min10_50ms_ksrs0_web30_rev', '_ebase3_min10_50ms_ksrs0kT70_web30_rev', '_ebase3_min10_50ms_ksrs0_web0_rev', '_ebase3_min10_50ms_ksrs280_web30_rev', '_ebase3_min10_50ms_ksrs280kT70_web30_rev','_ebase3_min10_50ms_ksrs280_web0_rev'};
vec_ksrs = [0 0 0 280 280 280];
vec_kT = [35 70 35 35 70 35];
n_est = length(filenames);

RMSEall = zeros(length(trials), n_est, 3);
RMSE50all = zeros(length(trials), n_est, 3);
avg_eb = zeros(length(trials), n_est);

for tr = 1:length(trials)
    RMSE = zeros(n_est, 3);
    RMSE50 = zeros(n_est, 3);
    
    trial = trials(tr)
    
    IKpath = [mainfolder,'\AMPC04\ModelsForForwardSimulation\SCP_gait10dof18musc_Trial' num2str(trial) '.mot'];
    IKsolution = importdata(IKpath);
    
    % PARAMETERS
    auxdata.NMuscles = M;
    
    load Fvparam
    Fvparam(1) = 1.475*Fvparam(1);
    Fvparam(2) = 0.25*Fvparam(2);
    Fvparam(3) = Fvparam(3) + 0.75;
    Fvparam(4) = Fvparam(4) - 0.027;
    auxdata.Fvparam = Fvparam;
    
    load Faparam
    auxdata.Faparam = Faparam;
    
    e0 = 0.6;
    kpe = 4;
    t50 = exp(kpe * (0.2 - 0.10e1) / e0);
    pp1 = (t50 - 0.10e1);
    t7 = exp(kpe);
    pp2 = (t7 - 0.10e1);
    Fpparam = [pp1;pp2];
    auxdata.Fpparam = Fpparam;
    
    
    load([mainfolder,'\AMPC04\PlatOnsets.dat'])
    load_simulation_parameters
    
    eval(['load Trial' num2str(trial) '_InputDataTime'])
    tOnset = PlatOnsets(trial);
    time = InputDataTime;
    
    indeces = find(time > tOnset);
    sampleBegin = indeces(1) - 1 + round(fs * 0.02); % start simulation 100 ms before perturbation
    t0 = InputDataTime(sampleBegin);
    sampleEnd = indeces(1) - 1 + round(fs * 0.12); % end simulation 50 ms after platform onset
    tf = InputDataTime(sampleEnd);
    
    
    % INPUT
    
    t_IK_all = IKsolution.data(:,1);
    ind_IK = find(t_IK_all >= t0 & t_IK_all <= tf);
    
    pelvis_tilt_exp = IKsolution.data(ind_IK(1),2);
    pelvis_tilt_mod = IKsolution.data(ind_IK(1),7) + IKsolution.data(ind_IK(1),6) +IKsolution.data(ind_IK(1),5);
    delta_pelvis_tilt = pelvis_tilt_exp + pelvis_tilt_mod;
    
    s0 = zeros(8,1);
    s0(2:4) = [-IKsolution.data(ind_IK(1),7)+delta_pelvis_tilt IKsolution.data(ind_IK(1),6) -IKsolution.data(ind_IK(1),5)]/180*pi;
    
    
    
    load([mainfolder,'\AMPC04\ExperimentalData\trial' num2str(trial)], 'data');
    
    platpos = data.analog.platemotion(:,2)/100;
    input.platacc = data.analog.platemotion(:,6);
    input.time = [1:1:length(data.analog.platemotion(:,6))]/data.analog.samplerate;
    
    % Import the OpenSim modeling classes
    import org.opensim.modeling.*
    
    % Read in the osim model
    osimModel = Model([mainfolder,'\AMPC04\ModelsForForwardSimulation\AMPC04_gait4dof9musc.osim']); %
    osimModel3D = Model([mainfolder, '\AMPC04\ModelsForForwardSimulation\AMPC04_gait10dof18musc.osim']); %
    nom = length(namen);
    params = zeros(5, nom);
    
    muscles = osimModel3D.getMuscles();
    
    for i = 1:nom
        muscle = muscles.get(namen{i});
        params(3,i) = muscle.getTendonSlackLength();
        params(2,i) = muscle.getOptimalFiberLength();
        params(1,i) = 2*muscle.getMaxIsometricForce(); % one leg model
        params(4,i) = muscle.getPennationAngleAtOptimalFiberLength();
    end
    params(5,:) = 10*params(2,:);
    auxdata.params = params;
    
    % Initialize the model (this builds the system and initialize the state)
    osimState = osimModel.initSystem();
    
    for k = 1:length(filenames)
        
        load(['musdyn_SRS_dlmt_Trial' num2str(trial) filenames{k}], 'solution')
        input.act =  max(solution.parameter(1:9), 0.01);
        avg_eb(tr,k) = mean(input.act);
        
        % Update model state with current values
        osimState.setTime(0);
        numVar = osimState.getNY();
        for i = 0:numVar-1
            osimState.updY().set(i, s0(i+1,1));
        end
        
        % Compute the isometric muscle length!
        
        % Get the number of states, coordinates, muscles and controls from the model;
        % in this case the number of controls equals the number of muscles
        Nstates       = osimModel.getNumStateVariables();
        Ncontrols     = osimModel.getNumControls();
        Ncoord        = osimModel.getNumCoordinates();
        % model_muscles = osimModel.getMuscles();
        % Nmuscles      = model_muscles.getSize();
        
        % Get the names of the states from the model
        stateNames = cell(Nstates,1);
        for i = 1:Nstates
            stateNames(i,1) = cell(osimModel.getStateVariableNames().getitem(i-1));
        end
        
        auxdata.NStates = Nstates;
        
        % Get the names of the controls/muscles from the model
        osimModel.computeStateVariableDerivatives(osimState); % To be able to access muscle fiber velocity.
        Actuators = osimModel.getActuators();
        % Muscles = osimModel.getMuscles();
        controlNames = cell(Ncontrols,1);
        controls_all = cell(Ncontrols,1);
        lMT = zeros(1,M);
        vMT = zeros(1,M);
        for i = 1:Ncontrols
            currentActuator = Actuators.get(i-1);
            if strcmp(currentActuator.getConcreteClassName(), 'PathActuator')
                shouldBePathAct = osimModel.getActuators().get(currentActuator.getName());
                currentPathAct = PathActuator.safeDownCast(shouldBePathAct);
                lMT(i) = currentPathAct.getLength(osimState);
                vMT(i) = currentPathAct.getLengtheningSpeed(osimState);
            end
            controls_all(i,1) = cell(currentActuator.getName());
            controlNames(i,1) = cell(currentActuator.getName());
        end
        
        % Get the names of the controls/muscles from the model
        osimModel.computeStateVariableDerivatives(osimState); % To be able to access muscle fiber velocity.
        Actuators = osimModel.getActuators();
        % Muscles = osimModel.getMuscles();
        controls_all = cell(Ncontrols,1);
        lMT = zeros(1,M);
        vMT = zeros(1,M);
        for i = 1:Ncontrols
            currentActuator = Actuators.get(i-1);
            if strcmp(currentActuator.getConcreteClassName(), 'PathActuator')
                shouldBePathAct = osimModel.getActuators().get(currentActuator.getName());
                currentPathAct = PathActuator.safeDownCast(shouldBePathAct);
                lMT(i) = currentPathAct.getLength(osimState);
                vMT(i) = currentPathAct.getLengtheningSpeed(osimState);
            end
            controls_all(i,1) = cell(currentActuator.getName());
        end
        
        t_input = [0;1];
        lMtilda_isom = zeros(M,1);
        fse_isom = zeros(M,1);
        for m=1:M
            A = [input.act(m); input.act(m)];
            LMT = [lMT(m); lMT(m)];
            VMT = [vMT(m); vMT(m)];
            
            [t,fse] = ode15s(@TendonForceOdeVecSRS,[0,1],0,[],t_input, A,LMT,VMT,params, m, Fvparam, Fpparam, Faparam, lMT(m), 0, vec_kT(k));
            
            fse_isom(m) = fse(end);
            
            FMo = ones(size(fse,1),1)*params(1,m);
            lMo = ones(size(fse,1),1)*params(2,m);
            lTs = ones(size(fse,1),1)*params(3,m);
            alphao = ones(size(fse,1),1)*params(4,m);
            
            FT = fse .* FMo;
            lTtilda = fse/35 + 1;
            lM = sqrt((lMo.*sin(alphao)).^2+(lMT(m)-lTs.*lTtilda).^2);
            lMtilda = lM./lMo;
            lMtilda_isom(m) = lMtilda(end);
        end
        
        input.lMtilda_isom = lMtilda_isom;
        
        auxdata.ksrs = vec_ksrs(k);
        auxdata.kT = vec_kT(k);
        [tFW,sFW] = ode15s(@compute_state_derivatives,[t0,tf],[s0;fse_isom],[], input, osimModel, osimState, auxdata);
        
        % Write model states to an OpenSim STO file
        StatesData.name = ['AMPC04_Trial' num2str(trial) '_' num2str(k) '_FW_States'];
        StatesData.nRows = length(tFW);
        StatesData.nColumns = Nstates + 1;
        StatesData.labels = [{'time'}; stateNames]';
        StatesData.inDegrees = false;
        StatesData.data = [tFW,sFW(:,1:Nstates)];
        writeOpenSimStatesFile(StatesData)
        
        % Create data structure for the controls file
        ControlData.name = ['AMPC04_Trial' num2str(trial) '_' num2str(k) '_FW_Controls'];
        ControlData.nRows = size(tFW);
        ControlData.nColumns = Ncontrols+1; %All the controls + time
        ControlData.inDegrees = false;
        ControlData.labels= [{'time'}; controlNames]';
        
        mass = 2000;
        acc = interp1(input.time, input.platacc, tFW)*9.81;
        controlFX = mass*acc/10000;
        ControlData.data = [tFW,sFW(:,1+Nstates:end), controlFX];
        writeOpenSimControlFile(ControlData)
        
        %
        
        positionsFW = sFW(:, 1:4);
        positionsFW(:, 2:4) = positionsFW(:, 2:4)*180/pi;
        qFW = [tFW positionsFW];
        positionsFW(:, 1) = 1000*positionsFW(:, 1);
        
        t_IK = IKsolution.data(ind_IK,1);
        s_IK = [-IKsolution.data(ind_IK,7)+delta_pelvis_tilt IKsolution.data(ind_IK,6) -IKsolution.data(ind_IK,5)];
        
        posFW = interp1(tFW, positionsFW(:, 2:4), t_IK(t_IK<tFW(end)));
        
        t50 = (t0+tf)/2;
        posFW50 = interp1(tFW, positionsFW(:, 2:4), t_IK(t_IK<t50));
        
        RMSE(k,:) = sqrt(mean((posFW - s_IK(t_IK<tFW(end),:)).^2));
        RMSE50(k,:) = sqrt(mean((posFW50 - s_IK(t_IK<t50,:)).^2));
        
    end
    
    RMSEall(tr,:,:) = RMSE;
    RMSE50all(tr,:,:) = RMSE50;
    
end

avg_avg_eb = mean(avg_eb)

avgRMSE50 = squeeze(mean(RMSE50all,1));
stdRMSE50 = squeeze(std(RMSE50all,1));

avgRMSE = squeeze(mean(RMSEall,1));
stdRMSE = squeeze(std(RMSEall,1));

subplot(4,7,[22:24])
x = 1:3;
bar_h = bar(x,avgRMSE50'); hold on;
set(gca,'XTickLabel',{'ankle' 'knee' 'hip'});
% set(gca,'YTick',[0 1 2 3])
colormap(color(1:6,:))
sb = 0.1325;
X = [x'-2.5*sb x'-1.5*sb x'-0.5*sb x'+0.5*sb x'+1.5*sb x'+2.5*sb];
h=errorbar(X,avgRMSE50',min(avgRMSE50', stdRMSE50') ,stdRMSE50','k'); set(h,'linestyle','none')
ylabel('RMSE [^o]')
axis([0.5 3.5 0 6])
set(gca,'YTick',[0 5])
box off

subplot(4,7,[26:28])
x = 1:3;
bar_h = bar(x,avgRMSE'); hold on;
set(gca,'XTickLabel',{'ankle' 'knee' 'hip'});
% set(gca,'YTick',[0 1 2 3])
colormap(color(1:6,:))
sb = 0.1325;
X = [x'-2.5*sb x'-1.5*sb x'-0.5*sb x'+0.5*sb x'+1.5*sb x'+2.5*sb];
h=errorbar(X,avgRMSE',min(avgRMSE', stdRMSE') ,stdRMSE','k'); set(h,'linestyle','none')
ylabel('RMSE [^o]')
axis([0.5 3.5 0 6])
set(gca,'YTick',[0 5])
box off

legend('no SRS - w_{opt}', 'no SRS - w_{opt} - k_T high', 'no SRS - w_{0}', 'SRS - w_{opt}', 'SRS - w_{opt} - k_T high', 'SRS  - w_{0}')

disp('S2')
avg50 = mean(avgRMSE50,2)
std50 = std(avgRMSE50, 0, 2)