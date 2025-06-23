% function state = forward_simulation(time, input, initial_state)

close all
clear all
global DATASET

% Used dataset
% ------------
DATASET = 'SCP_simple';

% Add paths for software and dataset
% ----------------------------------
warning off
voeg_paden_toe(DATASET);
warning on

load_simulation_parameters

trials = [15 16 28 35 36 37 42];
vec_ksrs = [0 280];

RMSEall = zeros(length(trials), 4, 3);
% without SRS 50ms, without SRS 100ms, with SRS 50ms, with SRS 100ms

%                  1          2          3          4          5          6          7        8           9        10            11        12
%                red     fade red     green    fade green    blue     fade blue   black     gray       yellow    orange        plum      lime
% color_palette= [255 0 0; 255 150 150; 0 114 54; 153 255 102; 0 0 255; 44 176 207; 0 0 0; 155 155 155; 255 153 0; 255 185 0; 153 51 102; 171 218 77]*(1/255);
% color=color_palette([12 3 11 5 10 1 11 5 10 4 3 5 12 7 8 9 1 10 4 6 11 3 2 5 12],:);

color_palette= [0 114 54; 75 171 80; 153 228 108; 153 51 102; 205 105 210; 221 160 230]*(1/255);
color=color_palette([1 2 3 4 5 6],:);

figure()
for tr = 1:length(trials)
    trial = trials(tr)
    
    IKpath = ['C:\Users\u0046458\Documents\Unief\code\GT\Models\SCP_gait10dof18musc_Trial' num2str(trial) '_ks.mot'];
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
    
    
    load C:\Users\u0046458\Documents\Unief\code\GT\OSIM_posturePerts_Friedl\matlab_files\PlatOnsets.dat
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
    
    
    
    load(['C:\Users\u0046458\Documents\Unief\code\GT\OSIM_posturePerts_Friedl\processedData\Trial' num2str(trial)], 'data');
    
    platpos = data.analog.platemotion(:,2)/100;
    input.platacc = data.analog.platemotion(:,6);
    input.time = [1:1:length(data.analog.platemotion(:,6))]/data.analog.samplerate;
    
    % Import the OpenSim modeling classes
    import org.opensim.modeling.*
    
    % Read in the osim model
    osimModel = Model('C:\Users\u0046458\Documents\Unief\code\GT\Models\SCP_gait4dof9musc.osim'); %
    osimModel3D = Model('C:\Users\u0046458\Documents\Unief\code\GT\Models\SCP_gait10dof18musc.osim'); %
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
    
    for k = 1:length(vec_ksrs)
        
        if k == 1
            load(['musdyn_SRS_dlmt_Trial' num2str(trial) '_ebase5_0_50ms_ksrs' num2str(vec_ksrs(k)) '_web3'], 'solution')
            % load(['musdyn_SRS_dlmt_Trial15_ebase3_0_50ms_ksrs' num2str(vec_ksrs(k)) '_web1'], 'solution')
        elseif k==2
            load(['musdyn_SRS_dlmt_Trial' num2str(trial) '_ebase5_0_50ms_ksrs' num2str(vec_ksrs(k)) '_web3'], 'solution')
        end
        input.act =  max(solution.parameter(1:9), 0.01);
        act = input.act
        
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
            controlNames(i,1) = cell(currentActuator.getName());
        end
        
        t_input = [0;1];
        lMtilda_isom = zeros(M,1);
        fse_isom = zeros(M,1);
        for m=1:M
            A = [input.act(m); input.act(m)];
            LMT = [lMT(m); lMT(m)];
            VMT = [vMT(m); vMT(m)];
            
            [t,fse] = ode15s(@TendonForceOdeVecSRS,[0,1],0,[],t_input, A,LMT,VMT,params, m, Fvparam, Fpparam, Faparam, lMT(m), 0, 35);
            
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
        auxdata.kT = 35;
        [tSRS,sSRS] = ode15s(@compute_state_derivatives,[t0,tf],[s0;fse_isom],[], input, osimModel, osimState, auxdata);
        
        % Write model states to an OpenSim STO file
        if k == 1
            StatesData.name = ['SCP_Trial' num2str(trial) '_Hill_FW_States'];
        else
            StatesData.name = ['SCP_Trial' num2str(trial) '_SRS_FW_States'];
        end
        StatesData.nRows = length(tSRS);
        StatesData.nColumns = Nstates + 1;
        StatesData.labels = [{'time'}; stateNames]';
        StatesData.inDegrees = false;
        StatesData.data = [tSRS, sSRS(:,1:Nstates)];
        writeOpenSimStatesFile(StatesData)
        
        % Create data structure for the controls file
        if k == 1
            ControlData.name = ['SCP_Trial' num2str(trial) '_Hill_FW_Controls'];
        else
            ControlData.name = ['SCP_Trial' num2str(trial) '_SRS_FW_Controls'];
        end
        ControlData.nRows = size(tSRS);
        ControlData.nColumns = Ncontrols+1; %All the controls + time
        ControlData.inDegrees = false;
        ControlData.labels= [{'time'}; controlNames]';
        
        mass = 2000;
        acc = interp1(input.time, input.platacc, tSRS)*9.81;
        controlFX = mass*acc/10000;
        ControlData.data = [tSRS, sSRS(:,1+Nstates:end), controlFX];
        writeOpenSimControlFile(ControlData)
        
        if k == 1
            tHill = tSRS;
            sHill = sSRS;
        end
        
    end
    
    positionsHill = sHill(:, 1:4);
    positionsHill(:, 2:4) = positionsHill(:, 2:4)*180/pi;
    qHill = [tHill positionsHill];
    positionsHill(:, 1) = 1000*positionsHill(:, 1);
    velocitiesHill = sHill(:, 5:8);
    velocitiesHill(:, 2:4) = velocitiesHill(:, 2:4)*180/pi;
    fseHill = sHill(:, 9:end);
    
    positionsSRS = sSRS(:, 1:4);
    positionsSRS(:, 2:4) = positionsSRS(:, 2:4)*180/pi;
    qSRS = [tSRS positionsSRS];
    positionsSRS(:, 1) = 1000*positionsSRS(:, 1);
    velocitiesSRS = sSRS(:, 5:8);
    velocitiesSRS(:, 2:4) = velocitiesSRS(:, 2:4)*180/pi;
    fseSRS = sSRS(:, 9:end);
    
    t_IK = IKsolution.data(sampleBegin:sampleEnd,1);
    s_IK = [-IKsolution.data(sampleBegin:sampleEnd,7)+delta_pelvis_tilt IKsolution.data(sampleBegin:sampleEnd,6) -IKsolution.data(sampleBegin:sampleEnd,5)];
    
    minHill = min(positionsHill(:, 2:4));
    minSRS = min(positionsSRS(:, 2:4));
    minIK = min(s_IK);
    minQ = min([minHill; minSRS; minIK]);
    
    maxHill = max(positionsHill(:, 2:4));
    maxSRS = max(positionsSRS(:, 2:4));
    maxIK = max(s_IK);
    maxQ = max([maxHill; maxSRS; maxIK]);
    
    avgQ = (minQ + maxQ)/2;
    
    % eval(['save data/', DATASET, '/results/forward_Trial' num2str(trial) '_ebase3_0_50ms_ksrs280_web0 tHill sHill tSRS sSRS'])
    
    save(['Trial' num2str(trial) '_kinematics_FW_Hill'], 'qHill', '-ascii')
    save(['Trial' num2str(trial) '_kinematics_FW_SRS'], 'qSRS', '-ascii')
    
    posHill = interp1(tHill, positionsHill(:, 2:4), t_IK(t_IK<tHill(end)));
    posSRS = interp1(tSRS, positionsSRS(:, 2:4), t_IK(t_IK<tSRS(end)));
    
    t50 = (t0+tf)/2;
    posHill50 = interp1(tHill, positionsHill(:, 2:4), t_IK(t_IK<t50));
    posSRS50 = interp1(tSRS, positionsSRS(:, 2:4), t_IK(t_IK<t50));
    
    RMSE_Hill = sqrt(mean((posHill - s_IK(t_IK<tHill(end),:)).^2));
    RMSE_SRS = sqrt(mean((posSRS - s_IK(t_IK<tSRS(end),:)).^2));
    
    RMSE_Hill50 = sqrt(mean((posHill50 - s_IK(t_IK<t50,:)).^2));
    RMSE_SRS50 = sqrt(mean((posSRS50 - s_IK(t_IK<t50,:)).^2));
    
    RMSEall(tr,1,:) = RMSE_Hill50;
    RMSEall(tr,2,:) = RMSE_SRS50;
    RMSEall(tr,3,:) = RMSE_Hill;
    RMSEall(tr,4,:) = RMSE_SRS;
    
    % RMSE_Hill_100
    % RMSE_SRS_100
    
    % Compute fiber length
    % Need to compute LMT first.
    
    % lMo = ones(size(fse,1),1)*params(2,:);
    % lTs = ones(size(fse,1),1)*params(3,:);
    % alphao = ones(size(fse,1),1)*params(4,:);
    %
    % lTtilda = fse/35 + 1;
    % lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilda).^2);
    % lMtilda = lM./lMo;
    %
    % lMo = ones(size(fseSRS,1),1)*params(2,:);
    % lTs = ones(size(fseSRS,1),1)*params(3,:);
    % alphao = ones(size(fseSRS,1),1)*params(4,:);
    %
    % lTtilda = fseSRS/35 + 1;
    % lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilda).^2);
    % lMtildaSRS = lM./lMo;
    
    dofNames = {'platform translation','ankle plantairflexion','knee extension','hip flexion'};
    
    % Plot dofs for all trials
    % for i = 1:4
    %     subplot(length(trials),4,(tr-1)*4+i)
    %     plot(1000*(tHill-t0), positionsHill(:,i), 'k', 'LineWidth', 2); hold on;
    %     plot(1000*(tSRS-t0), positionsSRS(:,i), 'k--', 'LineWidth', 2); hold on;
    %     if tr == length(trials)
    %     xlabel('t [ms]')
    %
    %     end
    %     if i == 1
    %         tplat = input.time(input.time>t0 & input.time<tf);
    %         pplat = platpos(input.time>t0 & input.time<tf);
    %         plot(1000*(tplat-t0), pplat*1000, 'r', 'LineWidth', 2);
    %         ylabel('t_{plat} [mm]')
    %
    %     elseif i == 2
    % %         ymin = min(-IKsolution.data(sampleBegin:sampleEnd,7)+delta_pelvis_tilt);
    % %         ymax = max(-IKsolution.data(sampleBegin:sampleEnd,7)+delta_pelvis_tilt);
    % %         ymid = (ymax-ymin)/2
    %         plot(1000*(IKsolution.data(sampleBegin:sampleEnd,1)-t0), -IKsolution.data(sampleBegin:sampleEnd,7)+delta_pelvis_tilt, 'r', 'LineWidth', 2); hold on;
    % %         axis([t0 tf 0 12.5])
    %         ylabel('q [^o]')
    %         axis([0 100 avgQ(i-1)-15 avgQ(i-1)+15])
    %         set(gca,'XTick',[0 50 100])
    %     elseif i == 3
    %         plot(1000*(IKsolution.data(sampleBegin:sampleEnd,1)-t0), IKsolution.data(sampleBegin:sampleEnd,6), 'r', 'LineWidth', 2); hold on;
    %         ylabel('q [^o]')
    %         axis([0 100 avgQ(i-1)-15 avgQ(i-1)+15])
    %         set(gca,'XTick',[0 50 100])
    %     else
    %         plot(1000*(IKsolution.data(sampleBegin:sampleEnd,1)-t0), -IKsolution.data(sampleBegin:sampleEnd,5), 'r', 'LineWidth', 2); hold on;
    %         ylabel('q [^o]')
    %         axis([0 100 avgQ(i-1)-15 avgQ(i-1)+15])
    %         set(gca,'XTick',[0 50 100])
    %         if tr == 1
    %         legend('Hill', 'Hill + SRS', 'experiments')
    %         end
    %     end
    %     if tr == 1
    %     title(dofNames{i})
    %     end
    %
    % end
    
    if trial == 15 | trial == 35
        if trial == 15
            extr = 1;
        else
            extr = 2;
        end
        
        for i = 1:3
            subplot(3,4,(extr-1)*4+i)
            plot(1000*(tHill-t0), positionsHill(:,i+1), 'Color', color(1,:), 'LineWidth', 2); hold on;
            plot(1000*(tSRS-t0), positionsSRS(:,i+1), 'Color', color(3,:), 'LineWidth', 2); hold on;
            if extr == 2
                xlabel('t [ms]')
            end
            if i == 1
                %         ymin = min(-IKsolution.data(sampleBegin:sampleEnd,7)+delta_pelvis_tilt);
                %         ymax = max(-IKsolution.data(sampleBegin:sampleEnd,7)+delta_pelvis_tilt);
                %         ymid = (ymax-ymin)/2
                plot(1000*(IKsolution.data(sampleBegin:sampleEnd,1)-t0), -IKsolution.data(sampleBegin:sampleEnd,7)+delta_pelvis_tilt, 'k', 'LineWidth', 3); hold on;
                %         axis([t0 tf 0 12.5])
                ylabel('q [^o]')
                axis([0 100 avgQ(i)-15 avgQ(i)+15])
                set(gca,'XTick',[0 50 100])
                box off
            elseif i == 2
                plot(1000*(IKsolution.data(sampleBegin:sampleEnd,1)-t0), IKsolution.data(sampleBegin:sampleEnd,6), 'k', 'LineWidth', 3); hold on;
                %         ylabel('q [^o]')
                axis([0 100 avgQ(i)-15 avgQ(i)+15])
                set(gca,'XTick',[0 50 100])
                box off
            else
                plot(1000*(IKsolution.data(sampleBegin:sampleEnd,1)-t0), -IKsolution.data(sampleBegin:sampleEnd,5), 'k', 'LineWidth', 3); hold on;
                %         ylabel('q [^o]')
                axis([0 100 avgQ(i)-15 avgQ(i)+15])
                set(gca,'XTick',[0 50 100])
                box off
                if extr == 1
                    legend('Hill', 'Hill + SRS', 'experiments')
                end
            end
            if extr == 1
                title(dofNames{i+1})
            end
            
        end
        
    end
    
    
    % muscleNames = {'hamstrings'; 'biceps femoris sh'; 'gluteus maximus'; 'iliopsoas'; 'rectus femoris'; 'vasti'; 'gastrocnemius'; 'soleus'; 'tibialis anterior'};
    % figure()
    % for i = 1:auxdata.NMuscles
    %     subplot(3,3,i)
    %     plot(t, fse(:,i), 'k', 'LineWidth', 2); hold on;
    %     plot(tSRS, fseSRS(:,i), 'k--', 'LineWidth', 2)
    %     xlabel('t [s]')
    %     ylabel('F_T^n [m]')
    %     title(muscleNames{i})
    %     axis tight
    % end
    
    % muscleNames = {'hamstrings'; 'biceps femoris sh'; 'gluteus maximus'; 'iliopsoas'; 'rectus femoris'; 'vasti'; 'gastrocnemius'; 'soleus'; 'tibialis anterior'};
    % figure()
    % for i = 1:auxdata.NMuscles
    %     subplot(3,3,i)
    %     plot(t, lMtilda(:,i), 'k', 'LineWidth', 2); hold on;
    %     plot(tSRS, lMtildaSRS(:,i), 'k--', 'LineWidth', 2); hold on;
    %     plot([t(1) t(end)], [input.lMtilda_isom(i), input.lMtilda_isom(i)], 'r', 'LineWidth', 2); hold on;
    %     xlabel('t [s]')
    %     ylabel('l_M^n [m]')
    %     title(muscleNames{i})
    %     axis tight
    % end
    
end

avgRMSE = squeeze(mean(RMSEall,1));
stdRMSE = squeeze(std(RMSEall,1));

subplot(3,4,[9:11])
bar_h = bar(avgRMSE'); hold on;
set(gca,'XTickLabel',{'ankle' 'knee' 'hip'});
% set(gca,'YTick',[0 1 2 3])
colormap(color([1 3 5 6],:))
x = 1:3;
X = [x'-0.27 x'-0.09 x'+0.09 x'+0.27];
h=errorbar(X,avgRMSE',min(avgRMSE', stdRMSE') ,stdRMSE','k'); set(h,'linestyle','none')
ylabel('RMSE [^o]')
box off

legend('Hill - 50ms', 'SRS - 50ms', 'Hill - 100ms', 'SRS - 100ms')
