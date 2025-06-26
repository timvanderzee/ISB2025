function[input, s0, t0, tf, models, muscle_names] = load_data(type, parms)

if strcmp(type, 'standing_balance')

    trial = 15;
    fs = 120;

    %% Initial state from inverse kinematics
    load('PlatOnsets.dat','PlatOnsets')
    eval(['load Trial' num2str(trial) '_InputDataTime'])
    tOnset = PlatOnsets(trial);
    time = InputDataTime;

    indeces = find(time > tOnset);
    sampleBegin = indeces(1) - 1 + round(fs * 0.02); % start simulation 100 ms before perturbation
    t0 = InputDataTime(sampleBegin);
    sampleEnd = indeces(1) - 1 + round(fs * 0.12); % end simulation 50 ms after platform onset
    tf = InputDataTime(sampleEnd);

    IKsolution = importdata(['SCP_gait10dof18musc_Trial' num2str(trial) '_ks.mot']);

    % INPUT
    pelvis_tilt_exp = IKsolution.data(sampleBegin,2);
    pelvis_tilt_mod = IKsolution.data(sampleBegin,7) + IKsolution.data(sampleBegin,6) +IKsolution.data(sampleBegin,5);
    delta_pelvis_tilt = pelvis_tilt_exp + pelvis_tilt_mod;

    s0 = zeros(8,1);
    s0(2:4) = [-IKsolution.data(sampleBegin,7)+delta_pelvis_tilt IKsolution.data(sampleBegin,6) -IKsolution.data(sampleBegin,5)]/180*pi;

    %% Force data and muscle excitation
    load(['Trial' num2str(trial)], 'data');

    platpos = data.analog.platemotion(:,2)/100;
    input.platacc = data.analog.platemotion(:,6);
    input.time = [1:1:length(data.analog.platemotion(:,6))]/data.analog.samplerate;

    % Muscle excitation
    load('musdyn_SRS_dlmt_Trial15_ebase3_min10_50ms_ksrs0_web30_rev.mat', 'solution')
    input.act =  max(solution.parameter(1:9), 0.01);

    models = {'SCP_gait4dof9musc.osim', 'SCP_gait10dof18musc.osim'};
    muscle_names = {'hamstrings_r'; 'bifemsh_r'; 'glut_max_r'; 'iliopsoas_r'; 'rect_fem_r'; 'vasti_r'; 'gastroc_r'; 'soleus_r'; 'tib_ant_r'};

elseif strcmp(type, 'pendulum_test')
    s0 = [];
    t0 = 0;
    tf = 1;

    models = {'leg39_path_actuators.osim', 'leg39_welded.osim'};
    muscle_names = {'bifemlh_r'; 'bifemsh_r'; 'glut_max2_r'; 'rect_fem_r'; 'vas_int_r'; 'med_gas_r'; 'soleus_r'; 'tib_ant_r'};
    
    input.act = parms.act * ones(1, length(muscle_names));
end

end