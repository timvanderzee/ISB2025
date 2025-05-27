function run_ramp(parms,sarcE)

% pCatoAct=load(strcat(parms.fp_custMyoSim,'Results/pCaActCurve/',num2str(parms.date),'/SimNo',num2str(parms.simNo),'pCaActCurve.mat'));
% 
% sarcE.pCa = interp1(pCatoAct.fracBoundNormE,pCatoAct.pCaList,sarcE.initial_act);
sarcE.pCa = 4.5;

ramp_vel = ((45)/100)*1300;
ramp_amp = ((7)/100)*1300;
ramp_rise_time_s = (ramp_amp/ramp_vel);
ramp_hold_time = 2;
pertStart = 2*1/parms.time_step;

no_of_ramp_steps = round(ramp_rise_time_s /parms.time_step);
dx = ramp_amp / no_of_ramp_steps;

hs_lengths = sarcE.hs_length;

if parms.generateResults
    tic
    delta_cdl=zeros(1,round(pertStart));

    delta_cdl = [delta_cdl , dx * ones(1,no_of_ramp_steps)];
    delta_cdl = [delta_cdl , 0 * ones(1,(ramp_hold_time/ parms.time_step))];

    t = [0:parms.time_step:(numel(delta_cdl)-1)*parms.time_step];
    sarcE.pCa = sarcE.pCa*ones(size(t));
    
%     actVec = sarcE.initial_act*ones(size(t));
%     sarcE.pCa = interp1(pCatoAct.fracBoundNormE,pCatoAct.pCaList,actVec/100);

    [hs,data] = musTenDriverForVdz4State(t,0,delta_cdl,sarcE);

    save(strcat(parms.fp_custMyoSim,'Results/ramp/',num2str(parms.date),...
        '/results_rampL0',num2str(hs_lengths),'_Act',num2str(sarcE.initial_act),'SimNo',num2str(parms.simNo),'.mat'),'hs','data')
    toc; beep;
end

if parms.plotFigs
    figure(8);hold on; grid on;
    ax(1)=subplot(211);hold on; grid on;
    ylabel('sarcomere length [nm]');
    title('3StateCustMyoSimWithCoopSRX for ramp SRS test')
    ax(2)=subplot(212);hold on; grid on;
    xlabel('time [s]');ylabel('sarcomere stresas [Nm^{-2}]')

        results_file = strcat(parms.fp_custMyoSim,'Results/ramp/',num2str(parms.date),...
        '/results_rampL0',num2str(hs_lengths),'_Act',num2str(sarcE.initial_act),'SimNo',num2str(parms.simNo),'.mat');

        % Load the simulation back in
        sim_output = load(results_file);

        figure(8);subplot(211);plot(sim_output.data.t,sim_output.data.hs_length);figurefyTalk
        subplot(212);plot(sim_output.data.t,sim_output.data.hs_force);figurefyTalk
        %         pause
end