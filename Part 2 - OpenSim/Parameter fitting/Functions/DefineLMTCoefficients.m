function [coeff_LMT_ma] = DefineLMTCoefficients(params_subject, info, muscles, bool_plot)

map_muscleanal        = [info.path 'MA_FakeMot_T',num2str(info.trial)];

%  1. Fakemot file, muscle analysis, data
% Theta
theta = (-140:0.5:20)' *pi/180; theta_sq = theta.^2; theta_th = theta.^3; theta_fo = theta.^4;

% Muscle analysis
name_LMT = [map_muscleanal,'_MuscleAnalysis_Length.sto'];
if params_subject.leg == 18
    name_MA = [map_muscleanal,'_MuscleAnalysis_MomentArm_knee_angle_l.sto'];
else
    name_MA = [map_muscleanal,'_MuscleAnalysis_MomentArm_knee_angle_r.sto'];
end

% LMT data
LMT_data  = importdata([name_LMT]);
if params_subject.leg == 18
    muscles_ind{1} = [char(muscles(1)),'l'];
    muscles_ind{2} = [char(muscles(2)),'l'];
else
    muscles_ind{1} = [char(muscles(1)),'r'];
    muscles_ind{2} = [char(muscles(2)),'r'];
end
col_ext        = find(strcmp(muscles_ind{1},LMT_data.colheaders));
col_flex       = find(strcmp(muscles_ind{2},LMT_data.colheaders));

LMT(:,1)       = LMT_data.data(:,col_ext);
LMT(:,2)       = LMT_data.data(:,col_flex);

% MA data
MA_data   = importdata([name_MA]);
if params_subject.leg == 18
    muscles_ind{1} = [char(muscles(1)),'l'];
    muscles_ind{2} = [char(muscles(2)),'l'];
else
    muscles_ind{1} = [char(muscles(1)),'r'];
    muscles_ind{2} = [char(muscles(2)),'r'];
end
col_ext        = find(strcmp(muscles_ind{1},MA_data.colheaders));
col_flex       = find(strcmp(muscles_ind{2},MA_data.colheaders));

MA(:,1)        = MA_data.data(:,col_ext);
MA(:,2)        = MA_data.data(:,col_flex);

%  2. Define LMT and MA coeefficients  - EXTENSOR
one      = ones(length(LMT),1);
zero     = zeros(length(MA),1);

A        = [one  theta theta_sq   theta_th ];
B        = [zero one   2*theta    3*theta_sq ];
C        = [A; B];

coeff_LMT_MA = zeros(4,2);

for i = 1:2
    LMT_ma_dat = [LMT(:,i); - MA(:,i)];
    coeff_LMT_MA(:,i) = C\LMT_ma_dat;
    
    res_LMT = sqrt(mean((A*coeff_LMT_MA(:,i) - LMT(:,i)).^2))
    res_ma  = sqrt(mean((-B*coeff_LMT_MA(:,i) - MA(:,i)).^2))
    
    if bool_plot == 1
        figure()
        subplot(211)
        plot(theta,LMT(:,i),'Color',[0.5 0.5 0.5],'LineWidth', 3); hold on;
        plot(theta,A*coeff_LMT_MA(:,i),'k','LineWidth', 1);
        ylabel('lMT')
        title('extensor')
        subplot(212)
        plot(theta,-MA(:,i),'Color',[0.5 0.5 0.5],'LineWidth', 3); hold on;
        plot(theta,B*coeff_LMT_MA(:,i),'k','LineWidth', 1);
        xlabel('theta')
        ylabel('moment arm')
        legend('input','fit')
    end
    
end
coeff_LMT_ma_ext  = coeff_LMT_MA(:,1);
coeff_LMT_ma_flex = coeff_LMT_MA(:,2);
coeff_LMT_ma      = [coeff_LMT_ma_ext coeff_LMT_ma_flex]; 

end