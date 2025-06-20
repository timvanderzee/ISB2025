function [data_exp] = ImportExperimentalData(info, bool_plot, params_subject, dt_spline)
%Load experimental data

%   1. Create name of trial you want to load
path        = info.path;
trial       = info.trial;
trialname   = [path, 'BK_Trial',num2str(trial),'.mat'];

%   2. Load data
load([trialname]);
q_exp = (data(:,2) - 90)*pi/180; % put 0 as horizontal and transform to radians ;
t_exp  = data(:,1);
t_span = [t_exp(1) t_exp(end)];

%   3. Export experimental angles in radians with 0 as horizontal
data_exp.t      = t_exp;
data_exp.q      = q_exp;
data_exp.qd     = diff(q_exp)./(t_exp(2)-t_exp(1)); 
data_exp.t_span = t_span;
data_exp.N      = length(t_exp);

% 4. Define whether trial is HR or MR (need to be known for SRS)
Nmr = params_subject.Nmr;

if trial < Nmr
    data_exp.srsbool= 1;
else
    data_exp.srsbool = 0;
end

% 5. Plot experimental data
if bool_plot == 1
    plot(data_exp.t,data_exp.q,'k','LineWidth',1.5)
end

% 6. Spline data 
[q_spline, qd_spline, N_spline, tvect_spline] = ExpDatAtDiscrTime(t_span,t_exp,q_exp, dt_spline);      

data_exp.tspline  = tvect_spline;
data_exp.qspline  = q_spline;
data_exp.qdspline = qd_spline;
data_exp.Nspline  = N_spline;

% 7. Calculate offset between IK and BK (used to calculate LMT and MA)
% x (BK) + offset = knee angle die nodig is voor lMT en ma te bepalen
offset          = CalculateOffset(info, params_subject, t_span); 
data_exp.offset = offset*pi/180; 
end