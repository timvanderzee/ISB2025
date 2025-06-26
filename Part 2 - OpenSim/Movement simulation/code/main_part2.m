clear all; close all; clc

mainfolder = 'C:\Users\u0167448\Documents\GitHub\ISB2025';
casadifolder = 'C:\Users\u0167448\Documents\GitHub\casadi-windows-matlabR2016a-v3.5.5';

addpath(genpath(mainfolder))
addpath(genpath(casadifolder))

import casadi.*;        % Import casadi libraries
opti = casadi.Opti();   % Initialise opti structure

%% fit model parameters
parms = fit_model_parameters(opti);

%% simulate movement
act = .05;
[FSE] = simulate_movement(act, parms, mainfolder);

%% first sing excursion
figure(4)
bar(categorical({'Hill', 'Biophysical'}), -FSE)
legend('No pre-movement', 'With pre-movement', 'location', 'bestoutside')
legend boxoff
ylabel('First swing excursion (rad)')
box off