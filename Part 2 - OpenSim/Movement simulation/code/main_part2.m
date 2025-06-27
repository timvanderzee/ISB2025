clear all; close all; clc

fullfilepath = which('main_part2.m');
filefolder = fullfilepath(1:end-12);

cd(filefolder)
cd ..
cd ..
cd ..

mainfolder = cd;
import casadi.*;        % Import casadi libraries

%% fit model parameters
opti = casadi.Opti();   % Initialise opti structure
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