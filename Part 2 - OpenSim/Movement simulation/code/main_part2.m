clear all; close all; clc

fullfilepath = which('main_part2.m');
filefolder = fullfilepath(1:end-12);

cd(filefolder)
cd ..
cd ..
cd ..

mainfolder = cd;
addpath(genpath(mainfolder));
import casadi.*;        % Import casadi libraries

%% Assignment 1: simulate standing balance
% load some parameters to start with
cd([mainfolder, '\Part 2 - OpenSim\Movement simulation\input\common'])
load('parms.mat','parms')

if ishandle(1), close(1); end; figure(1)
simulate_movement('standing_balance', {'reguluar'}, [], parms, mainfolder);

%% Assignment 2: fit model parameters
% load some parameters to start with
cd([mainfolder, '\Part 2 - OpenSim\Movement simulation\input\common'])
load('new_parms.mat','parms')
load('thix_data.mat','xm','ym','ys')

% initialise opti structure
opti = casadi.Opti(); 

% define desired force-velocity properties
vmax = 5; % maximal shortening velocity [L0/s]

% define desired SRS properties
RT = 0.1; % recovery time (s)
SRS_rel = 0.7; % relative SRS compared to no pre-movement (-) for the above recovery time
V_rel = 0.5; % relative velocity at which SRS is evaluated (relative to vmax) 

% define weigth vector
w = [10 100 1]; % weights for fitting (1) force-velocity, (2) SRS, (3) regularization

% specify parameters to be fitted
optparms = {'f', 'k11', 'k22', 'k21', 'JF'};

% run the fitting function
if ishandle(1), close(1); end; figure(1)
[newparms, out] = fit_model_parameters(opti, optparms, w, vmax, RT, SRS_rel, V_rel, parms);

%% Assignment 2: evaluate model force-velocity and SRS properties
if ishandle(2), close(2); end; figure(2)
eval_model_behavior(vmax, RT, SRS_rel, V_rel, newparms, out)

% need to include this in figure and need to base SRS_rel on this data
subplot(122)
errorbar(xm, ym, ys, 'o')

%% Assignment 3: simulate standing balance
if ishandle(3), close(3); end; figure(3)
simulate_movement('standing_balance', {'reguluar'}, [], newparms, mainfolder);

%% Assignment 4: simulate pendulum test
if ishandle(4), close(4); end; figure(4)
act = 1e-2; % activation
simulate_movement('pendulum_test', {'reguluar', 'pre_movement'}, act, newparms, mainfolder);

%% Assignment 5: simulate pendulum test with different properties
if ishandle(5), close(5); end; figure(5)
act = 1e-1; % activation
simulate_movement('pendulum_test', {'reguluar', 'pre_movement'}, act, newparms, mainfolder);


