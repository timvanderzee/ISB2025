%% Assignment 0: Preparation
clear all; close all; clc

filename = 'main_part2.m';
fullfilepath = which(filename);
filefolder = fullfilepath(1:end-length(filename));

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

% load SRS data
cd([mainfolder, '\Part 2 - OpenSim\Movement simulation\input\common'])
load('thix_data.mat','SRSdata')

% initialise opti structure
opti = casadi.Opti(); 

% define weigth vector
w1 = 10;    % weight for fitting force-velocity
w2 = 100;   % weight for fitting short-range stiffness
w3 = 1; 	% weight for regularization
w = [w1 w2 w3];

% specify parameters to be fitted
optparms = {'f', 'k11', 'k22', 'k21'};

% run the fitting function
if ishandle(1), close(1); end; figure(1)
[newparms, out] = fit_model_parameters(opti, optparms, w, SRSdata, parms);

%% Assignment 2: compare the rate functions
if ishandle(2), close(2); end; figure(2)
compare_rate_funcs(newparms, parms)

%% Assignment 2: evaluate model force-velocity and SRS properties
if ishandle(3), close(3); end; figure(3)
eval_model_behavior(newparms, out, SRSdata)

%% Assignment 3: simulate standing balance
if ishandle(4), close(4); end; figure(4)
simulate_movement('standing_balance', {'reguluar'}, [], newparms, mainfolder);

%% Assignment 4: simulate pendulum test
if ishandle(5), close(5); end; figure(5)
act = 1e-1; % activation
simulate_movement('pendulum_test', {'reguluar', 'pre_movement'}, act, newparms, mainfolder);

%% Assignment 5: simulate pendulum test with different properties
if ishandle(6), close(6); end; figure(6)
act = 1e-1; % activation
simulate_movement('pendulum_test', {'reguluar', 'pre_movement'}, act, newparms, mainfolder);


