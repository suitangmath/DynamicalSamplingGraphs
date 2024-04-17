% This code is for test RIP propoerty of the dynamical sampling matrix
% Last-update  Sui Tang 05-21-2021


addpath(genpath('..'))
if ispc
    SAVE_DIR =  [getenv('USERPROFILE'), '\DataAnalyses\LearningGraphSignal'];
else
    SAVE_DIR = [getenv('HOME'), '/DataAnalyses/LearningGraphSignal'];
end
% Please keep this fixed, simply create a symlink ~/DataAnalyses pointing wherever you like
VERBOSE                         = 1;                                                                % indicator to print certain output
time_stamp                      = datestr(now, 30);
if ~exist('Params','var')
    Params = [];
end
if ~exist(SAVE_DIR,'dir')
    mkdir(SAVE_DIR);
end

%% Get example parameters
Example                        = Community2_def();
sysInfo                        = Example.sysInfo;
obsInfo                        = Example.obsInfo;
samplingInfo                   = Example.samplingInfo;
obsInfo.VERBOSE                = VERBOSE;
obsInfo.SAVE_DIR               = SAVE_DIR;




%% Construct bandlimited signal
N = sysInfo.N;
tau = obsInfo.tau;
bwidth = sysInfo.bwidth;
[V,D] = eigs(sysInfo.A,sysInfo.N); % get eigenbasis and eigenvalues

trials =250;
result =cell(trials,1);

% total number of samples
samples = 20:20:500;
prob1 =zeros(length(samples),1);
prob2 = prob1;
prob3=  prob1;
Delta_rips1 = zeros(length(samples),trials);
Delta_rips2 = zeros(length(samples),trials);
Delta_rips3 = zeros(length(samples),trials);

for j=1:length(samples)
    
    
    obsInfo.total_samples=samples(j);
    
    for k=1:trials
        rng('shuffle')
        
        %% Construct sampling operators
        samplingInfo = construct_sampling(sysInfo,obsInfo,samplingInfo);
        
        
        
        B1 =samplingInfo.PiAL'*samplingInfo.M1'*sqrt(samplingInfo.W1)*pinv(samplingInfo.Pomega1)*sqrt(samplingInfo.W1)*samplingInfo.M1*samplingInfo.PiAL;
        Delta_rip1 = min(svds(V(:,1:bwidth)'*B1*V(:,1:bwidth),bwidth));
        Delta_rips1(j,k) = Delta_rip1;
        
        
        B2 =samplingInfo.PiAL'*samplingInfo.M2'*sqrt(samplingInfo.W2)*pinv(samplingInfo.Pomega2)*sqrt(samplingInfo.W2)*samplingInfo.M2*samplingInfo.PiAL;
        Delta_rip2  = min(svds(V(:,1:bwidth)'*B2*V(:,1:bwidth),bwidth));
        Delta_rips2(j,k) = Delta_rip2;
        
        
        B3 =samplingInfo.PiAL'*samplingInfo.M3'*sqrt(samplingInfo.W3)*pinv(samplingInfo.Pomega3)*sqrt(samplingInfo.W3)*samplingInfo.M3*samplingInfo.PiAL;
        Delta_rip3  = min(svds(V(:,1:bwidth)'*B3*V(:,1:bwidth),bwidth));
        Delta_rips3(j,k) = Delta_rip3;
        
        prob1 (j)= prob1(j) + double(Delta_rip1 >0.005)./trials;
        prob2 (j)= prob2(j) + double(Delta_rip2 >0.005)./trials;
        prob3(j) = prob3(j) + double(Delta_rip3>0.005)./trials;
        
    end
end


if samplingInfo.opt
    filename=strcat('Opt',sysInfo.name,'m',num2str(samples(end)),'k',num2str(bwidth),num2str(time_stamp),'.mat');
else
    filename=strcat('Unif',sysInfo.name,'m',num2str(samples(end)),'k',num2str(bwidth),num2str(time_stamp),'.mat');
end
fullFileName = fullfile(SAVE_DIR, filename);
%save(fullFileName,'sysInfo','trials','samples','prob1','prob2','prob3','-v7.3');
save(fullFileName,'sysInfo','trials','samples','prob1','prob2','prob3','Delta_rips1','Delta_rips2','Delta_rips3','-v7.3');
