function Example = graph_Lac_d_Emosson_def()

% Longxiu Huang -update 06-22-2021
load('graph_Lac_d_Emosson2.mat');

% System
sysInfo.name  = 'graph_graph_Lac_d_Emosson';  % name of the graph
N = size(L,1);
sysInfo.N = N;      % # of nodes
%sysInfo.G = G;
sysInfo.L = L;
deltat = 1;% sampling time gap
sysInfo.A =  expm(-deltat*L);% heat diffusion matrix
sysInfo.bwidth = 1000;
sysInfo.X = X;



% Observations will be up to this time
obsInfo.tau  = 5;
%obsInfo.noise       = 1e-3;
obsInfo.total_samples = ceil(N);


samplingInfo.opt = 0;

% Define the sampling distribution, the defaul is the uniform
samplingInfo.p1 = ones(1,N)./N; % random fixed sensors

%random sensors at consective time

samplingInfo.p2 =[];

for l=1:obsInfo.tau
    samplingInfo.p2 = [samplingInfo.p2 ones(1,N)./N];
end

samplingInfo.p3 = ones (1,N*obsInfo.tau)./(N*obsInfo.tau);% random space-time sensors

% if we use the optimal distribution
if samplingInfo.opt
    [samplingInfo] =  construct_opt_sampling_distribution(sysInfo,obsInfo,samplingInfo);
end

% package the data
Example.sysInfo         = sysInfo;
Example.obsInfo         = obsInfo;
Example.samplingInfo    = samplingInfo;

end