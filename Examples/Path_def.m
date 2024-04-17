function Example = Path_def()

% Sui Tang -update 05-21


N = 50;
% G = gsp_path(N);
% 
% 
% %% update the weight matrix
% % G = Graph_update_weights(G,0.8);
% 
% %% construct the heat diffusion matrix
% D = diag(G.d);
% L = eye(N)-pinv(sqrt(D))*G.W*pinv(sqrt(D)); % normalized weighted graph Laplacian
% deltat = 10;% sampling time gap
% A = expm(-deltat*L);% heat diffusion matrix
% 
% 
% save('path_graph50.mat','G','A','L');
load('path_graph50.mat');

% System
sysInfo.name            = 'Path';  % name of the graph
sysInfo.N               = N;      % # of nodes
sysInfo.G = G;
sysInfo.L = L;
sysInfo.bwidth = 10;
sysInfo.A =A;% heat diffusion matrix


% Observations will be up to this time
obsInfo.tau  = 20;
obsInfo.noise       = 1e-3;
obsInfo.total_samples = 300;

samplingInfo.opt =0;
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