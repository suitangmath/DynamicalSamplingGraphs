function Example = Mineosta_Road_def()

% Sui Tang -update 05-21


% G =gsp_minnesota;

%% update the weight matrix
% G = Graph_update_weights(G,0.8);

%% construct the heat diffusion matrix
% D = diag(G.d);
% N = size(D,1);
% L = eye(N)-pinv(sqrt(D))*G.W*pinv(sqrt(D)); % normalized weighted graph Laplacian
% deltat = 30;% sampling time gap
% A = expm(-deltat*L);% heat diffusion matrix

% save('minesota.mat','G','A','L');
load('minesota.mat');

% System
sysInfo.name            = 'Monesota';  % name of the graph
N = size(A,1);
sysInfo.N               = N;      % # of nodes
sysInfo.G = G;
sysInfo.L = L;
sysInfo.A = A;% heat diffusion matrix
sysInfo.bwidth = 100;



% Observations will be up to this time
obsInfo.tau  = 1;
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