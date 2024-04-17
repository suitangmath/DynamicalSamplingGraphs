function Example = Community5_def()

% Sui Tang -update 05-21
% param_sim.band_limit = 10; % Band limit

N = 1000;


% 
% type =1;
% 
% switch type
%     case 1
%         
%         param_graph.com_sizes = 100*ones(1,10);
%         
%     case 2
%         
%         param_graph.com_sizes =[50 105*ones(1,8) 110];
%         
%     case 3
%         
%         param_graph.com_sizes =[25 108*ones(1,8) 111];
%         
%     case 4
%         
%         param_graph.com_sizes =[17 109*ones(1,8) 111];
%     case 5
%         
%         param_graph.com_sizes =[13 109*ones(1,8) 115];
% end
% 
% 
% param_graph.com_lims = [0, cumsum(param_graph.com_sizes)];
% % Probability of connections intercluster
% param_graph.world_density = 1/N;
% % Construct graph
% G = gsp_community(N, param_graph);

%  save('community_graph1.mat','G');
load('community_graph5.mat');
%% update the weight matrix
%G = Graph_update_weights(G,0.5);

%% construct the heat diffusion matrix
D = diag(G.d);
N = size(D,1);
L = eye(N)-pinv(sqrt(D))*G.W*pinv(sqrt(D)); % normalized weighted graph Laplacian
deltat = 4;% sampling time gap
A = expm(-deltat*L);% heat diffusion matrix


% System
sysInfo.name            = 'Community5';  % name of the graph
sysInfo.N               = N;      % # of nodes
sysInfo.G = G;
sysInfo.L = L;
sysInfo.A = A;% heat diffusion matrix
sysInfo.bwidth = 100;
bwidth = sysInfo.bwidth;

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