function Example = Erdos_renyi_def()

% Sui Tang -update 05-21


N = 300;
G = gsp_erdos_renyi(N,0.4);


%% update the weight matrix
G = Graph_update_weights(G,0.5);

%% construct the heat diffusion matrix
D = diag(G.d);
L = eye(N)-pinv(sqrt(D))*G.W*pinv(sqrt(D)); % normalized weighted graph Laplacian
deltat = 0.2;% sampling time gap
A = expm(-deltat*L);% heat diffusion matrix


% System
sysInfo.name            = 'Rand_erdos_renyi';  % name of the graph
sysInfo.N               = N;      % # of nodes
sysInfo.G = G;
sysInfo.L = L;
sysInfo.A = expm(-deltat*L);% heat diffusion matrix



% Observations will be up to this time
obsInfo.tau  = 40;                                                                   
obsInfo.noise       = 1e-3;
obsInfo.total_samples = 300;


% Define the sampling distribution, the defaul is the uniform 
samplingInfo.p1 = ones(1,N)./N; % random fixed sensors 
samplingInfo.m1 = ceil(obsInfo.total_samples/obsInfo.tau);

%random sensors at consective time

samplingInfo.p2 = [];
samplingInfo.m2=[];
for l=1:obsInfo.tau
samplingInfo.p2 = [samplingInfo.p2 ones(1,N)/N];
samplingInfo.m2 = [samplingInfo.m2 ceil(obsInfo.total_samples/obsInfo.tau)];
end

samplingInfo.p3 = ones (1,N*obsInfo.tau)./(N*obsInfo.tau);% random space-time sensors

% package the data
Example.sysInfo         = sysInfo;
Example.obsInfo         = obsInfo;
Example.samplingInfo    = samplingInfo;

end