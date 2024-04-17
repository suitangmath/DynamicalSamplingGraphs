function Example = Sensor_def()

% Sui Tang -update 05-21


N = 200;
G = gsp_sensor(N);


% update the weight matrix
G = Graph_update_weights(G,0.5);

% construct the heat diffusion matrix
D = diag(G.d);
L = eye(N)-pinv(sqrt(D))*G.W*pinv(sqrt(D)); % normalized weighted graph Laplacian
deltat = 20;% sampling time gap
A = expm(-deltat*L);% heat diffusion matrix
save('sensor200.mat','G','D','L','A');
load('sensor200.mat');

% System
sysInfo.name            = 'Sensor';  % name of the graph
sysInfo.N               = N;      % # of nodes
sysInfo.G = G;
sysInfo.L = L;
sysInfo.A = A;% heat diffusion matrix
sysInfo.bwidth = 10;


% Observations will be up to this time
obsInfo.tau  = 1;                                                                   
obsInfo.noise       = 0;
obsInfo.total_samples = N;

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