function Example = graph_SeaSurfaceTemperature_def()

% Longxiu Huang -update 06-22-2021
load('SeaSurfaceTemperature.mat');
load('/home/huangL3/DSonGraphs-master/codes/CurrentPaper/compared_methods/Timevarying_GS_Reconstruction/experiments/SeaSurfaceTemperature/paramAWD.mat');

% System
sysInfo.name  = 'SeaSurfaceTemperature';  % name of the graph
%sysInfo.G = G;
%sysInfo.L = L;


W = zeros(100);
load('SeaSurfaceTemperature.mat')


sigma =100;

for i=1:100
    for j=1:100
        W(i,j)=exp(-norm(Position(i,:)-Position(j,:))^2/sigma);
    end
end

L = eye(100)-sqrt(pinv(diag(sum(W,2))))*W*sqrt(pinv(diag(sum(W,2))));
tau =0.4;
A = expm(-tau*L);




Temp = Data(:,1:20); % the time length is 600.
[N,T] = size(Temp);
sysInfo. L =L;
sysInfo.N = N;    % # of nodes
sysInfo.A = A;%Temp(:,2:end)*pinv(Temp(:,1:end-1));
sysInfo.Data = Temp;
sysInfo.bwidth = 20;


% Observations will be up to this time
obsInfo.tau  = T;
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