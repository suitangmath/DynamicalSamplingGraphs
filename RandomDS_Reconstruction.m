% Longxiu Huang, Deanna Needell, and Sui Tang. Robust recovery of bandlimited graph signals via randomized dynamical sampling Information and Inference: A Journal of the IMA, 2024.

% This code is a demo code for test RIP propoerty of the dynamical sampling matrix on community graph
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
%%% choose examples
Example           =  Community1_def(); % other choice include Bunny_def, Sensor_def, etc in Exampls folder
sysInfo           =  Example.sysInfo;
obsInfo           =  Example.obsInfo;
samplingInfo      =  Example.samplingInfo;
obsInfo.VERBOSE   =  VERBOSE;
obsInfo.SAVE_DIR  =  SAVE_DIR;




%% Construct 10 bandlimited signals
N = sysInfo.N; % number of nodes
tau = obsInfo.tau; % observation times 
bwidth = sysInfo.bwidth;% bandwidth of signal 
[V,D] = eigs(sysInfo.A,sysInfo.N); % get eigenbasis and eigenvalues

trials =10; % number of learning trials 
result =cell(trials,1);

% total number of space-time samples
samples = 200;

sigmaset = 0;% noise 
%sigmaset= 5*1e-2;
gammaset =10.^(-3:0.2:2); % regularization parameters


errors1 = zeros(length(sigmaset),length(gammaset),3); % reconstruction error for regime 1
errors2 = zeros(length(sigmaset),length(gammaset),3);% reconstruction error for regime 2
errors3 = zeros(length(sigmaset),length(gammaset),3); % reconstruction error for regime 3

    
    
obsInfo.total_samples=samples;
samplingInfo.opt = 1;   %% optimal sampling distribution (samplingInfo.opt=1) uniform sampling distribution (sampligInfo.opt = 0) 
if samplingInfo.opt == 1
    samplingInfo = construct_opt_sampling_distribution(sysInfo,obsInfo,samplingInfo);
end

for k=1:trials
    rng('shuffle');
    f=V(:,1:bwidth)*randn(bwidth,1); % a bandlimited signal with randomized gaussian coefficients

    %% Construct sampling operators
    samplingInfo = construct_sampling(sysInfo,obsInfo,samplingInfo);

    %% noisy samples

    for s =1:length(sigmaset)
        y1 = samplingInfo.M1*samplingInfo.PiAL*f;
        noise = sigmaset(s)*randn(size(y1));
        y1 = y1+noise;
        y2 = samplingInfo.M2*samplingInfo.PiAL*f+noise;
        y3 = samplingInfo.M3*samplingInfo.PiAL*f+noise;

        B1 =samplingInfo.PiAL'*samplingInfo.M1'*sqrt(samplingInfo.W1)*pinv(samplingInfo.Pomega1)*sqrt(samplingInfo.W1)*samplingInfo.M1*samplingInfo.PiAL;
        B2 =samplingInfo.PiAL'*samplingInfo.M2'*sqrt(samplingInfo.W2)*pinv(samplingInfo.Pomega2)*sqrt(samplingInfo.W2)*samplingInfo.M2*samplingInfo.PiAL;
        B3 =samplingInfo.PiAL'*samplingInfo.M3'*sqrt(samplingInfo.W3)*pinv(samplingInfo.Pomega3)*sqrt(samplingInfo.W3)*samplingInfo.M3*samplingInfo.PiAL;

        for g =1:length(gammaset)
            B1_reg = B1+ gammaset(g)*(sysInfo.L)^4;
            B2_reg = B2+ gammaset(g)*(sysInfo.L)^4;
            B3_reg = B3+ gammaset(g)*(sysInfo.L)^4;


            f1_reg = pinv(B1_reg)*samplingInfo.PiAL'*samplingInfo.M1'*sqrt(samplingInfo.W1)*pinv(samplingInfo.Pomega1)*sqrt(samplingInfo.W1)*y1;
            alpha_f1_reg = V(:,1:bwidth)*V(:,1:bwidth)'*f1_reg;
            beta_f1_reg = f1_reg - alpha_f1_reg;

            f2_reg = pinv(B2_reg)*samplingInfo.PiAL'*samplingInfo.M2'*sqrt(samplingInfo.W2)*pinv(samplingInfo.Pomega2)*sqrt(samplingInfo.W2)*y2;
            alpha_f2_reg = V(:,1:bwidth)*V(:,1:bwidth)'*f2_reg;
            beta_f2_reg = f2_reg - alpha_f2_reg;

            f3_reg = pinv(B3_reg)*samplingInfo.PiAL'*samplingInfo.M3'*sqrt(samplingInfo.W3)*pinv(samplingInfo.Pomega3)*sqrt(samplingInfo.W3)*y3;
            alpha_f3_reg = V(:,1:bwidth)*V(:,1:bwidth)'*f3_reg;
            beta_f3_reg = f3_reg - alpha_f3_reg;

            % errors with regularization
            errors1 (s,g,1)= errors1 (s,g,1)+norm(f1_reg-f,2)./norm(f,2);
            errors1 (s,g,2)= errors1 (s,g,2)+norm(alpha_f1_reg-f,2)./norm(f,2);
            errors1 (s,g,3)=errors1 (s,g,3)+norm(beta_f1_reg,2)./norm(f,2);

            errors2 (s,g,1)= errors2 (s,g,1)+norm(f2_reg-f,2)./norm(f,2);
            errors2 (s,g,2)=errors2 (s,g,2)+norm(alpha_f2_reg-f,2)./norm(f,2);
            errors2 (s,g,3)=errors2 (s,g,3)+norm(beta_f2_reg,2)./norm(f,2);

            errors3 (s,g,1) = errors3 (s,g,1)+norm(f3_reg-f,2)./norm(f,2);
            errors3 (s,g,2) = errors3 (s,g,2)+norm(alpha_f3_reg-f,2)./norm(f,2);
            errors3 (s,g,3) = errors3 (s,g,3)+norm(beta_f3_reg,2)./norm(f,2);
        end
    end


end

errors1= errors1/trials;
errors2= errors2/trials;
errors3= errors3/trials;

%% start to plot results 
line_width = 1.5;
marker_size = 6;
font_size = 30;
width = 680;
heigth = 460;

figure
plot(log10(errors1(:,:,1)),'ro--','LineWidth',line_width,'MarkerSize',marker_size);
hold on
plot(log10(errors2(:,:,1)),'m+:','LineWidth',line_width,'MarkerSize',marker_size);
plot(log10(errors3(:,:,1)),'bs-','LineWidth',line_width,'MarkerSize',marker_size);
legend('Regime 1','Regime 2','Regime 3')
if samplingInfo.opt == 0
    title('Uniform sampling')
else
    title('Optimal sampling')
end
hold off


% % 
% if samplingInfo.opt
%     filename=strcat('RecOpt',sysInfo.name,'m',num2str(samples(end)),'k',num2str(bwidth),num2str(time_stamp),'.mat');
% else
%     filename=strcat('RecUnif',sysInfo.name,'m',num2str(samples(end)),'k',num2str(bwidth),num2str(time_stamp),'.mat');
% end
% fullFileName = fullfile(SAVE_DIR, filename);
% save(fullFileName,'sysInfo','trials','samples','errors1','errors2','errors3','-v7.3');

 



