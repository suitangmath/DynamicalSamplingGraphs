function [samplingInfo]=construct_sampling(sysInfo,obsInfo,samplingInfo)
% Longxiu Huang, Deanna Needell, and Sui Tang. Robust recovery of bandlimited graph signals via randomized dynamical sampling Information and Inference: A Journal of the IMA, 2024.

% Construct sampling operators associated with three sampling regimes  
N = sysInfo.N;
A = sysInfo.A;
tau = obsInfo.tau;
total = obsInfo.total_samples;
noise = obsInfo.noise*randn(N*tau,1);% add possible noise 
I = eye(N);


% embedding operator 
PiAL= zeros (N*tau,N);
for l=1:obsInfo.tau
     PiAL((l-1)*N+1:l*N,:) = A^(l-1);
end



%% construct sampling operators for random fixed sensors
 m1 = ceil(obsInfo.total_samples/obsInfo.tau);
 M1 = zeros(tau*m1,N*tau);
 W1 = eye(tau*m1,tau*m1)./m1;
 Pomega1 = zeros (tau*m1,tau*m1);
 p1 = samplingInfo.p1;
 noise1 = zeros(m1*tau,1);

 T1=randsample(N,m1,true,p1);
%  T1 = gendist(p1,1,m1); % random samplimg set at t = 0 and then fix it 
 Omega1 = zeros(N,1);
 Omega1(unique(T1))=1;
 
 for l =1:tau
     
 M1((l-1)*m1+1:l*m1,(l-1)*N+1:l*N) = I(T1,:);
 Pomega1((l-1)*m1+1:l*m1,(l-1)*m1+1:l*m1) = diag(p1(T1));
 noise1(m1*(l-1)+1:l*m1) = noise((l-1)*N+T1); 
 
 end
 
 
%% construct sampling operator for random sensors at consective times 

m2=zeros(obsInfo.tau,1);
for l=1:obsInfo.tau
  m2(l) =ceil(obsInfo.total_samples/obsInfo.tau);
end

p2 = samplingInfo.p2;
M2 = zeros(sum(m2),N*tau);
W2 = eye(sum(m2),sum(m2));
Pomega2 = zeros(sum(m2),sum(m2));
noise2 = zeros(sum(m2),1);
rows=0;
  
 for l=1:tau
     T2 = randsample(N,m2(l),true,p2((l-1)*N+1:l*N));
     M2(1+rows:rows+m2(l),(l-1)*N+1:l*N) = I(T2,:);
     W2(1+rows:rows+m2(l),:)= W2(1+rows:rows+m2(l),:)./m2(l);
     Pomega2(1+rows:rows+m2(l),1+rows:rows+m2(l)) = diag(p2(rows+T2));
     noise2(rows+1:rows+m2(l)) = noise(rows+T2); 
     rows = rows+m2(l);
     
 end
 

 %% construct sampling operator for random space-time samples 

 p3 = samplingInfo.p3;
 W3 = eye(total)./total;
%  T3 = gendist(p3,1,total);

T3 = randsample(N*tau,total,true,p3);
 
 M3 = zeros(length(T3),N*tau);
 for i=1:length(T3)
     M3(i,T3(i))=1;
 end
 
 Pomega3 = diag(p3(T3));
 noise3 = noise(T3);

 


 %% pack all the sampling operators into samplingInfo
 
 samplingInfo.PiAL = PiAL;
 samplingInfo.M1 = M1;
 samplingInfo.M2 = M2;
 samplingInfo.M3 = M3;
 samplingInfo.Pomega1 = Pomega1;
 samplingInfo.Pomega2 = Pomega2;
 samplingInfo.Pomega3 = Pomega3;
 samplingInfo.W1 = W1;
 samplingInfo.W2 = W2;
 samplingInfo.W3 = W3;
 samplingInfo.noise1 = noise1;
 samplingInfo.noise2 = noise2;
 samplingInfo.noise3 = noise3;
 
end

