function [samplingInfo] =  construct_opt_sampling_distribution(sysInfo,obsInfo,samplingInfo)
% Longxiu Huang, Deanna Needell, and Sui Tang. Robust recovery of bandlimited graph signals via randomized dynamical sampling Information and Inference: A Journal of the IMA, 2024.

% This code is for compute the optimal sampling distributions

A = sysInfo.A;
n = sysInfo.N;
tau = obsInfo.tau;
[V,D] = eigs(A,n);
d = diag(D);
bwidth = sysInfo.bwidth;
d = d(1:bwidth);
Vk = V(:,1:bwidth);

%%%%% need to avoid d=1
Ftau = zeros(bwidth,bwidth);

for i=1:bwidth
    if d(i)==1
        Ftau(i,i)=1/sqrt(tau);
    else   
        Ftau(i,i)=sqrt((1-d(i)^2)./(1-d(i)^(2*tau)));
    end
end

    
% Uk_tau = zeros(tau*n,bwidth);
% for i = 1:tau
%     Uk_tau((i-1)*n+1:i*n,:) = Vk*diag(d.^(i-1))*Ftau;
% end

Vand_k_tau = (bsxfun(@power,d,0:tau-1))'; %%% generate a tau-by-k vandermonde matrix
%%%%%%%%%%%%% define the first optimal distribution p1 %%%%%%%%%%%%%%%%%%%%
p1 = zeros(1,n);
for i = 1:n
    p1(i) = norm(Vand_k_tau*Ftau*diag(Vk(i,:)))^2;
end
samplingInfo.p1 = p1/sum(p1);

%%%%%%%%%%%%%%%% define the second optimal distribution p2 %%%%%%%%%%%%%%%%
p2 = zeros(n,tau);
for i = 1:n
    p2(i,:) = sum((Vand_k_tau*Ftau*diag(Vk(i,:))).^2,2);
end
samplingInfo.p2 = reshape(p2./sum(p2,1),[],1)';

%%%%%%%%%%%%%%%% define the second optimal distribution p3 %%%%%%%%%%%%%%%%
p3 = zeros(1,tau*n);
for t = 1:tau
    temp = Vk*diag(d.^(t-1))*Ftau;
    p3(1,(t-1)*n+1:t*n) = (sum(temp.^2,2));
end
%p3 = (sum(Uk_tau.^2,2))';
samplingInfo.p3 = p3/sum(p3);