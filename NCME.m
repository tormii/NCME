function [IFmset,IA,smset] = NCME(s,fs,eIF,lambda,beta,tol)
% Authors: Xiaotong Tu and Fucai Li
% email:xttu@xmu.edu.cn, tormiier@gmail.com;
% https://www.researchgate.net/profile/Xiaotong_Tu2
% https://tormii.github.io/
%%%%%%%%%%%%%%%%%%%%%%%  input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s: measured signal,a row vector
% fs: sampling frequency
% eIF: initial instantaneous frequency (IF) time series for all the signal modes; each row of eIF corresponds to the IF of each mode
% lambda: allows for a tuning of the relative weighting between the cost function and the desired degree of sparseness of the solution
% beta: penalty parameter controling the smooth degree of the IF increment during iterations;the smaller the beta is, the more smooth the IF increment would be
% tol: tolerance of convergence criterion; 
%%%%%%%%%%%%%%%%%%%%%%% output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IFmset: the collection of the obtained IF time series of all the signal modes at each iteration
% smset: the collection of the obtained signal modes at each iteration
% IA: the finally estimated instantaneous amplitudes of the obtained signal modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% When using this code, please cite our paper:
% Xiaotong Tu, Johan Sward, Andreas Jakobsson, and Fucai Li, Estimating Nonlinear Chirp Modes Exploiting Sparsity, Signal Processing
% Author: Xiaotong Tu£¨Jan.,2021£©

iternum = 50; %the maximum allowable iterations

[K,N] = size(eIF);%K is the number of the components£¬N is thenumber of the samples
t = (0:N-1)/fs;%time
e = ones(N,1);
e2 = -2*e;
e2(1) = -1;e2(end) = -1;
oper = spdiags([e e2 e], -1:1, N, N);%the modified second-order difference matrix
opedoub = oper'*oper;%
sinm = zeros(K,N);cosm = zeros(K,N);%
IFsetiter = zeros(K,N,iternum+1); IFsetiter(:,:,1) = eIF; %the collection of the obtained IF time series of all the signal modes at each iteration
ssetiter = zeros(K,N,iternum+1); %the collection of the obtained signal modes at each iteration
sDif = tol + 1;%
iter = 1;% iteration counter

for i = 1:K
    sinm(i,:) = sin(2*pi*(cumtrapz(t,eIF(i,:))));
    cosm(i,:) = cos(2*pi*(cumtrapz(t,eIF(i,:))));
    Bm = spdiags(sinm(i,:)', 0, N, N);
    Am = spdiags(cosm(i,:)', 0, N, N);
end
%% CVX LASSO
while ( sDif > tol*1e-2 &&  iter <= iternum ) % 
    betathr = 10^(iter/36-10);%gradually increase the parameter beta during the iterations
    if betathr>beta
        betathr = beta; 
    end
    %ADMM solver
    A=[Am Bm]; D=blkdiag(oper,oper);b=s';rho=1;
    [CVX_u_v, ~] = gene_lasso(A,D, b, lambda, rho, 1);  
    [~,n]=size(s);
    CVX_xm=CVX_u_v(1:n)';
    CVX_ym=CVX_u_v(n+1:end)';
    %%%%%%%%%%%%%  update the IFs  %%%%%%%%%%%%%%%%%%%%%%%%  
       i=1;% for i = 1:K
       ybar = Differ(CVX_ym(i,:),1/fs); xbar = Differ(CVX_xm(i,:),1/fs);%compute the derivative of the functions
       deltaIF = (CVX_xm(i,:).*ybar - CVX_ym(i,:).*xbar)./(CVX_xm(i,:).^2 + CVX_ym(i,:).^2)/2/pi;% obtain the frequency increment by arctangent demodulation
       deltaIF = (2/betathr*opedoub + speye(N))\deltaIF';% smooth the frequency increment by low pass filtering

       eIF(i,:) = eIF(i,:) - 0.5*deltaIF';% update the IF
       IFsetiter(:,:,iter+1) = eIF;
       %%%%%%%%%%%%%  update cos and sin functions  %%%%%%%%%%%%%%%%%%%%%%%% 
       sinm(i,:) = sin(2*pi*(cumtrapz(t,eIF(i,:))));
       cosm(i,:) = cos(2*pi*(cumtrapz(t,eIF(i,:))));
       %%%%%%%%%%%%% update sums %%%%%%%%%%%%%%%%%       
%        sum_x = CVX_xm(i,:).*cosm(i,:); %
%        sum_y = CVX_ym(i,:).*sinm(i,:);%
       ssetiter(i,:,iter+1) = CVX_xm(i,:).*cosm(i,:) + CVX_ym(i,:).*sinm(i,:);%
       %%%%%%%%%%%%% update the two matrices A and B %%%%%%%%%%%%%%%%%%%%%%%%%%%         
       Bm = spdiags(sinm(i,:)', 0, N, N);Bdoubm = spdiags((sinm(i,:).^2)', 0, N, N);
       Am = spdiags(cosm(i,:)', 0, N, N);Adoubm = spdiags((cosm(i,:).^2)', 0, N, N);
       %%%%%%%%%%%%%  compute the convergence index %%%%%%%%%%%%%%%%%%  
       sDif = 0;
       sDif = sDif + (norm(ssetiter(i,:,iter+1) - ssetiter(i,:,iter))/norm(ssetiter(i,:,iter))).^2
%      lambda=abs(lambda* ((ssetiter(i,:,iter+1)*ssetiter(i,:,iter+1)')/(ssetiter(i,:,iter+1)*s')))
        weight=abs( ((ssetiter(i,:,iter+1)*ssetiter(i,:,iter+1)')/(ssetiter(i,:,iter+1)*s')))
%       lambda=(0.5+iter/20)
%       if weight <1
%           lambda=lambda/weight
%       end
       iter = iter + 1 
       
   
end
 IFmset= IFsetiter(:,:,1:iter);
 smset = ssetiter(:,:,1:iter);
 IA = sqrt(CVX_xm.^2 + CVX_ym.^2);



end