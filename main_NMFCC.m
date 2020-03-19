function [V_x,V] = main_NMFCC(X, label,sigma,alpha,gamma,mu, options)
% sum_p||X^p-U^p(V^p)'||+2*alpha*sum_p sum_w(tr((V_p)'V_w))
%  +sigma*sum_p (1'(U^p)'U^p1-tr((U^p)'U^p))+gamma*sum_p tr((V^p)'L^pV^p)
%  + 0.5*mu*sum_p||(V^p)'V^p-I||
% s.t. X,U,V>=0.

optionDefault.maxIter = 200;
optionDefault.error = 1e-6;
optionDefault.nRepeat = 2; % 不能为1
optionDefault.minIter = 100;
if ~exist('option','var')
   options=optionDefault;
else
    options=mergeOption(options,optionDefault);  % 混合两个结构体
end 
options.beta=alpha;
options.sigma=sigma; % 基
options.alpha=alpha; % 视角间多样性
options.gamma=gamma;  % 图
options.mu=mu; % 视角内的多样性


K=length(unique(label));
viewNum = length(X);
% Construct graph Laplacian
nSmp=size(X{1},2);
if options.gamma>0
for v = 1:length(X)
    options.WeightMode='Binary';
    W{v}=constructW_cai(X{v}',options); 
    DCol = full(sum(W{v},2));
    D{v} = spdiags(DCol,0,nSmp,nSmp);
    L{v} = D{v} - W{v};
    if isfield(options,'NormW') && options.NormW
        D_mhalf = spdiags(DCol.^-.5,0,nSmp,nSmp) ;
        L{v} = D_mhalf*L{v}*D_mhalf;
    end
end
 options.L=L;   
  options.D=D;
 options.S=W;
end

U_ = [];  %要保持为空的，因为NMF的初始化矩阵用到
V_ = [];  %要保持为空的

A = 0;
log = 0;
    for i = 1:viewNum 
        [U{i},V{i}]=UV_Init(X{i},U_,V_,K);
    end
    
    for i = 1 : viewNum
        options.nv_i=i;   % 当前的视角
        [U{i}, V{i}] = NMFCC(X{i}, K, options, U_, V_, U{i}, V{i}, X, U, V);        
        A=A+V{i} ;
    end
    V_x=A/viewNum;
