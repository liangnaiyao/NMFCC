function [V_x,V] = main_NMFCC(X, label,sigma,alpha,gamma,mu, options)
% sum_p||X^p-U^p(V^p)'||+2*alpha*sum_p sum_w(tr((V_p)'V_w))
%  +sigma*sum_p (1'(U^p)'U^p1-tr((U^p)'U^p))+gamma*sum_p tr((V^p)'L^pV^p)
%  + 0.5*mu*sum_p||(V^p)'V^p-I||
% s.t. X,U,V>=0.

optionDefault.maxIter = 200;
optionDefault.error = 1e-6;
optionDefault.nRepeat = 2; % ����Ϊ1
optionDefault.minIter = 100;
if ~exist('option','var')
   options=optionDefault;
else
    options=mergeOption(options,optionDefault);  % ��������ṹ��
end 
options.beta=alpha;
options.sigma=sigma; % ��
options.alpha=alpha; % �ӽǼ������
options.gamma=gamma;  % ͼ
options.mu=mu; % �ӽ��ڵĶ�����


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

U_ = [];  %Ҫ����Ϊ�յģ���ΪNMF�ĳ�ʼ�������õ�
V_ = [];  %Ҫ����Ϊ�յ�

A = 0;
log = 0;
    for i = 1:viewNum 
        [U{i},V{i}]=UV_Init(X{i},U_,V_,K);
    end
    
    for i = 1 : viewNum
        options.nv_i=i;   % ��ǰ���ӽ�
        [U{i}, V{i}] = NMFCC(X{i}, K, options, U_, V_, U{i}, V{i}, X, U, V);        
        A=A+V{i} ;
    end
    V_x=A/viewNum;
