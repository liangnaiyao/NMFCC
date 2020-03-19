function [U_final, V_final, objhistory_final] = NMFCC(X, k, options, U_, V_, U, V, A, B, C)

% �������ͣ� tryNo������ǰ���Դ���
%            nRepeat�����ܹ��ĳ��Դ���
%            selectInit�������ž���õ��ķ�����1Ϊ���ݳ��Դ���nRepeatĿ�꺯��ֵ��С����ѡ�񣬴ﵽ��С���������������Ϊ0�Ļ��ο�Դ�ļ��������޸ġ�
%           AΪ����X�� BΪ����U�� CΪ����V��

differror = options.error;
maxIter = options.maxIter;
nRepeat = options.nRepeat;
minIterOrig = options.minIter;
minIter = minIterOrig-1;
objhistory_final=0;
gamma=options.gamma;

% [U,V] = NormalizeUV(U, V, 0, 1);   % ��һ��U��V

selectInit = 1;   % selectInit һ��ҪΪ1������ѡ�����Ŀ�꺯��ֵ�Ĺ���

tryNo = 0;
while tryNo < nRepeat   
    tryNo = tryNo+1;  % ���Դ���
    nIter = 0;
    maxErr = 1;
    while(maxErr > differror)  %   =====  ��ʼһ�γ���  ===== 
        nIter = nIter + 1;   
  
        [U, V]=update(X, U, V, options, A, B, C);  % ���µ���
       
% ======== ���ﵽ��С�����������������ͱ��  ==============
        if nIter > minIter   %  ������minIter��ʼ�ж��������������û���ж�������������
            if selectInit
                objhistory = CalculateObj(X, U, V);
                maxErr = 0;
            else    %  nRepeat �����һ��ִ��
                    maxErr = 1;
                    if nIter >= maxIter
                        maxErr = 0;
                        objhistory = 0;  % ���0��������һ�ν����ֵ��U_final
                    end
            end
        end
    end       %  ======  ����һ�γ���  =======
    
% ===========  �ԱȲ�ͬ���Ե�Ŀ�꺯��ֵ������Ŀ�꺯����С�ľ��󵽺�׺Ϊ_final�ı�����  ==============
    if tryNo == 1       % ��һ�ε������뱻�Ƚ϶���
        U_final = U;
        V_final = V;
        objhistory_final = objhistory;
    else
       if objhistory(end) < objhistory_final(end) % �����ǵ�һ�ε��������ｲ�Ҽ���nRepeat��Ŀ�꺯����С���Ǵε����ķֽ���
           U_final = U;
           V_final = V;
           objhistory_final = objhistory;
       end
    end

% ================ �ж��ǲ��ǵ������һ�γ��� ========================    
    if selectInit  %��selectInit��Ϊ0��ʱ��϶���ִ��һ�Σ����������¿�ʼ����������һ�αȽ�Ŀ�꺯���Ķ���
        if tryNo < nRepeat %nRepeatΪ2�ǣ���һ�ε������ʱ��Ϊ1���ڶ��ε�����Ͳ���ִ���ˣ�ִ��159��else
            %re-start
          [U,V]=UV_Init(X,U_,V_,k); % ����U,V��ֵ�����������г���     
        else
            tryNo = tryNo - 1;  % ����tryNo-1��Ϊ����ִ��һ��
            minIter = 0;        % ����ִ�У�ֻ������һ��
            selectInit = 0;  % selectInitΪ0��ʱ��ڶ��ξͲ���ִ������
            U = U_final;
            V = V_final;     % �����ڶ���V_final�ḳֵ��U��V ȥ������һ�γ���           
        end
    end
end

function [U, V]=update(X, U, V, options, A, B, C)
%           AΪ����X�� BΪ����U�� CΪ����V��
       alpha=options.alpha;
       beta=options.beta;  
       sigma=options.sigma;
       gamma=options.gamma;
       mu=options.mu;
       if gamma>0
         L=options.L;
         DD=options.D;
         S=options.S;
       end
       nv_i=options.nv_i;   % ��ǰ�ӽ�
        % ===================== update V ========================
        ts=0;  tx=0;
        if gamma > 0
            ts=S{nv_i}*V; 
            tx= DD{nv_i}*V;
        end       
        
        XU = X'*U+mu*V+gamma*ts;  % mnk or pk (p<<mn)
        UU = U'*U;  % mk^2
        
       K_di = 0; % ����     
       if alpha > 0
         for k=1:length(A)
            if (k==nv_i) 
                continue;
            end
            K_di =  K_di + C{k};  
         end
        end    
        VUU = V*UU+0.5*alpha*K_di+mu*V*(V'*V)+gamma*tx; % nk^2
        
        V = V.*(XU./max(VUU,eps));
        V=double(V);
        
        % ===================== update U ========================
        XV = X*V+sigma*U;   % mnk or pk (p<<mn)
        VV = V'*V;  % nk^2
        UVV = U*VV+sigma*repmat(sum(U,2),1,size(U,2)); % mk^2f
        
        U = U.*(XV./max(UVV,eps)); % 3mk
        U=double(U);

     

        