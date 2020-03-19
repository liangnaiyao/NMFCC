function [U_final, V_final, objhistory_final] = NMFCC(X, k, options, U_, V_, U, V, A, B, C)

% 变量解释： tryNo——当前尝试次数
%            nRepeat——总共的尝试次数
%            selectInit——最优矩阵得到的方法，1为根据尝试次数nRepeat目标函数值最小进行选择，达到最小迭代比跳出。如果为0的话参考源文件，进行修改。
%           A为所有X； B为所有U； C为所有V；

differror = options.error;
maxIter = options.maxIter;
nRepeat = options.nRepeat;
minIterOrig = options.minIter;
minIter = minIterOrig-1;
objhistory_final=0;
gamma=options.gamma;

% [U,V] = NormalizeUV(U, V, 0, 1);   % 归一化U、V

selectInit = 1;   % selectInit 一定要为1，才有选择最低目标函数值的功能

tryNo = 0;
while tryNo < nRepeat   
    tryNo = tryNo+1;  % 尝试次数
    nIter = 0;
    maxErr = 1;
    while(maxErr > differror)  %   =====  开始一次尝试  ===== 
        nIter = nIter + 1;   
  
        [U, V]=update(X, U, V, options, A, B, C);  % 更新迭代
       
% ======== 当达到最小迭代次数，给出解释标号  ==============
        if nIter > minIter   %  迭代到minIter开始判断收敛情况，这里没有判断收敛，输出结果
            if selectInit
                objhistory = CalculateObj(X, U, V);
                maxErr = 0;
            else    %  nRepeat 的最后一次执行
                    maxErr = 1;
                    if nIter >= maxIter
                        maxErr = 0;
                        objhistory = 0;  % 这个0负责把最后一次结果赋值回U_final
                    end
            end
        end
    end       %  ======  结束一次尝试  =======
    
% ===========  对比不同尝试的目标函数值，保存目标函数较小的矩阵到后缀为_final的变量中  ==============
    if tryNo == 1       % 第一次迭代存入被比较对象
        U_final = U;
        V_final = V;
        objhistory_final = objhistory;
    else
       if objhistory(end) < objhistory_final(end) % 当不是第一次迭代，这里讲找几次nRepeat中目标函数最小的那次迭代的分解结果
           U_final = U;
           V_final = V;
           objhistory_final = objhistory;
       end
    end

% ================ 判断是不是到了最后一次尝试 ========================    
    if selectInit  %当selectInit不为0的时候肯定会执行一次，这里是重新开始迭代好有下一次比较目标函数的对象
        if tryNo < nRepeat %nRepeat为2是，第一次到这里的时候为1，第二次到这里就不会执行了，执行159的else
            %re-start
          [U,V]=UV_Init(X,U_,V_,k); % 打乱U,V的值，以重新运行程序     
        else
            tryNo = tryNo - 1;  % 这里tryNo-1是为了再执行一次
            minIter = 0;        % 重新执行，只迭代了一次
            selectInit = 0;  % selectInit为0的时候第二次就不会执行这里
            U = U_final;
            V = V_final;     % 倒数第二次V_final会赋值回U、V 去再运行一次程序           
        end
    end
end

function [U, V]=update(X, U, V, options, A, B, C)
%           A为所有X； B为所有U； C为所有V；
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
       nv_i=options.nv_i;   % 当前视角
        % ===================== update V ========================
        ts=0;  tx=0;
        if gamma > 0
            ts=S{nv_i}*V; 
            tx= DD{nv_i}*V;
        end       
        
        XU = X'*U+mu*V+gamma*ts;  % mnk or pk (p<<mn)
        UU = U'*U;  % mk^2
        
       K_di = 0; % 清零     
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

     

        