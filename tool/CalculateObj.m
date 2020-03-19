function [obj, dV] = CalculateObj(X, U, V, deltaVU, dVordU)
    if ~exist('deltaVU','var')
        deltaVU = 0;
    end
    if ~exist('dVordU','var')
        dVordU = 1;
    end
    dV = [];
    maxM = 62500000;
    [mFea, nSmp] = size(X);
    mn = numel(X);
    nBlock = floor(mn*3/maxM);

    if mn < maxM
        dX = U*V'-X;
        obj_NMF = sum(sum(dX.^2));
        if deltaVU
            if dVordU
                dV = dX'*U;
            else
                dV = dX*V;
            end
        end
    else
        obj_NMF = 0;
        if deltaVU
            if dVordU
                dV = zeros(size(V));
            else
                dV = zeros(size(U));
            end
        end
        for i = 1:ceil(nSmp/nBlock)
            if i == ceil(nSmp/nBlock)
                smpIdx = (i-1)*nBlock+1:nSmp;
            else
                smpIdx = (i-1)*nBlock+1:i*nBlock;
            end
            dX = U*V(smpIdx,:)'-X(:,smpIdx);
            obj_NMF = obj_NMF + sum(sum(dX.^2));
            if deltaVU
                if dVordU
                    dV(smpIdx,:) = dX'*U;
                else
                    dV = dU+dX*V(smpIdx,:);
                end
            end
        end
        if deltaVU
            if dVordU
                dV = dV ;
            end
        end
    end
   %obj_Lap = alpha*sum(sum((L*V).*V));
   
    obj = obj_NMF;