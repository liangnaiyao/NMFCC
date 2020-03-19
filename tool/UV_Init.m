function [U,V]=UV_Init(X,U_,V_,k)
[mFea,nSmp]=size(X);
if isempty(U_)
    U = abs(rand(mFea,k));
    norms = sqrt(sum(U.^2,1));
    norms = max(norms,1e-10);
    U = U./repmat(norms,mFea,1);
    if isempty(V_)
        V = abs(rand(nSmp,k));
        V = V/sum(sum(V));
    else
        V = V_;
    end
else
    U = U_;
    if isempty(V_)
        V = abs(rand(nSmp,k));
        V = V/sum(sum(V));
    else
        V = V_;
    end
end