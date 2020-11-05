function y=GenSeq(n,states,alpha,TPM,samples)
% GenSeq(n,states,alpha,TPM,samples)
calpha=[0 cumsum(alpha,2)];
cTPM=[zeros(length(TPM),1) cumsum(TPM,2)];

for i=1:samples
    fprintf('%d :',i)
    r=rand(1);
    x=sum(calpha<r);
    for j=2:n+1
        r=rand(1);
        y(i,j-1)=sum(cTPM(x,:)<r);
        x=y(i,j-1);
        fprintf('.');
    end
    clear x;
    clc;
end
y=states(y);
