function Table=PrScan_N1(ch,n,k,r,m,P)
% -------------------------------------------------------------------------
%
%       Table=PrScan_N1(ch,n,k,r,m,P)
%
%           This function returns the probabilities of number of windows
%           of length at most k containing r successes (occurence of 'ch')
%           where recounting starts from scratch each time r successes
%           encountered in n trials with statespace S={0,1,...,m} and
%           transition probability matrix P, assuming X0=0.
%
%       Example:
%       Input
%       n=20; k=4; r=2; m=1; ch=1;
%       P=[0.4,0.6;0.2,0.8];
%       Val_Prob=PrScan_N1(ch,n,k,r,m,P)
%
%       Output
%       Val_Prob =
%            0    0.0000
%       1.0000    0.0001
%       2.0000    0.0010
%       3.0000    0.0064
%       4.0000    0.0276
%       5.0000    0.0857
%       6.0000    0.1900
%       7.0000    0.2900
%       8.0000    0.2748
%       9.0000    0.1158
%      10.0000    0.0086
%
%       Reference:
%       Koutras, M. V. and Alexandrou, V. A. (1995) Runs, scans and urn
%       model distributions: Aunified Markov chain approach, Annals of the
%       Institute of Statistical Mathematics, 47, 743–766.
%
%       See also:
%           PrScan_N1, PrScan_N2, PrScan_M, PrScan_E and PrScan_ML
%
%       EmpScans Toolbox: Toolbox for Scan statistics with MATLAB
%       Authors:   Patil Manoj C. and Shinde Ramakrishna L., 2019.
% -------------------------------------------------------------------------

mat=allposs(k,0:m);
mat=mat(:,k:-1:1);
m1=sum(mat==ch,2);        m2=sum(mat(:,1:k-1)==ch,2);
S1=mat(m1<=r-1,:);    S2=mat((m2==r-1).*(mat(:,k)==ch)==1,:);
S=[S1;S2];

nS=size(S,1);       nS1=size(S1,1);
A=zeros(nS);B=zeros(nS);
for i=1:nS1;
    for j=1:nS
        for z=0:m
            if min([S(i,[2:k]) z]==S(j,:))==1
                A(i,j)=P(S(i,k)+1,z+1);
            end
        end
    end
end

A(nS1+1:nS,1:m+1)=repmat(P(ch+1,:),nS-nS1,1);

B=A;
B(:,1:nS1)=0;
A=A-B;

C(:,:,1,1)=A;
C(:,:,1,2)=B;

for j=2:n
    C(:,:,j,1)=  C(:,:,j-1,1)*A;
    C(:,:,1,j)=  zeros(nS);
end
for j=2:n
    for i=2:n
        C(:,:,j,i)=  C(:,:,j-1,i)*A + C(:,:,j-1,i-1)*B;
    end
end

Max_N=floor(n/r);
for i=1:Max_N+1
    prob(i)=eye(1,nS)*C(:,:,n,i)*ones(nS,1);
end
Table=[[0:Max_N]' prob'];