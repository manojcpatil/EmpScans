function Table=PrScan_WN1(ch,u,k,r,m,P)
% -------------------------------------------------------------------------
%
%       Table=PrScan_WN1(ch,u,k,r,m,P)
%
%       This function returns the probabilities of trials required to 
%       get u scanning windows of length at most k containing r successes
%       (occurence of state 'ch') where recounting starts from scratch each  
%       time r successes encountered in trials with statespace S={0,1,...,m}
%       and transition probability matrix P, assuming X0=0.
%
%       Example:
%       Input 
%         u=2; k=4; r=2; m=1; ch=1;  P=[0.4,0.6;0.2,0.8];
%         Val_Prob=PrScan_WN1(ch,u,k,r,m,P)
%
%       Output
%         Val_Prob =
% 
%             4.0000    0.3072
%             5.0000    0.2611
%             6.0000    0.1805
%             7.0000    0.0972
%             8.0000    0.0604
%             9.0000    0.0378
%            10.0000    0.0235
%            11.0000    0.0137
%            12.0000    0.0079
%            13.0000    0.0046
%            14.0000    0.0026
%            15.0000    0.0015
%            16.0000    0.0009
%            17.0000    0.0005
%            18.0000    0.0003
%            19.0000    0.0002
%            20.0000    0.0001
%            21.0000    0.0000
%            22.0000    0.0000
%            23.0000    0.0000
%            24.0000    0.0000
%            25.0000    0.0000
%            26.0000    0.0000
%            27.0000    0.0000
%            28.0000    0.0000
%            29.0000    0.0000
%            30.0000    0.0000
%            31.0000    0.0000
%            32.0000    0.0000
%
%       Reference:
%       Koutras, M. V. and Alexandrou, V. A. (1995) Runs, scans and urn
%       model distributions: Aunified Markov chain approach, Annals of the
%       Institute of Statistical Mathematics, 47, 743–766.
%
%       See also:
%           PrScan_WN1, PrScan_WN2, PrScan_WM, PrScan_WE and PrScan_WML
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

prob_WN1(1:u*r-1)=0;j=u*r-1;cp=0;
 while cp<0.9999999
     j=j+1;
     prob_WN1(j)=eye(1,nS)*Coef(j-1,u-1,A,B)*B*ones(nS,1);
     cp=cp+prob_WN1(j);
 end

Table=[(u*r:j)' prob_WN1(u*r:j)'];

% x=[u*r:j+2]';
% px=prob_WN1(u*r:j+2)'
% [x,px]

% 
% clc;clear all;
% n=30; k=4; r=2; m=1; ch=1; u=3;
% %P=[0.2,0.2,0.6;0.3,0.5,0.2;0.5,0.3,0.2];
% P=[0.4,0.6;0.2,0.8];
% 
% samples=GenSeq(n,0:m,eye(1,m+1),P,100000);
% y=tabulate(Scans_WN1(samples,1,n,k,r,u))
% 
% % Stored in sampleseq.txt
