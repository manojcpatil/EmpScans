function Table=PrScan_WM(ch,u,k,r,m,P)
% -------------------------------------------------------------------------
%
%       Table=PrScan_WM(ch,u,k,r,m,P)
%
%       This function returns the probabilities of trials required to get u
%       overlapping k-tuples containing at least r successes (occurence of
%       state 'ch') in trials with statespace S={0,1,...,m} and transition
%       probability matrix P, assuming X0=0.
%
%       Example:
%       Input 
%         u=2; k=4; r=2; m=1; ch=1;  P=[0.4,0.6;0.2,0.8];
%         Val_Prob=PrScan_WM(ch,u,k,r,m,P)
% 
%       Output
%         Val_Prob =
% 
%             5.0000    0.8515
%             6.0000    0.0584
%             7.0000    0.0399
%             8.0000    0.0233
%             9.0000    0.0128
%            10.0000    0.0065
%            11.0000    0.0035
%            12.0000    0.0019
%            13.0000    0.0010
%            14.0000    0.0006
%            15.0000    0.0003
%            16.0000    0.0002
%            17.0000    0.0001
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
mat=allposs(k,0:m);L=k-1;
s=sum((m+1).^(1:k));
A=zeros(s);B=zeros(s);

jkl=1;
g=0;
for jkl=1:k-2
    for i=1:((m+1)^jkl)
        for z=0:m
            ej=[sum(mat(1:(m+1)^(jkl+1),1:jkl+1)==repmat([mat(i,1:jkl) z],(m+1)^(jkl+1),1),2)==(jkl+1)]'*[1:(m+1)^(jkl+1)]';
            A(g+i,g+(m+1)^(jkl)+ej)=P(mat(i,jkl)+1,z+1);
        end
    end
    g=g+(m+1)^jkl;
end

jkl=k-1;
for i=1:((m+1)^jkl)
    for z=0:m
        ej=[sum(mat(1:(m+1)^(jkl+1),1:jkl+1)==repmat([mat(i,1:jkl) z],(m+1)^(jkl+1),1),2)==(jkl+1)]'*[1:(m+1)^(jkl+1)]';
        if sum([mat(i,1:jkl) z]==ch)>=r
            B(g+i,g+(m+1)^(jkl)+ej)=P(mat(i,jkl)+1,z+1);
        else
            A(g+i,g+(m+1)^(jkl)+ej)=P(mat(i,jkl)+1,z+1);
        end
    end
end
g=g+i;gL=sum((m+1).^(1:L));

jkl=k;
for i=1:((m+1)^k)
    if sum(mat(i,:)==ch)>=r
        for z=0:m
            if L==0
                ej=[sum(mat(1:(m+1)^(L+1),1:L+1)==repmat([zeros(1,L) z],(m+1)^(L+1),1),2)==L+1]'*[1:(m+1)^(L+1)]';
                A(g+i,gL+ej)=P(mat(i,k)+1,z+1);
            elseif L>=1&&L<=k-2
                ej=[sum(mat(1:(m+1)^(L+1),1:L+1)==repmat([mat(i,k-L+1:k) z],(m+1)^(L+1),1),2)==L+1]'*[1:(m+1)^(L+1)]';
                A(g+i,gL+ej)=P(mat(i,k)+1,z+1);
            elseif L==k-1
                ej=[sum(mat==repmat([mat(i,2:k) z],(m+1)^k,1),2)==k]'*[1:(m+1)^k]';
                if sum([mat(i,2:k) z]==ch)>=r
                    B(g+i,gL+ej)=P(mat(i,k)+1,z+1);
                else
                    A(g+i,gL+ej)=P(mat(i,k)+1,z+1);
                end
            end
        end
    else
        for z=0:m
            ej=[sum(mat==repmat([mat(i,2:k) z],(m+1)^k,1),2)==k]'*[1:(m+1)^k]';
            if sum([mat(i,2:k) z]==ch)>=r
                B(g+i,g+ej)=P(mat(i,k)+1,z+1);
            else
                A(g+i,g+ej)=P(mat(i,k)+1,z+1);
            end
        end
    end
end

ur=k+(u-1)*(k-L);
p=[P(1,:)';zeros(s-m-1,1)];
prob_WML(1:ur-1)=0;j=ur-1;cp=0;
while cp<0.9999
    j=j+1;
    prob_WML(j)=p'*Coef(j-2,u-1,A,B)*B*ones(s,1);
    cp=cp+prob_WML(j);
end

Table=[(ur:j)' prob_WML(ur:j)'];

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
