function Table=PrScan_E(ch,n,k,r,m,P)
% -------------------------------------------------------------------------
%
%       Table=PrScan_E(ch,n,k,r,m,P)
%
%       This function returns the probabilities of number of windows of length k
%       containing exactly r occurrences of state ‘ch’ where recounting starts
%       only after the completion of k-tuple of successive trials containing
%       exactly r occurrences of state ‘ch’ in n trials with statespace 
%       S={0,1,...,m} and transition probability matrix P, assuming X0=0.
%
%       Example:
%         Input 
%         n=20; k=4; r=2; m=1; ch=1;
%         P=[0.4,0.6;0.2,0.8];
%         Val_Prob=PrScan_N2(ch,n,k,r,m,P)
% 
%         Output
%         Val_Prob =
%                  0    0.0000
%             1.0000    0.0002
%             2.0000    0.0032
%             3.0000    0.0419
%             4.0000    0.3412
%             5.0000    0.6135
%
%       Reference:
%       Mood, A. M. (1940) The Distribution Theory of Runs, Annals of
%       Mathematical Statistics, 11, 367–392.
%
%       See also:
%           Scans_N1, Scans_N2, Scans_M, Scans_E and Scans_ML
%
%       EmpScans Toolbox: Toolbox for Scan statistics with MATLAB
%       Authors:   Patil Manoj C. and Shinde Ramakrishna L., 2019.
% -------------------------------------------------------------------------
% 

mat=allposs(k,0:m);L=0;
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
        if sum([mat(i,1:jkl) z]==ch)==r
            B(g+i,g+(m+1)^(jkl)+ej)=P(mat(i,jkl)+1,z+1);
        else
            A(g+i,g+(m+1)^(jkl)+ej)=P(mat(i,jkl)+1,z+1);
        end
    end
end
g=g+i;gL=sum((m+1).^(1:L));

jkl=k;
for i=1:((m+1)^k)
    if sum(mat(i,:)==ch)==r
        for z=0:m
            if L==0
                ej=[sum(mat(1:(m+1)^(L+1),1:L+1)==repmat([zeros(1,L) z],(m+1)^(L+1),1),2)==L+1]'*[1:(m+1)^(L+1)]';
                A(g+i,gL+ej)=P(mat(i,k)+1,z+1);
            elseif L>=1&&L<=k-2
                ej=[sum(mat(1:(m+1)^(L+1),1:L+1)==repmat([mat(i,k-L+1:k) z],(m+1)^(L+1),1),2)==L+1]'*[1:(m+1)^(L+1)]';
                A(g+i,gL+ej)=P(mat(i,k)+1,z+1);
            elseif L==k-1
                ej=[sum(mat==repmat([mat(i,2:k) z],(m+1)^k,1),2)==k]'*[1:(m+1)^k]';
                if sum([mat(i,2:k) z]==ch)==r
                    B(g+i,gL+ej)=P(mat(i,k)+1,z+1);
                else
                    A(g+i,gL+ej)=P(mat(i,k)+1,z+1);
                end
            end
        end
    else
        for z=0:m
            ej=[sum(mat==repmat([mat(i,2:k) z],(m+1)^k,1),2)==k]'*[1:(m+1)^k]';
            if sum([mat(i,2:k) z]==ch)==r
                B(g+i,g+ej)=P(mat(i,k)+1,z+1);
            else
                A(g+i,g+ej)=P(mat(i,k)+1,z+1);
            end
        end
    end
end

C(:,:,1,1)=A;
C(:,:,1,2)=B;
C(:,:,2,1)=C(:,:,1,1)*A;

for j=3:n+1
    C(:,:,j,1)=  C(:,:,j-1,1)*A;
    C(:,:,1,j)=  zeros(s);
end

for j=2:n
    for i=2:n+1
        C(:,:,j,i)=  C(:,:,j-1,i)*A + C(:,:,j-1,i-1)*B;
    end
end
p=[P(1,:)';zeros(s-m-1,1)];
Max_E=floor(n/k);
for i=1:Max_E+1
    prob(i)=p'*C(:,:,n-1,i)*ones(s,1);
end
Table=[[0:Max_E]' prob'];
