function Table=PrScan_ML(ch,n,k,L,r,m,P)
% -------------------------------------------------------------------------
%
%       Table=PrScan_ML(ch,n,k,L,r,m,P)
%
%           This function returns the probabilities of number of
%           L-overlapping k-tuples which contain at least r occurences of
%           state 'ch' in n trials with statespace S={0,1,...,m}
%           and transition probability matrix P, assuming X0=0.
%           
%
%       Example:
%         Input 
%         n=20; k=4; r=2; m=1; ch=1; L=1;
%         P=[0.4,0.6;0.2,0.8];
%         Val_Prob=PrScan_M(ch,n,k,r,m,P)
% 
%         Output
%         Val_Prob =
%              0    0.0000
%         1.0000    0.0001
%         2.0000    0.0008
%         3.0000    0.0067
%         4.0000    0.0449
%         5.0000    0.2335
%         6.0000    0.7139
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
Max_N=floor((n-L)/(k-L));
for i=1:Max_N+1
    prob(i)=p'*C(:,:,n-1,i)*ones(s,1);
end
Table=[[0:Max_N]' prob'];