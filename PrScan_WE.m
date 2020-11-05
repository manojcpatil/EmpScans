function Table=PrScan_WE(ch,u,k,r,m,P)
% -------------------------------------------------------------------------
%
%       Table=PrScan_WE(ch,u,k,r,m,P)
%
%       This function returns the probabilities of number trials required
%       to get u scanning windows of length k containing exactly r
%       occurrences of state ‘ch’ where recounting starts only after the
%       completion of k-tuple of successive trials containing exactly r
%       occurrences of state ‘ch’ in n trials with statespace S={0,1,...,m}
%       and transition probability matrix P, assuming X0=0.
%       
%       Example:
%       Input: 
%       u=2; k=4; r=2; m=1; ch=1;
%       P=[0.4,0.6;0.2,0.8];
%       Val_Prob=PrScan_WE(ch,u,k,r,m,P)
%       
%       Output:
%       Val_Prob =
% 
%         8.0000    0.0449
%         9.0000    0.0449
%        10.0000    0.0486
%        11.0000    0.0488
%        12.0000    0.0481
%        13.0000    0.0461
%        14.0000    0.0441
%        15.0000    0.0421
%        16.0000    0.0401
%        17.0000    0.0381
%        18.0000    0.0362
%        19.0000    0.0344
%        20.0000    0.0325
%        21.0000    0.0307
%        22.0000    0.0290
%        23.0000    0.0273
%        24.0000    0.0256
%        25.0000    0.0241
%        26.0000    0.0226
%        27.0000    0.0211
%        28.0000    0.0198
%        29.0000    0.0185
%        30.0000    0.0172
%        31.0000    0.0161
%        32.0000    0.0150
%        33.0000    0.0139
%        34.0000    0.0130
%        35.0000    0.0120
%        36.0000    0.0112
%        37.0000    0.0104
%        38.0000    0.0096
%        39.0000    0.0089
%        40.0000    0.0083
%        41.0000    0.0076
%        42.0000    0.0071
%        43.0000    0.0065
%        44.0000    0.0060
%        45.0000    0.0056
%        46.0000    0.0051
%        47.0000    0.0048
%        48.0000    0.0044
%        49.0000    0.0040
%        50.0000    0.0037
%        51.0000    0.0034
%        52.0000    0.0032
%        53.0000    0.0029
%        54.0000    0.0027
%        55.0000    0.0025
%        56.0000    0.0023
%        57.0000    0.0021
%        58.0000    0.0019
%        59.0000    0.0018
%        60.0000    0.0016
%        61.0000    0.0015
%        62.0000    0.0014
%        63.0000    0.0013
%        64.0000    0.0012
%        65.0000    0.0011
%        66.0000    0.0010
%        67.0000    0.0009
%
%       Reference:
%       Mood, A. M. (1940) The Distribution Theory of Runs, Annals of
%       Mathematical Statistics, 11, 367–392.
%
%       See also:
%           PrScan_WN1, PrScan_WN2, PrScan_WM, PrScan_WE and PrScan_WML
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

ur=u*k;
p=[P(1,:)';zeros(s-m-1,1)];
prob_WE(1:ur-1)=0;j=ur-1;cp=0;
while cp<0.99
    j=j+1;
    prob_WE(j)=p'*Coef(j-2,u-1,A,B)*B*ones(s,1);
    cp=cp+prob_WE(j);
end

Table=[(ur:j)' prob_WE(ur:j)']