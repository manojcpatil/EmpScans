function no=Scans_ML(samples,ch,n,k,L,r)
% -------------------------------------------------------------------------
%
%       ML = Scans_ML(samples,ch,n,k,L,r)
%           This function returns the number of scanning windows of length k
%           containing at least r successes (occurence of 'ch') which may
%           have overlapping part of length L with the previous scanning
%           window that has been counted.
%
%       Example
%           Input:
%               sample_seq = [1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0];
%               statistic = Scans_ML(sample_seq,1,16,4,1,2)
%               stat_wait = Scans_WML(sample_seq,1,16,4,1,2,3)
%
%           Output:
%               statistic = 4
%               stat_wait = 10
%
%       Reference:
%       Aki, S. and Hirano, K. (2000) Number of success runs of specified
%       length until certain stopping time rules and generalized binomial
%       distributions of order k, Annals of the Institute of Statistical
%       Mathematics, 52, 767–777.
%
%       See also:
%           Scans_N1, Scans_N2, Scans_M, Scans_E and Scans_ML
%
%       EmpScans Toolbox: Toolbox for Scan statistics with MATLAB
%       Authors:   Patil Manoj C. and Shinde Ramakrishna L., 2019.
% -------------------------------------------------------------------------
if nargin<6
    error('Scans_ML:TooFewInputs','Requires at least six input arguments.');
else
    n_ch=length(ch);n_k=length(k);n_r=length(r);,n_L=length(L);
    n_max=max([n_ch,n_k,n_r,n_L]);
    if n_max>1
        if n_ch==1
            ch=ch*ones(n_max,1);
        end
        if n_k==1
            k=k*ones(n_max,1);
        end
        if n_r==1
            r=r*ones(n_max,1);
        end
        if n_L==1
            L=L*ones(n_max,1);
        end
    end
    [m,n]=size(samples);
    f=length(r);
    for h=1:f
        samples1=(samples==ch(h));
        for i=1:m
            temp=samples1(i,:);
            j=k(h)-1;    const=0;
            while j<n
                j=j+1;
                temp1=sum(temp(j-k(h)+1:j));
                const=const+(temp1>=r(h));
                j=j+(k(h)-1-L(h))*(temp1>=r(h));
            end
            no(i,h)=const;
        end
    end
end