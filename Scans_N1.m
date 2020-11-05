function no=Scans_N1(samples,ch,n,k,r)
% -------------------------------------------------------------------------
%
%       N1 = Scans_N1(samples,ch,n,k,r)
%
%           This function returns the Number of windows of length
%           at most k containing r successes where recounting starts from
%           scratch each time r successes encountered.
%
%       Example:
%       Input:
%           sample_seq = [1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0];
%           statistic  = Scans_N1(sample_seq,1,16,4,2)
%           stat_wait  = Scans_WN1(sample_seq,1,16,4,2,3)
%
%       Output
%           statistic  = 4
%           stat_wait  = 11
%
%       Reference:
%       Koutras, M. V. and Alexandrou, V. A. (1995) Runs, scans and urn
%       model distributions: Aunified Markov chain approach, Annals of the
%       Institute of Statistical Mathematics, 47, 743–766.
%
%       See also:
%           Scans_N1, Scans_N2, Scans_M, Scans_E and Scans_ML
%
%       EmpScans Toolbox: Toolbox for Scan statistics with MATLAB
%       Authors:   Patil Manoj C. and Shinde Ramakrishna L., 2019.
% -------------------------------------------------------------------------
if nargin<5
    error('Scans_N1:TooFewInputs','Requires at least five input arguments.');
else
    n_ch=length(ch);n_k=length(k);n_r=length(r);
    n_max=max([n_ch,n_k,n_r]);
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
    end

    [m,n]=size(samples);
    f=length(r);
    no=zeros(m,f);
    for h=1:f
        if r(h)==1
            no(:,h)=sum((samples==ch(h)),2);
        else
            samples1=(samples==ch(h));
            for i=1:m
                temp=[zeros(k(h)-r(h),1);samples1(i,:)'];
                const=0;j=k(h);
                while j<=n+k(h)-r(h)
                    temp1=sum(temp((j-k(h)+1:j)));
                    const=const+(temp1==r(h));
                    if temp1==r(h)
                        temp(j-k(h)+1:j)=0;    j=j+k(h)-r(h)-1;
                        % As recounting starts only after completion of r success in
                        % atmost k tuples
                    end
                    j=j+1;
                end
                no(i,h)=const;clear temp temp1 const;
            end
        end
    end
end