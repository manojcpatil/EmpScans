function no=Scans_E(samples,ch,n,k,r)
% -------------------------------------------------------------------------
%
%       E = Scans_E(samples,ch,n,k,r)
%
%       This function returns the number of scanning windows of length k
%       containing exactly r occurrences of state ‘ch’ where recounting starts
%       only after the completion of k-tuple of successive trials containing
%       exactly r occurrences of state ‘ch’.
%
%       Example:
%           Input:
%               sample_seq = [1,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0];
%               statistic = Scans_E(sample_seq,1,16,4,2)
%               stat_wait = Scans_WE(sample_seq,1,16,4,2,3)
%
%           Output:
%               statistic = 3
%               stat_wait = 14
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
if nargin<5
    error('Scans_E:TooFewInputs','Requires at least five input arguments.');
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
    for h=1:f
        samples1=(samples==ch(h));
        for i=1:m
            temp=samples1(i,:);
            j=k(h)-1;    const=0;
            while j<n
                j=j+1;
                temp1=sum(temp(j-k(h)+1:j));
                const=const+(temp1==r(h));
                if temp1==r(h)    % exactly r successes in k tuples
                    j=j+k(h)-1;
                end
            end
            no(i,h)=const;
        end
    end
end