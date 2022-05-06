function [eiglist,func] = getEigsCyl(a,b,nmax,typebounds,dx,cutoff)
% get eigenvalues for cylindrically symmetric diffusion problem
% in annular region between radius a and b
% nmax is the number of eigenvalues to calculate
% typebounds = for inner, outer boundary, 0 for absorbing; 1 for reflecting
% dx is the step in initial x values for finding zeros
% default is 1

if (~exist('dx','var'))
    dx = 5/b;
end
if (~exist('cutoff','var'))
    cutoff = 2*pi*5;
end

if (typebounds(1)==0 & typebounds(2)==1)
    % inner absorbing, outer reflecting
    func =@(x) real(besselj(1,x*b).*bessely(0,x*a) - besselj(0,x*a).*bessely(1,x*b));
elseif(typebounds(1)==1 & typebounds(2)==0)
    % inner reflecting, outer absorbing
    func = @(x) real(besselj(1,x*a).*bessely(0,x*b) - besselj(0,x*b).*bessely(1,x*a));
elseif (typebounds(1)==0 & typebounds(2)==0)
    % both bounds absorbing
    func = @(x) real(besselj(0,x*b).*bessely(0,x*a) - besselj(0,x*a).*bessely(0,x*b));    
end

eiglist = [];
ct = 0;

x0 = [dx/2,3*dx/2];
while length(eiglist)<nmax    
    if (func(x0(1))*func(x0(2)) > 0)        
        % same sign, just use this starting point
        val =  fzero(func,mean(x0));
    else % different sign, search in this interval
        val = fzero(func, x0);
    end
    %[x0,val]
    if (isempty(eiglist))
        ct=ct+1;
        eiglist(ct) = val;
        if (val>cutoff)
            x0 = val+pi/(b-a);
            x0 = [x0-dx,x0+dx];
        end
    elseif (val > eiglist(end)+1e-8/b)
        ct =ct+1;
        eiglist(ct) = val;
        if (val>cutoff)
            x0 = val+pi/(b-a);
            x0 = [x0-dx,x0+dx]; % try an interval
        end
    else
        x0 = x0 + dx;
    end
   
end

%plot(xlist,flist,eiglist,zeros(size(eiglist)),'*')
end