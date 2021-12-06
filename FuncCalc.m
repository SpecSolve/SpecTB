function [zi,vi] = FuncCalc(A,psi,a,b,delta,varargin)
% This function computes nodes zi and vectors vi such that f(A)psi is
% approximately \sum_{i}f(zi)vi for suitable holomorphic functions f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   A: the operator truncated to suitable rectangular matrix
%   psi: initial vector
%   a: min real part along contour
%   b: max real part along contour
%   delta: distance to real axis above and below spectrum
% OPTIONAL LABELLED INPUTS
%   N1: number of quadrature points for vertical contours
%   N2: number of quadrature points for horizontal contours
%   Parallel: parfor (on) or normal for (off) loop, default is "off"
%   Prog: display progress of the algorithm, default is "on". NB: turn this
%   off when running multiple instances of the algorithm at the same time
% OUTPUTS
%   zi: quadrature nodes
%   vi: quadrature weights (vectors)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: currently not set up to provide error control.
% Collect the optional inputs
p = inputParser;
validPar = {'on','off'};
checkPar = @(x) any(validatestring(x,validPar));
addParameter(p,'N1',ceil(10*max(1,delta)),@(x) x==round(x));
addParameter(p,'N2',ceil(10*(b-a)),@(x) x==round(x));
addParameter(p,'Parallel','off',checkPar)
addParameter(p,'Prog','on',checkPar)
p.CaseSensitive = false;
parse(p,varargin{:})
   
% Perform the quadrature
N1=p.Results.N1; N2=p.Results.N2;
zi=zeros(1,2*(N1+N2)); vi=zeros(size(A,2),2*(N1+N2));

if p.Results.Prog=="on" % display progress
    pf = parfor_progress(2*(N1+N2));
    pfcleanup = onCleanup(@() delete(pf));
else
    pf=[];
end

[zi(1:N1),vi(:,1:N1)]=line_int2(a-delta*1i,a+delta*1i,A,psi,N1,p,pf);
[zi(N1+1:(N1+N2)),vi(:,N1+1:(N1+N2))]=line_int2(a+delta*1i,b+delta*1i,A,psi,N2,p,pf);
[zi((N1+N2+1):(2*N1+N2)),vi(:,(N1+N2+1):(2*N1+N2))]=line_int2(b+delta*1i,b-delta*1i,A,psi,N1,p,pf);
[zi((2*N1+N2+1):end),vi(:,(2*N1+N2+1):end)]=line_int2(b-delta*1i,a-delta*1i,A,psi,N2,p,pf);
end

function  [x,w]=lgwt2(N,a,b)
[x,w]=legpts(N);
% Linear map from[-1,1] to [a,b]
x=a+(b-a)*(x+1)/2;    
% Compute the weights
w=(b-a)/2*w;
end

function [z,V]=line_int2(z0,z1,A,v0,N,p,pf)
    [x,w]=lgwt2(N,0,abs(z1-z0));
    z=z0+(z1-z0)*transpose(x(:))/abs(z1-z0);
    V=zeros(size(A,2),N);
    if p.Results.Parallel=="off"
        for j=1:length(x)
            V(:,j)=((A-z(j)*speye(size(A,1),size(A,2)))\v0)*w(j)/(2*pi*1i)*(z1-z0)/abs(z1-z0);
            if p.Results.Prog=="on" % display progress
                parfor_progress(pf);
            end
        end
    else
        prog=p.Results.Prog;
        parfor j=1:length(x)
            V(:,j)=((A-z(j)*speye(size(A,1),size(A,2)))\v0)*w(j)/(2*pi*1i)*(z1-z0)/abs(z1-z0);
            if prog=="on" % display progress
                parfor_progress(pf);
            end
        end
    end
end

