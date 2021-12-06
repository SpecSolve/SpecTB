function [spec,dist,maxdist] = CompSpec(A,I,varargin)
% This code approximates the spectrum of the operator A using the algorithm
% of M.J. Colbrook, B. Roman, and A.C. Hansen, "How to compute spectra with
% error control", Physical Review Letters 122 (2019). The code below is
% written for operators given by sparse self-adjoint matrices. For other
% settings (including speed ups with orderings of the basis), see the paper
% and other papers about the Solvability Complexity Index.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   A: the operator truncated to suitable rectangular matrix
%   I: the discrete grid of points where we approximate the spectrum
% OPTIONAL LABELLED INPUTS
%   Tol: tolerance for svds
%   Parallel: parfor (on) or normal for (off) loop, default is "off"
%   Prog: display progress of the algorithm, default is "on". NB: turn this
%   off when running multiple instances of the algorithm at the same time
% OUTPUTS
%   spec: approximation of the spectrum
%   dist: the distance function computed over I
%   maxdist: the max of E1 over spec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect the optional inputs
p = inputParser;
validPar = {'on','off'};
checkPar = @(x) any(validatestring(x,validPar));
addParameter(p,'Tol',0.1*(I(2)-I(1)),@(x) x>0);
addParameter(p,'Parallel','off',checkPar)
addParameter(p,'Prog','on',checkPar)
p.CaseSensitive = false;
parse(p,varargin{:})

% Compute the singular values
B=speye(size(A,1),size(A,2));
dist=zeros(length(I),1);

if p.Results.Prog=="on" % display progress
    pf = parfor_progress(length(I));
    pfcleanup = onCleanup(@() delete(pf));
end
if p.Results.Parallel=="off"
    for j=1:length(I)
        dist(j)=svds(A-I(j)*B,1,'smallest','Tolerance',p.Results.Tol);
        if p.Results.Prog=="on" % display progress
            parfor_progress(pf);
        end
    end
else % perform the computation in parallel
    tol=p.Results.Tol;
    prog=p.Results.Prog;
    parfor j=1:length(I)
        dist(j)=svds(A-I(j)*B,1,'smallest','Tolerance',tol,'MaxIterations',10000);
        if prog=="on" % display progress
            parfor_progress(pf);
        end
    end
end

% Perform the local search routine
spec=[];
maxdist=0;
for j=1:length(I)
    if dist(j)<0.5
        J=find(abs(I-I(j))<=dist(j));
        a=I(J);
        b=dist(J);
        maxdist=max(maxdist,min(b));
        spec=union(spec(:),a(b==min(b)));
    end
end
end

