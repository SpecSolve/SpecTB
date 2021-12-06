function [P] = pvm2(A,b,mode_range,nnodes,epsilon,f,varargin)
% This code applies an approximate spectral projector of a finite (possibly
% rectangular) matrix A, associated with real interval [a,b] to vector b. 
% To apply to infinite discrete matrices, rectangular truncations must be 
% used so that A is f(n)\times n, where f describes the sparsity structure 
% of off diagonal decay.
% NB: REQUIRES TWICE AS MANY SOLVES AS THE SCALAR REAL VALUED CASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% A: finite rectangular matrix (sparse or full)
% b: vector in \mathbb{C}^{f(n)}
% mode_range: interval [a b] input to PVM
% nnodes: number of quadrature nodes, weights for approximate projection
% epsilon: the smoothing parameter
% f: sparsity function for A

% OPTIONAL LABELLED INPUTS
% PoleType: where to place poles, default is equispaced
% Parallel: parfor (on) or normal for (off) loop, default is "off"
% Order: order of kernel used, default is 2
% DiscMin: minimum truncation size to begin adaptive truncation
% DiscMax: maximum truncation size used during adaptive truncation
% ContType: Integration contour for smoothed PVM, default is 'deform'

% OUTPUTS
% P: approximate spectral projection of b onto 'states' assoc. with
%    the interval specified in mode_range

% REQUIRES
% Chebfun: https://www.chebfun.org/ (we use the command chebpts() below)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect the optional inputs
p = inputParser;
addRequired(p,'A',@isnumeric);
addRequired(p,'b',@isnumeric);
addRequired(p,'mode_range',@isnumeric);
addRequired(p,'nnodes',@isnumeric);
addRequired(p,'epsilon',@isnumeric);
addRequired(p,'f',@isnumeric);

validPole = {'cheb','roots','extrap','equi'};
checkPole = @(x) any(validatestring(x,validPole));
validPar = {'on','off'};
checkPar = @(x) any(validatestring(x,validPar));
validCont = {'flat','deform'};
checkCont = @(x) any(validatestring(x,validCont));
validEndpts = {'include','exclude'};
checkEndpts = @(x) any(validatestring(x,validEndpts));

addParameter(p,'DiscMin',min(2^10,size(A,2)),@(x) x==floor(x))
addParameter(p,'DiscMax',size(A,2),@(x) x==floor(x))
addParameter(p,'PoleType','equi',checkPole)
addParameter(p,'Parallel','off',checkPar)
addParameter(p,'ContType','deform',checkCont)
addParameter(p,'Order',2,@(x) x==floor(x))
addParameter(p,'EndPts','include',checkEndpts)

p.CaseSensitive = false;
parse(p,A,b,mode_range,nnodes,epsilon,f,varargin{:})

%% SpecSolve projection algorithm

%compute quadrature weights and nodes
if p.Results.ContType == "deform"
    [nodes,weights]=legpts(nnodes,[0 1]);
    z0=mean(mode_range); r=diff(mode_range)/2;
    nodes=z0+r*exp(1i*pi*nodes); weights=1i*pi*(nodes-z0).*transpose(weights);
elseif p.Results.ContType=="flat"
    [nodes,weights]=legpts(nnodes,mode_range);
else
    fprintf('Invalid argument in ContType\n')
end

%include or exclude endpoint contributions
if p.Results.EndPts=="include"
    weightEP=[epsilon/2 epsilon/2];
elseif p.Results.EndPts=="exclude"
    weightEP=[-epsilon/2 -epsilon/2];
end

%compute poles and residues for rational kernel
[poles,res]=rational_kernel(p.Results.Order,p.Results.PoleType);
tol=epsilon^(p.Results.Order+1); X=nodes;

%compute resolvent adaptively at quadrature nodes
N=p.Results.DiscMin;
errIndx=1:length(X);
P=zeros(N,1); mu=zeros(N,length(X));
while ~isempty(errIndx) && N<=p.Results.DiscMax
    
    AT=A(1:f(N),1:N); bt=b(1:f(N)); I=speye(size(AT));  %truncate A and b
    temp1=zeros(N,length(errIndx)); temp1(1:size(mu,1),1:size(mu,2))=mu;
    temp2=zeros(N,1); temp2(1:length(P))=P; mu=zeros(N,length(errIndx));
    if issparse(AT)
        if p.Results.Parallel=="off"
            for j=1:length(X(errIndx))
                a=zeros(N,1);
                for mm=1:p.Results.Order
                    a=a+weights(errIndx(j))*res(mm)*((AT-(X(errIndx(j))-epsilon*poles(mm))*I)\bt);
                    a=a-conj(weights(errIndx(j))*res(mm))*((AT-conj(X(errIndx(j))-epsilon*poles(mm))*I)\bt);
                end
                mu(:,errIndx(j))=a/(2*pi*1i);
            end
            temp3=errIndx; errIndx=find(vecnorm(mu-temp1)>tol*vecnorm(mu));
            goodIndx=setdiff(temp3,errIndx); P=temp2-sum(mu(:,goodIndx),2);
        else
            XTemp=X(errIndx); WTemp=weights(errIndx); muTemp=mu(:,errIndx);
            parfor j=1:length(XTemp)
                a=zeros(N,1);
                for mm=1:p.Results.Order
                    a=a+WTemp(j)*res(mm)*((AT-(XTemp(j)-epsilon*poles(mm))*I)\bt);
                    a=a-conj(WTemp(j)*res(mm))*((AT-conj(XTemp(j)-epsilon*poles(mm))*I)\bt);
                end
                muTemp(:,j)=a/(2*pi*1i);
            end
            mu(:,errIndx)=muTemp;
            temp3=errIndx; errIndx=find(vecnorm(mu-temp1)>tol*vecnorm(mu));
            goodIndx=setdiff(temp3,errIndx); P=temp2-sum(mu(:,goodIndx),2);
        end
    else
        if p.Results.Parallel=="off"
            for j=1:length(errIndx)
                a=zeros(N,1);
                for mm=1:p.Results.Order
                    a=a+weights(errIndx(j))*res(mm)*((AT-(X(errIndx(j))-epsilon*poles(mm))*I)\bt);
                    a=a-conj(weights(errIndx(j))*res(mm))*((AT-conj(X(errIndx(j))-epsilon*poles(mm))*I)\bt);
                end
                mu(:,errIndx(j))=a/(2*pi*1i);
            end
            temp3=errIndx; errIndx=find(vecnorm(mu-temp1)>tol*vecnorm(mu));
            goodIndx=setdiff(temp3,errIndx); P=temp2-sum(mu(:,goodIndx),2);
        else
            XTemp=X(errIndx); WTemp=weights(errIndx); muTemp=mu(:,errIndx);
            parfor j=1:length(errIndx)
                a=zeros(N,1);
                for mm=1:p.Results.Order
                    a=a+WTemp(j)*res(mm)*((AT-(XTemp(j)-epsilon*poles(mm))*I)\bt);
                    a=a-conj(WTemp(j)*res(mm))*((AT-conj(XTemp(j)-epsilon*poles(mm))*I)\bt);
                end
                muTemp(:,j)=a/(2*pi*1i);
            end
            mu(:,errIndx)=muTemp;
            temp3=errIndx; errIndx=find(vecnorm(mu-temp1)>tol*vecnorm(mu));
            goodIndx=setdiff(temp3,errIndx); P=temp2-sum(mu(:,goodIndx),2);
        end
    end
    N=3*N;
end
P=P-sum(mu(:,errIndx),2);
temp4=zeros(size(A,2),1); temp4(1:length(P))=P; P=temp4;

% add endpoint contributions using Poisson w/ max discretization size
I=speye(size(A));
endpt1=(A-(mode_range(1)-epsilon*1i)*I)\b-(A-(mode_range(1)+epsilon*1i)*I)\b;
endpt2=(A-(mode_range(2)-epsilon*1i)*I)\b-(A-(mode_range(2)+epsilon*1i)*I)\b;
P=P-weightEP(1)*endpt1/2i-weightEP(2)*endpt2/2i;
end