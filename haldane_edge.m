function [H,pos,mIndex,nIndex,f] = haldane_edge(N_trunc,t,M,tprime,PHI)
% Computes full edge Hamiltonian s.t.
% INPUT
% N_trunc = number of basis points
% t,M,tprime,PHI the physical parameters
% OUTPUT
% H = Hamiltonian
% pos = vector of positions on hexagonal lattice (complex format)
% mIndex = vector of m index for each basis
% nIndex = vector of n index for each basis
%% construct the relevant submatrices
load('Pre_computed_mats/haldane_lattice.mat','f','n1','n2','pos','ord')
N=(2*n1+1)*(2*n2+1)*2;

d1=spdiags(ones(2*n1,1),-1,2*n1+1,2*n1+1); % m to m-1
d2=spdiags(ones(2*n2,1),-1,2*n2+1,2*n2+1); % n to n-1

H1 = [sparse(round(N/2),round(N/2)), speye(round(N/2))+kron(speye(2*n2+1),d1)+kron(d2,speye(2*n1+1)); % n1 in second kron input DELIBERATE (dont change)
      speye(round(N/2))+kron(speye(2*n2+1),d1')+kron(d2',speye(2*n1+1)),sparse(round(N/2),round(N/2))];

H2 = [-speye(round(N/2)), sparse(round(N/2),round(N/2));
      sparse(round(N/2),round(N/2)), speye(round(N/2))];

H3 = [kron(speye(2*n2+1),d1)+kron(d2',speye(2*n1+1))+kron(d2,d1'), sparse(round(N/2),round(N/2));
      sparse(round(N/2),round(N/2)),kron(speye(2*n2+1),d1')+kron(d2,speye(2*n1+1))+kron(d2',d1)];

%% construct the Hamiltonian
H=t*H1+M*H2+tprime*(exp(1i*PHI)*H3+exp(-1i*PHI)*H3');
H=H(ord,ord);

%% truncate position vectors
mIndex=kron(zeros(1,2*n2+1)+1,-n1:n1);
mIndex=[mIndex(:);mIndex(:)];
mIndex=mIndex(ord);

nIndex=kron(-n2:n2,zeros(1,2*n1+1)+1);
nIndex=[nIndex(:);nIndex(:)];
nIndex=nIndex(ord);

I=find(mIndex>0);
pos(I)=[];
mIndex(I)=[];
nIndex(I)=[];
H(I,:)=[]; H(:,I)=[];

pos=pos(1:N_trunc);
mIndex=mIndex(1:N_trunc);
nIndex=nIndex(1:N_trunc);
H=H(1:f(N_trunc),1:N_trunc);

f=0*f; % compute new f
for j=1:N_trunc
    if j==1
        f(j)=max([max(find(abs(H(1:end,j))>0)),j]);
    else
        f(j)=max([max(find(abs(H(1:end,j))>0)),f(j-1),j]);
    end
end
H=H(1:f(N_trunc),1:N_trunc);
end