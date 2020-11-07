function [K,rho,feas]=LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc)
% Computes, using LMIs, the distributed "state feedback" control law for the continuous-time system, with reference to the control
% information structure specified by 'ContStruc'.
%
% Inputs:
% - Atot: system matrix.
% - Bdec: input matrices (i.e., Bdec{1},..., Bdec{N} are the input matrices of the decomposed system, one for each channel).
% - Cdec: output matrices  (i.e., Cdec{1},..., Cdec{N} are the output matrices of the decomposed system, one for each channel, where [Cdec{1}',...,
% Cdec{N}']=I).
% - N: number of subsystems.
% - ContStruc: NxN matrix that specifies the information structure
% constraints (ContStruc(i,j)=1 if communication is allowed between channel
% j to channel i, ContStruc(i,j)=0 otherwise).
%
% Output:
% - K: structured control gain
% - rho: spectral abscissa of matrix (A+B*K) - note that [Cdec{1}',...,
% Cdec{N}']=I
% - feas: feasibility of the LMI problem (=0 if yes)

Btot=[];
for i=1:N
    m(i)=size(Bdec{i},2);
    n(i)=size(Cdec{i},1);
    Btot=[Btot,Bdec{i}];
end
ntot=size(Atot,1);
mtot=sum(m);

yalmip clear

if ContStruc==ones(N,N)
    % Centralized design
    Y=sdpvar(ntot);
    L=sdpvar(mtot,ntot);
else
    % Decentralized/distributed design
    Y=[];
    L=sdpvar(mtot,ntot);
    minc=0;
    for i=1:N
        Y=blkdiag(Y,sdpvar(n(i)));
        ninc=0;
        for j=1:N
            if ContStruc(i,j)==0
                L(minc+1:minc+m(i),ninc+1:ninc+n(j))=zeros(m(i),n(j));
            end
            ninc=ninc+n(j);
        end
        minc=minc+m(i);
    end  
end

LMIconstr=[Y*Atot'+Atot*Y+Btot*L+L'*Btot'<-1e-2*eye(ntot)]+[Y>1e-2*eye(ntot)];
options=sdpsettings('solver','sedumi');
J=solvesdp(LMIconstr,[],options);
feas=J.problem;
L=double(L);
Y=double(Y);

K=L/Y;
rho=max(real(eig(Atot+Btot*K)));
