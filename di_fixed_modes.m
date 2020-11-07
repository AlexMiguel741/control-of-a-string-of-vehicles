function [Difm]=di_fixed_modes(Atot,Baggr,Caggr,N,ContStruc,rounding_n)
% Computes the fixed modes of the system, with reference to the control
% information structure specified by 'ContStruc'.
%
% Inputs:
% - Atot: system matrix - either continuous-time or discrete-time.
% - Baggr: input matrices - either continuous-time or discrete-time -
% (i.e., Baggr{1},..., Baggr{N} are the input matrices of the decomposed system, one for each channel).
% - Caggr: output matrices - either continuous-time or discrete-time -
% (i.e., Caggr{1},..., Caggr{N} are the output matrices of the decomposed
% system, one for each channel).
% - N: number of subsystems.
% - ContStruc: NxN matrix that specifies the information structure
% constraints (ContStruc(i,j)=1 if communication is allowed between channel
% j to channel i, ContStruc(i,j)=0 otherwise).
% - rounding_n: the eigenvalues are rounded at the nth decimal place 
%
% Output:
% - Difm: vector containing the fixed modes (if empty, there are no FMs)

Btot=[];
Ctot=[];
for i=1:N
    m(i)=size(Baggr{i},2);
    p(i)=size(Caggr{i},1);
    Btot=[Btot,Baggr{i}];
    Ctot=[Ctot
        Caggr{i}];
end

m_tot=size(Btot,2);
p_tot=size(Ctot,1);

Difm=round(eig(Atot),rounding_n);
nelD=length(Difm);
 
kend=1000;
k=0;
while (nelD~=0)&&(k<=kend)
    k=k+1;
    K=zeros(m_tot,p_tot);
    m_inc=0;
    for i=1:N
        p_inc=0;
        for j=1:N
            if ContStruc(i,j)~=0
                K(m_inc+1:m_inc+m(i),p_inc+1:p_inc+p(j))=100*randn(m(i),p(j));
            end
            p_inc=p_inc+p(j);
        end
        m_inc=m_inc+m(i);
    end
    eF=round(eig(Atot+Btot*K*Ctot),rounding_n);
    C=intersect(Difm,eF);
    Difm=C;
    nelD=length(Difm);
end

