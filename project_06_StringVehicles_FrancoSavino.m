clc
clear all

% String of vehicles
% Marcello Farina, 19/12/2018

Atot=[0 1 0 0 0 0 0 0
    0 -1 0 0 0 0 0 0
    0 -1 0 1 0 0 0 0
    0 0 0 -1 0 0 0 0
    0 0 0 -1 0 1 0 0
    0 0 0 0 0 -1 0 0
    0 0 0 0 0 -1 0 1
    0 0 0 0 0 0 0 -1];
Bdec{1}=[0 1 0 0 0 0 0 0]';
Bdec{2}=[0 0 0 1 0 0 0 0]';
Bdec{3}=[0 0 0 0 0 1 0 0]';
Bdec{4}=[0 0 0 0 0 0 0 1]';
Ctot=eye(8);
for i = 1:4
    Cdec{i}=Ctot((i-1)*2+1:2*i,:);
end
 
%--------------------------------------------------%
%1. Generate system matrices in CT and DT


h = 0.1;                  %DT sampling time

Btot = [[0 1 0 0 0 0 0 0]', [0 0 0 1 0 0 0 0]', [0 0 0 0 0 1 0 0]', [0 0 0 0 0 0 0 1]'];

systemCT = ss(Atot,Btot,Ctot,0);
systemDT = c2d(systemCT,h);
[F,G,H,L,h] = ssdata(systemDT);
for i = 1:4
    Gdec{i} = G(:,i);
end
 
Hdec = Cdec;
 
% a. eigenvalues of CT
 
eigenCT = eig(Atot);
spec_abscCT = max(eigenCT);
 
% b. eigenvalues of DT
 
eigenDT = eig(F);
spec_abscDT = max(abs(eigenDT));
 
%--------------------------------------------% 
  
CentrContStruc = ones(4,4);                 %Centralized
 
DecentrContStruc = diag(ones(4,1));         %Decentralized
  
%%%%%%%% STRUTTURE DI CONTROLLO DISTRIBUITO 
       % - 1) ( NODO STELLA caso unidirezionale e bidirezionale );
       
            DistrContStarUnidirect = [1 0 0 0;    %star with unidirectional communication
                                      1 1 0 0;
                                      1 0 1 0;
                                      1 0 0 1];
       
            DistrContStarBidirect = [1 1 1 1;     %star with bidirectional communication
                                     1 1 0 0;
                                     1 0 1 0;
                                     1 0 0 1];
       
       % - 2) ( Ogni veicolo passa info al successivo: -2.1) unidirezionale e -2.2) bidirezionale );
            
            CycleContStruc_Unidir = [1 0 0 1;
                                     1 1 0 0;
                                     0 1 1 0;
                                     0 0 1 1];
            
            CycleContStruc_Bidir =  [1 1 0 1;
                                     1 1 1 0;
                                     0 1 1 1;
                                     1 0 1 1];
                                 
       
       % - 3) ( Ogni veicolo passa info al precedente, unidirezionale );
         
            Inv_CycleContStruc_Unidir = [1 1 0 0;
                                         0 1 1 0;
                                         0 0 1 1;
                                         1 0 0 1];
         
         
% a. CT Fixed modes
  
[CFM_CT] = di_fixed_modes(Atot,Bdec,Cdec,4,CentrContStruc,3);
[DFM_CT] = di_fixed_modes(Atot,Bdec,Cdec,4,DecentrContStruc,3);

[DSUM_CT] = di_fixed_modes(Atot,Bdec,Cdec,4,DistrContStarUnidirect,3);
[DSBM_CT] = di_fixed_modes(Atot,Bdec,Cdec,4,DistrContStarBidirect,3);

[CycleUnidir_FM_CT] = di_fixed_modes(Atot,Bdec,Cdec,4,CycleContStruc_Unidir,3);
[CycleBidir_FM_CT] = di_fixed_modes(Atot,Bdec,Cdec,4,CycleContStruc_Bidir,3);

[Inv_CycleUnidir_FM_CT] = di_fixed_modes(Atot,Bdec,Cdec,4,Inv_CycleContStruc_Unidir,3);
  
 
% b. DT Fixed Modes
  
% There are no CFM and DFM in CT so we do not compute in DT
  
% [CFM_DT] = di_fixed_modes(F,Gdec,Hdec,4,CentrContStruc,3);
% [DFM_DT] = di_fixed_modes(F,Gdec,Hdec,4,DecentrContStruc,3);

% [DSUM_DT] = di_fixed_modes(F,Gdec,Hdec,4,DistrContStarUnidirect,3);
% [DSBM_DT] = di_fixed_modes(F,Gdec,Hdec,4,DistrContStarBidirect,3);

% [CycleUnidir_FM_DT] = di_fixed_modes(F,Gdec,Hdec,4,CycleContStruc_Unidir,3);
% [CycleBidir_FM_DT] = di_fixed_modes(F,Gdec,Hdec,4,CycleContStruc_Bidir,3);

% [Inv_CycleUnidir_FM_DT] = di_fixed_modes(F,Gdec,Hdec,4,Inv_CycleContStruc_Unidir,3);  
 
  
% c. Compute CT control gain
  
[K_Centr_CT, spec_abscCT_Centr_closed, feas_Centr_CT] = LMI_CT_DeDicont(Atot,Bdec,Cdec,4,CentrContStruc);          %Centralized
[K_Decentr_CT, spec_abscCT_Decentr_closed, feas_Decentr_CT] = LMI_CT_DeDicont(Atot,Bdec,Cdec,4,DecentrContStruc);  %Decentralized

[K_Dist_Star_Uni_CT,spec_abscCT_Dist_Star_Uni_closed,feas_Dist_Star_Uni_CT] = LMI_CT_DeDicont(Atot,Bdec,Cdec,4,DistrContStarUnidirect);   %distributed Star Unidirectional  
[K_Dist_Star_Bi_CT,spec_abscCT_Dist_Star_Bi_closed,feas_Dist_Star_Bi_CT] = LMI_CT_DeDicont(Atot,Bdec,Cdec,4,DistrContStarBidirect);       %distributed Star Bidirectional

[K_CycleUnidir_CT, spec_abscCT_CycleUnidir_closed, feas_CycleUnidir_CT] = LMI_CT_DeDicont(Atot,Bdec,Cdec,4,CycleContStruc_Unidir);  %Cycle Unidirectional 
[K_CycleBidir_CT, spec_abscCT_CycleBidir_closed, feas_CycleBidir_CT] = LMI_CT_DeDicont(Atot,Bdec,Cdec,4,CycleContStruc_Bidir);      %Cycle Bidirectional

[K_Inv_CycleUnidir_CT, spec_abscCT_Inv_CycleUnidir_closed, feas_Inv_CycleUnidir_CT] = LMI_CT_DeDicont(Atot,Bdec,Cdec,4,Inv_CycleContStruc_Unidir);  %Ineverted Cycle Unidirectional
  
% d. Compute DT control gain
  
[K_Centr_DT, spec_abscDT_Centr_closed, feas_Centr_DT] = LMI_DT_DeDicont(F,Gdec,Hdec,4,CentrContStruc);          %Centralized
[K_Decentr_DT, spec_abscDT_Decentr_closed, feas_Decentr_DT] = LMI_DT_DeDicont(F,Gdec,Hdec,4,DecentrContStruc);  %Decentralized

[K_Dist_Star_Uni_DT,spec_abscDT_Dist_Star_Uni_closed,feas_Dist_Star_Uni_DT] = LMI_DT_DeDicont(F,Gdec,Hdec,4,DistrContStarUnidirect);  %distributed Star Unidirectional  
[K_Dist_Star_Bi_DT,spec_abscDT_Dist_Star_Bi_closed,feas_Dist_Star_Bi_DT] = LMI_DT_DeDicont(F,Gdec,Hdec,4,DistrContStarBidirect);      %distributed Star Bidirectional

[K_CycleUnidir_DT, spec_abscDT_CycleUnidir_closed, feas_CycleUnidir_DT] = LMI_DT_DeDicont(F,Gdec,Hdec,4,CycleContStruc_Unidir);  %Cycle Unidirectional
[K_CycleBidir_DT, spec_abscDT_CycleBidir_closed, feas_CycleBidir_DT] = LMI_DT_DeDicont(F,Gdec,Hdec,4,CycleContStruc_Bidir);      %Cycle Bidirectional

[K_Inv_CycleUnidir_DT, spec_abscDT_Inv_CycleUnidir_closed, feas_Inv_CycleUnidir_DT] = LMI_DT_DeDicont(F,Gdec,Hdec,4,Inv_CycleContStruc_Unidir);     %Ineverted Cycle Unidirectional

%-------------------------------------------------------%

% e. Closed loop system trajectories

x0 = [5 0 5 0 5 0 5 0]';              %Init. position
xref = [0 0 0 0 0 0 0 0]';            %Reference position

x_Centr = x0;
x_Decentr = x0;

x_Dist_Star_Uni = x0;
x_Dist_Star_Bi = x0;

x_CycleUnidir = x0;
x_CycleBidir = x0;

x_Inv_CycleUnidir = x0;



A_Centr_CL = Atot+Btot*K_Centr_CT;        %CL centralized CT
F_Centr_CL = F+G*K_Centr_DT;              %CL centralized DT

A_Decentr_CL = Atot+Btot*K_Decentr_CT;    %CL decentralized CT
F_Decentr_CL = F+G*K_Decentr_DT;          %CL decentralized DT

A_Dist_Star_Uni_CL = Atot+Btot*K_Dist_Star_Uni_CT;    %CL distributed star unidirectional CT
F_Dist_Star_Uni_CL = F+G*K_Dist_Star_Uni_DT;          %CL distributed star unidirectional DT

A_Dist_Star_Bi_CL = Atot+Btot*K_Dist_Star_Bi_CT;      %CL distributed star bidirectional CT
F_Dist_Star_Bi_CL = F+G*K_Dist_Star_Bi_DT;            %CL distributed star bidirectional DT

A_CycleUnidir_CL = Atot+Btot*K_CycleUnidir_CT;      %CL Cycle Unidirectional CT
F_CycleUnidir_CL = F+G*K_CycleUnidir_DT;            %CL Cycle Unidirectional DT

A_CycleBidir_CL = Atot+Btot*K_CycleBidir_CT;        %CL Cycle Unidirectional CT
F_CycleBidir_CL = F+G*K_CycleBidir_DT;              %CL Cycle Unidirectional DT

A_Inv_CycleUnidir_CL = Atot+Btot*K_Inv_CycleUnidir_CT;    %CL Inverted Cycle Unidirectional CT
F_Inv_CycleUnidir_CL = F+G*K_Inv_CycleUnidir_DT;          %CL Inverted Cycle Unidirectional DT


T_final = 20;
h_c = 0.01;                                   %CT plotting step
t = 0 : h_c : T_final;                        %CT timeline
steps = 0 : T_final/h;                        %DT timeline


%DT Closed loops

for k = 1 :(length(t)-1)*h
    
    x_Centr(:,k+1) = F_Centr_CL* x_Centr(:,k) - G*K_Centr_DT*xref;
    x_Decentr(:,k+1) = F_Decentr_CL*x_Decentr(:,k) - G*K_Decentr_DT*xref;
    
    x_Dist_Star_Uni(:,k+1) = F_Dist_Star_Uni_CL* x_Dist_Star_Uni(:,k) - G*K_Dist_Star_Uni_DT*xref;
    x_Dist_Star_Bi(:,k+1) = F_Dist_Star_Bi_CL* x_Dist_Star_Bi(:,k) - G*K_Dist_Star_Bi_DT*xref;
    
    x_CycleUnidir(:,k+1) = F_CycleUnidir_CL* x_CycleUnidir(:,k) - G*K_CycleUnidir_DT*xref;
    x_CycleBidir(:,k+1) = F_CycleBidir_CL* x_CycleBidir(:,k) - G*K_CycleBidir_DT*xref;

    x_Inv_CycleUnidir(:,k+1) = F_Inv_CycleUnidir_CL* x_Inv_CycleUnidir(:,k) - G*K_Inv_CycleUnidir_DT*xref;
end





%% RELATIVE TRAJECTORIES
%----------------------------------------------------%



figure ('Name', 'Closed loop CENTRALIZED control relative trajectories')



subplot(2,2,1)

plot(out.V1_centralized.Time,out.V1_centralized.Data)
hold on 
grid on
plot(out.V2_centralized.Time,out.V2_centralized.Data)
plot(out.V3_centralized.Time,out.V3_centralized.Data)
plot(out.V4_centralized.Time,out.V4_centralized.Data)


title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')

legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

for k  = 1:2:7
    hold on
    plot(steps*h,x_Centr(k,:),'o')
end
grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')

subplot(2,2,3:4)

plot(out.V1_centralized.Time,out.V1_centralized.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_centralized.Time,out.V2_centralized.Data,'LineWidth',2)
plot(out.V3_centralized.Time,out.V3_centralized.Data,'LineWidth',2)
plot(out.V4_centralized.Time,out.V4_centralized.Data,'LineWidth',2)
plot(steps*h,x_Centr(1,:),'b.')
plot(steps*h,x_Centr(3,:),'r.')
plot(steps*h,x_Centr(5,:),'y.')
plot(steps*h,x_Centr(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')




figure ('Name', 'Closed loop DECENTRALIZED control relative trajectories')



subplot(2,2,1)

plot(out.V1_decentralized.Time,out.V1_decentralized.Data)
hold on 
grid on
plot(out.V2_decentralized.Time,out.V2_decentralized.Data)
plot(out.V3_decentralized.Time,out.V3_decentralized.Data)
plot(out.V4_decentralized.Time,out.V4_decentralized.Data)

title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')
grid on
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

for k = 1:2:7
    hold on
    plot(steps*h,x_Decentr(k,:),'o')
end
grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,3:4)

plot(out.V1_decentralized.Time,out.V1_decentralized.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_decentralized.Time,out.V2_decentralized.Data,'LineWidth',2)
plot(out.V3_decentralized.Time,out.V3_decentralized.Data,'LineWidth',2)
plot(out.V4_decentralized.Time,out.V4_decentralized.Data,'LineWidth',2)
plot(steps*h,x_Decentr(1,:),'b.')
plot(steps*h,x_Decentr(3,:),'r.')
plot(steps*h,x_Decentr(5,:),'y.')
plot(steps*h,x_Decentr(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')




figure ('Name', 'Closed loop DISTRIBUTED STAR WITH UNIDIRECTIONAL COMMUNICATION control relative trajectories')

subplot(2,2,1)

plot(out.V1_dist_uni.Time,out.V1_dist_uni.Data)
hold on 
grid on
plot(out.V2_dist_uni.Time,out.V2_dist_uni.Data)
plot(out.V3_dist_uni.Time,out.V3_dist_uni.Data)
plot(out.V4_dist_uni.Time,out.V4_dist_uni.Data)

title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')
grid on
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

for k=1:2:7
    hold on
    plot(steps*h,x_Dist_Star_Uni(k,:),'o')
end
grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,3:4)

plot(out.V1_dist_uni.Time,out.V1_dist_uni.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_dist_uni.Time,out.V2_dist_uni.Data,'LineWidth',2)
plot(out.V3_dist_uni.Time,out.V3_dist_uni.Data,'LineWidth',2)
plot(out.V4_dist_uni.Time,out.V4_dist_uni.Data,'LineWidth',2)
plot(steps*h,x_Dist_Star_Uni(1,:),'b.')
plot(steps*h,x_Dist_Star_Uni(3,:),'r.')
plot(steps*h,x_Dist_Star_Uni(5,:),'y.')
plot(steps*h,x_Dist_Star_Uni(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')




figure ('Name', 'Closed loop DISTRIBUTED STAR WITH BIDIRECTIONAL COMMUNICATION control relative trajectories')

subplot(2,2,1)

plot(out.V1_dist_bi.Time,out.V1_dist_bi.Data)
hold on 
grid on
plot(out.V2_dist_bi.Time,out.V2_dist_bi.Data)
plot(out.V3_dist_bi.Time,out.V3_dist_bi.Data)
plot(out.V4_dist_bi.Time,out.V4_dist_bi.Data)

title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')
grid on
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

for k=1:2:7
    hold on
    plot(steps*h,x_Dist_Star_Bi(k,:),'o')
end
grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,3:4)

plot(out.V1_dist_bi.Time,out.V1_dist_bi.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_dist_bi.Time,out.V2_dist_bi.Data,'LineWidth',2)
plot(out.V3_dist_bi.Time,out.V3_dist_bi.Data,'LineWidth',2)
plot(out.V4_dist_bi.Time,out.V4_dist_bi.Data,'LineWidth',2)
plot(steps*h,x_Dist_Star_Bi(1,:),'b.')
plot(steps*h,x_Dist_Star_Bi(3,:),'r.')
plot(steps*h,x_Dist_Star_Bi(5,:),'y.')
plot(steps*h,x_Dist_Star_Bi(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')





figure ('Name', 'Closed loop CYCLE UNIDIRECTIONAL control relative trajectories')



subplot(2,2,1)

plot(out.V1_cycle_uni.Time,out.V1_cycle_uni.Data)
hold on 
grid on
plot(out.V2_cycle_uni.Time,out.V2_cycle_uni.Data)
plot(out.V3_cycle_uni.Time,out.V3_cycle_uni.Data)
plot(out.V4_cycle_uni.Time,out.V4_cycle_uni.Data)


title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')

legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

for k = 1:2:7
    hold on
    plot(steps*h,x_CycleUnidir(k,:),'o')
end
grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')

subplot(2,2,3:4)

plot(out.V1_cycle_uni.Time,out.V1_cycle_uni.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_cycle_uni.Time,out.V2_cycle_uni.Data,'LineWidth',2)
plot(out.V3_cycle_uni.Time,out.V3_cycle_uni.Data,'LineWidth',2)
plot(out.V4_cycle_uni.Time,out.V4_cycle_uni.Data,'LineWidth',2)
plot(steps*h,x_CycleUnidir(1,:),'b.')
plot(steps*h,x_CycleUnidir(3,:),'r.')
plot(steps*h,x_CycleUnidir(5,:),'y.')
plot(steps*h,x_CycleUnidir(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')




figure ('Name', 'Closed loop CYCLE BIDIRECTIONAL control relative trajectories')



subplot(2,2,1)

plot(out.V1_cycle_bi.Time,out.V1_cycle_bi.Data)
hold on 
grid on
plot(out.V2_cycle_bi.Time,out.V2_cycle_bi.Data)
plot(out.V3_cycle_bi.Time,out.V3_cycle_bi.Data)
plot(out.V4_cycle_bi.Time,out.V4_cycle_bi.Data)


title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')

legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

for k = 1:2:7
    hold on
    plot(steps*h,x_CycleBidir(k,:),'o')
end
grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')

subplot(2,2,3:4)

plot(out.V1_cycle_bi.Time,out.V1_cycle_bi.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_cycle_bi.Time,out.V2_cycle_bi.Data,'LineWidth',2)
plot(out.V3_cycle_bi.Time,out.V3_cycle_bi.Data,'LineWidth',2)
plot(out.V4_cycle_bi.Time,out.V4_cycle_bi.Data,'LineWidth',2)
plot(steps*h,x_CycleBidir(1,:),'b.')
plot(steps*h,x_CycleBidir(3,:),'r.')
plot(steps*h,x_CycleBidir(5,:),'y.')
plot(steps*h,x_CycleBidir(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')




figure ('Name', 'Closed loop INVERTED CYCLE UNIDIRECTIONAL control relative trajectories')



subplot(2,2,1)

plot(out.V1_inv_cycle_uni.Time,out.V1_inv_cycle_uni.Data)
hold on 
grid on
plot(out.V2_inv_cycle_uni.Time,out.V2_inv_cycle_uni.Data)
plot(out.V3_inv_cycle_uni.Time,out.V3_inv_cycle_uni.Data)
plot(out.V4_inv_cycle_uni.Time,out.V4_inv_cycle_uni.Data)


title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')

legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

for k = 1:2:7
    hold on
    plot(steps*h,x_Inv_CycleUnidir(k,:),'o')
end
grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')

subplot(2,2,3:4)

plot(out.V1_inv_cycle_uni.Time,out.V1_inv_cycle_uni.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_inv_cycle_uni.Time,out.V2_inv_cycle_uni.Data,'LineWidth',2)
plot(out.V3_inv_cycle_uni.Time,out.V3_inv_cycle_uni.Data,'LineWidth',2)
plot(out.V4_inv_cycle_uni.Time,out.V4_inv_cycle_uni.Data,'LineWidth',2)
plot(steps*h,x_Inv_CycleUnidir(1,:),'b.')
plot(steps*h,x_Inv_CycleUnidir(3,:),'r.')
plot(steps*h,x_Inv_CycleUnidir(5,:),'y.')
plot(steps*h,x_Inv_CycleUnidir(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')






%% ABSOLUTE TRAJECTORIES
%-------------------------------------------------------%




figure ('Name', 'Closed loop CENTRALIZED control absolute trajectories')



subplot(2,2,1)

plot(out.V1_centralized.Time,out.V1_centralized.Data)
hold on 
grid on
plot(out.V2_centralized.Time,out.V1_centralized.Data+out.V2_centralized.Data)
plot(out.V3_centralized.Time,out.V1_centralized.Data+out.V2_centralized.Data+out.V3_centralized.Data)
plot(out.V4_centralized.Time,out.V1_centralized.Data+out.V2_centralized.Data+out.V3_centralized.Data+out.V4_centralized.Data)


title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')

legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

plot(steps*h,x_Centr(1,:),'o')
hold on
plot(steps*h,x_Centr(1,:)+x_Centr(3,:),'o')
plot(steps*h,x_Centr(1,:)+x_Centr(3,:)+x_Centr(5,:),'o')
plot(steps*h,x_Centr(1,:)+x_Centr(3,:)+x_Centr(5,:)+x_Centr(7,:),'o')
grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')

subplot(2,2,3:4)

plot(out.V1_centralized.Time,out.V1_centralized.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_centralized.Time,out.V1_centralized.Data+out.V2_centralized.Data,'LineWidth',2)
plot(out.V3_centralized.Time,out.V1_centralized.Data+out.V2_centralized.Data+out.V3_centralized.Data,'LineWidth',2)
plot(out.V4_centralized.Time,out.V1_centralized.Data+out.V2_centralized.Data+out.V3_centralized.Data+out.V4_centralized.Data,'LineWidth',2)
plot(steps*h,x_Centr(1,:),'b.')
plot(steps*h,x_Centr(1,:)+x_Centr(3,:),'r.')
plot(steps*h,x_Centr(1,:)+x_Centr(3,:)+x_Centr(5,:),'y.')
plot(steps*h,x_Centr(1,:)+x_Centr(3,:)+x_Centr(5,:)+x_Centr(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')




figure ('Name', 'Closed loop DECENTRALIZED control absolute trajectories')



subplot(2,2,1)

plot(out.V1_decentralized.Time,out.V1_decentralized.Data)
hold on 
grid on
plot(out.V2_decentralized.Time,out.V1_decentralized.Data+out.V2_decentralized.Data)
plot(out.V3_decentralized.Time,out.V1_decentralized.Data+out.V2_decentralized.Data+out.V3_decentralized.Data)
plot(out.V4_decentralized.Time,out.V1_decentralized.Data+out.V2_decentralized.Data+out.V3_decentralized.Data+out.V4_decentralized.Data)

title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')
grid on
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

plot(steps*h,x_Decentr(1,:),'o')
hold on
plot(steps*h,x_Decentr(1,:)+x_Decentr(3,:),'o')
plot(steps*h,x_Decentr(1,:)+x_Decentr(3,:)+x_Decentr(5,:),'o')
plot(steps*h,x_Decentr(1,:)+x_Decentr(3,:)+x_Decentr(5,:)+x_Decentr(7,:),'o')

grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,3:4)

plot(out.V1_decentralized.Time,out.V1_decentralized.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_decentralized.Time,out.V1_decentralized.Data+out.V2_decentralized.Data,'LineWidth',2)
plot(out.V3_decentralized.Time,out.V1_decentralized.Data+out.V2_decentralized.Data+out.V3_decentralized.Data,'LineWidth',2)
plot(out.V4_decentralized.Time,out.V1_decentralized.Data+out.V2_decentralized.Data+out.V3_decentralized.Data+out.V4_decentralized.Data,'LineWidth',2)
plot(steps*h,x_Decentr(1,:),'b.')
plot(steps*h,x_Decentr(1,:)+x_Decentr(3,:),'r.')
plot(steps*h,x_Decentr(1,:)+x_Decentr(3,:)+x_Decentr(5,:),'y.')
plot(steps*h,x_Decentr(1,:)+x_Decentr(3,:)+x_Decentr(5,:)+x_Decentr(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')




figure ('Name', 'Closed loop DISTRIBUTED STAR WITH UNIDIRECTIONAL COMMUNICATION control absolute trajectories')

subplot(2,2,1)

plot(out.V1_dist_uni.Time,out.V1_dist_uni.Data)
hold on 
grid on
plot(out.V2_dist_uni.Time,out.V1_dist_uni.Data+out.V2_dist_uni.Data)
plot(out.V3_dist_uni.Time,out.V1_dist_uni.Data+out.V2_dist_uni.Data+out.V3_dist_uni.Data)
plot(out.V4_dist_uni.Time,out.V1_dist_uni.Data+out.V2_dist_uni.Data+out.V3_dist_uni.Data+out.V4_dist_uni.Data)

title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')
grid on
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

plot(steps*h,x_Dist_Star_Uni(1,:),'o')
hold on
plot(steps*h,x_Dist_Star_Uni(1,:)+x_Dist_Star_Uni(3,:),'o')
plot(steps*h,x_Dist_Star_Uni(1,:)+x_Dist_Star_Uni(3,:)+x_Dist_Star_Uni(5,:),'o')
plot(steps*h,x_Dist_Star_Uni(1,:)+x_Dist_Star_Uni(3,:)+x_Dist_Star_Uni(5,:)+x_Dist_Star_Uni(7,:),'o')
grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,3:4)

plot(out.V1_dist_uni.Time,out.V1_dist_uni.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_dist_uni.Time,out.V1_dist_uni.Data+out.V2_dist_uni.Data,'LineWidth',2)
plot(out.V3_dist_uni.Time,out.V1_dist_uni.Data+out.V2_dist_uni.Data+out.V3_dist_uni.Data,'LineWidth',2)
plot(out.V4_dist_uni.Time,out.V1_dist_uni.Data+out.V2_dist_uni.Data+out.V3_dist_uni.Data+out.V4_dist_uni.Data,'LineWidth',2)
plot(steps*h,x_Dist_Star_Uni(1,:),'b.')
plot(steps*h,x_Dist_Star_Uni(1,:)+x_Dist_Star_Uni(3,:),'r.')
plot(steps*h,x_Dist_Star_Uni(1,:)+x_Dist_Star_Uni(3,:)+x_Dist_Star_Uni(5,:),'y.')
plot(steps*h,x_Dist_Star_Uni(1,:)+x_Dist_Star_Uni(3,:)+x_Dist_Star_Uni(5,:)+x_Dist_Star_Uni(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')




figure ('Name', 'Closed loop DISTRIBUTED STAR WITH BIDIRECTIONAL COMMUNICATION control absolute trajectories')

subplot(2,2,1)

plot(out.V1_dist_bi.Time,out.V1_dist_bi.Data)
hold on 
grid on
plot(out.V2_dist_bi.Time,out.V1_dist_bi.Data+out.V2_dist_bi.Data)
plot(out.V3_dist_bi.Time,out.V1_dist_bi.Data+out.V2_dist_bi.Data+out.V3_dist_bi.Data)
plot(out.V4_dist_bi.Time,out.V1_dist_bi.Data+out.V2_dist_bi.Data+out.V3_dist_bi.Data+out.V4_dist_bi.Data)

title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')
grid on
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

plot(steps*h,x_Dist_Star_Bi(1,:),'o')
hold on
plot(steps*h,x_Dist_Star_Bi(1,:)+x_Dist_Star_Bi(3,:),'o')
plot(steps*h,x_Dist_Star_Bi(1,:)+x_Dist_Star_Bi(3,:)+x_Dist_Star_Bi(5,:),'o')
plot(steps*h,x_Dist_Star_Bi(1,:)+x_Dist_Star_Bi(3,:)+x_Dist_Star_Bi(5,:)+x_Dist_Star_Bi(7,:),'o')
grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,3:4)

plot(out.V1_dist_bi.Time,out.V1_dist_bi.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_dist_bi.Time,out.V1_dist_bi.Data+out.V2_dist_bi.Data,'LineWidth',2)
plot(out.V3_dist_bi.Time,out.V1_dist_bi.Data+out.V2_dist_bi.Data+out.V3_dist_bi.Data,'LineWidth',2)
plot(out.V4_dist_bi.Time,out.V1_dist_bi.Data+out.V2_dist_bi.Data+out.V3_dist_bi.Data+out.V4_dist_bi.Data,'LineWidth',2)
plot(steps*h,x_Dist_Star_Bi(1,:),'b.')
plot(steps*h,x_Dist_Star_Bi(1,:)+x_Dist_Star_Bi(3,:),'r.')
plot(steps*h,x_Dist_Star_Bi(1,:)+x_Dist_Star_Bi(3,:)+x_Dist_Star_Bi(5,:),'y.')
plot(steps*h,x_Dist_Star_Bi(1,:)+x_Dist_Star_Bi(3,:)+x_Dist_Star_Bi(5,:)+x_Dist_Star_Bi(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')





figure ('Name', 'Closed loop CYCLE UNIDIRECTIONAL control absolute trajectories')



subplot(2,2,1)

plot(out.V1_cycle_uni.Time,out.V1_cycle_uni.Data)
hold on 
grid on
plot(out.V2_cycle_uni.Time,out.V1_cycle_uni.Data+out.V2_cycle_uni.Data)
plot(out.V3_cycle_uni.Time,out.V1_cycle_uni.Data+out.V2_cycle_uni.Data+out.V3_cycle_uni.Data)
plot(out.V4_cycle_uni.Time,out.V1_cycle_uni.Data+out.V2_cycle_uni.Data+out.V3_cycle_uni.Data+out.V4_cycle_uni.Data)


title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')

legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

plot(steps*h,x_CycleUnidir(1,:),'o')
hold on
plot(steps*h,x_CycleUnidir(1,:)+x_CycleUnidir(3,:),'o')
plot(steps*h,x_CycleUnidir(1,:)+x_CycleUnidir(3,:)+x_CycleUnidir(5,:),'o')
plot(steps*h,x_CycleUnidir(1,:)+x_CycleUnidir(3,:)+x_CycleUnidir(5,:)+x_CycleUnidir(7,:),'o')
grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')

subplot(2,2,3:4)

plot(out.V1_cycle_uni.Time,out.V1_cycle_uni.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_cycle_uni.Time,out.V1_cycle_uni.Data+out.V2_cycle_uni.Data,'LineWidth',2)
plot(out.V3_cycle_uni.Time,out.V1_cycle_uni.Data+out.V2_cycle_uni.Data+out.V3_cycle_uni.Data,'LineWidth',2)
plot(out.V4_cycle_uni.Time,out.V1_cycle_uni.Data+out.V2_cycle_uni.Data+out.V3_cycle_uni.Data+out.V4_cycle_uni.Data,'LineWidth',2)
plot(steps*h,x_CycleUnidir(1,:),'b.')
plot(steps*h,x_CycleUnidir(1,:)+x_CycleUnidir(3,:),'r.')
plot(steps*h,x_CycleUnidir(1,:)+x_CycleUnidir(3,:)+x_CycleUnidir(5,:),'y.')
plot(steps*h,x_CycleUnidir(1,:)+x_CycleUnidir(3,:)+x_CycleUnidir(5,:)+x_CycleUnidir(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')




figure ('Name', 'Closed loop CYCLE BIDIRECTIONAL control absolute trajectories')



subplot(2,2,1)

plot(out.V1_cycle_bi.Time,out.V1_cycle_bi.Data)
hold on 
grid on
plot(out.V2_cycle_bi.Time,out.V1_cycle_bi.Data+out.V2_cycle_bi.Data)
plot(out.V3_cycle_bi.Time,out.V1_cycle_bi.Data+out.V2_cycle_bi.Data+out.V3_cycle_bi.Data)
plot(out.V4_cycle_bi.Time,out.V1_cycle_bi.Data+out.V2_cycle_bi.Data+out.V3_cycle_bi.Data+out.V4_cycle_bi.Data)


title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')

legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

plot(steps*h,x_CycleBidir(1,:),'o')
hold on
plot(steps*h,x_CycleBidir(1,:)+x_CycleBidir(3,:),'o')
plot(steps*h,x_CycleBidir(1,:)+x_CycleBidir(3,:)+x_CycleBidir(5,:),'o')
plot(steps*h,x_CycleBidir(1,:)+x_CycleBidir(3,:)+x_CycleBidir(5,:)+x_CycleBidir(7,:),'o')
grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')

subplot(2,2,3:4)

plot(out.V1_cycle_bi.Time,out.V1_cycle_bi.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_cycle_bi.Time,out.V1_cycle_bi.Data+out.V2_cycle_bi.Data,'LineWidth',2)
plot(out.V3_cycle_bi.Time,out.V1_cycle_bi.Data+out.V2_cycle_bi.Data+out.V3_cycle_bi.Data,'LineWidth',2)
plot(out.V4_cycle_bi.Time,out.V1_cycle_bi.Data+out.V2_cycle_bi.Data+out.V3_cycle_bi.Data+out.V4_cycle_bi.Data,'LineWidth',2)
plot(steps*h,x_CycleBidir(1,:),'b.')
plot(steps*h,x_CycleBidir(1,:)+x_CycleBidir(3,:),'r.')
plot(steps*h,x_CycleBidir(1,:)+x_CycleBidir(3,:)+x_CycleBidir(5,:),'y.')
plot(steps*h,x_CycleBidir(1,:)+x_CycleBidir(3,:)+x_CycleBidir(5,:)+x_CycleBidir(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')




figure ('Name', 'Closed loop INVERTED CYCLE UNIDIRECTIONAL control absolute trajectories')



subplot(2,2,1)

plot(out.V1_inv_cycle_uni.Time,out.V1_inv_cycle_uni.Data)
hold on 
grid on
plot(out.V2_inv_cycle_uni.Time,out.V1_inv_cycle_uni.Data+out.V2_inv_cycle_uni.Data)
plot(out.V3_inv_cycle_uni.Time,out.V1_inv_cycle_uni.Data+out.V2_inv_cycle_uni.Data+out.V3_inv_cycle_uni.Data)
plot(out.V4_inv_cycle_uni.Time,out.V1_inv_cycle_uni.Data+out.V2_inv_cycle_uni.Data+out.V3_inv_cycle_uni.Data+out.V4_inv_cycle_uni.Data)


title ('Continuous time trajectories')
xlabel('t')
ylabel('distance')

legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')


subplot(2,2,2)

plot(steps*h,x_Inv_CycleUnidir(1,:),'o')
hold on
plot(steps*h,x_Inv_CycleUnidir(1,:)+x_Inv_CycleUnidir(3,:),'o')
plot(steps*h,x_Inv_CycleUnidir(1,:)+x_Inv_CycleUnidir(3,:)+x_Inv_CycleUnidir(5,:),'o')
plot(steps*h,x_Inv_CycleUnidir(1,:)+x_Inv_CycleUnidir(3,:)+x_Inv_CycleUnidir(5,:)+x_Inv_CycleUnidir(7,:),'o')
grid on
title ('Discrete time trajectories')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')

subplot(2,2,3:4)

plot(out.V1_inv_cycle_uni.Time,out.V1_inv_cycle_uni.Data,'LineWidth',2)
hold on 
grid on
plot(out.V2_inv_cycle_uni.Time,out.V1_inv_cycle_uni.Data+out.V2_inv_cycle_uni.Data,'LineWidth',2)
plot(out.V3_inv_cycle_uni.Time,out.V1_inv_cycle_uni.Data+out.V2_inv_cycle_uni.Data+out.V3_inv_cycle_uni.Data,'LineWidth',2)
plot(out.V4_inv_cycle_uni.Time,out.V1_inv_cycle_uni.Data+out.V2_inv_cycle_uni.Data+out.V3_inv_cycle_uni.Data+out.V4_inv_cycle_uni.Data,'LineWidth',2)
plot(steps*h,x_Inv_CycleUnidir(1,:),'b.')
plot(steps*h,x_Inv_CycleUnidir(1,:)+x_Inv_CycleUnidir(3,:),'r.')
plot(steps*h,x_Inv_CycleUnidir(1,:)+x_Inv_CycleUnidir(3,:)+x_Inv_CycleUnidir(5,:),'y.')
plot(steps*h,x_Inv_CycleUnidir(1,:)+x_Inv_CycleUnidir(3,:)+x_Inv_CycleUnidir(5,:)+x_Inv_CycleUnidir(7,:),'m.')
title ('Comparison')
xlabel('t')
ylabel('distance')
legend('vehicle 1', 'vehicle 2', 'vehicle 3', 'vehicle 4')