clear
close all
clc

%% Constant Parameters

g=9.81;                                 % gravity acceleration [m/s^2]
BallMass=0.06157;                       % mass [kg]
BallDiam=0.062;                         % diameter [m]
EMsDist=0.098;                          % distance between magnets [m]
dt=0.001;                               % Simulink step [s]
Drag=0.5*1.293*0.47*pi*(BallDiam/2)^2;  % Sphere drag coeff.

%% RL Identification (Upper Magnet)

fileInductance = {'Induttanza/Induttanza1', 'Induttanza/Induttanza2', 'Induttanza/Induttanza3', 'Induttanza/Induttanza4', 'Induttanza/Induttanza5', 'Induttanza/Induttanza6', 'Induttanza/Induttanza7'};

Res = 4;  % Resistance
V = 11; % Voltage for the experiment (1 duty cycle)

for i=1:7
    experim = load(fileInductance{i});

    position = mean(experim.MLS2EMExpData.signals(1).values);
    curr_ss = mean(experim.MLS2EMExpData.signals(3).values(770:800));
    curr_traj_fit = fit(experim.MLS2EMExpData.time(702:720),experim.MLS2EMExpData.signals(3).values(702:720)','smoothingspline');
    curr_traj_vect = curr_traj_fit(7.01:0.001:7.19);

    for k=1:length(curr_traj_vect)
        if curr_traj_vect(k)>=curr_ss*0.632
            tau = k*0.001;
            break
        end
    end

    L(i) = tau*Res;
    pos_L(i) = position;

end

figure(1)
plot(pos_L,L)
grid on
title('Lsup(x)');
xlabel('Position [m]');
ylabel('Inductance [H]');


L1 = 0.116;
L0x0 = 0.00006;

pos_L_plot = flip(pos_L);
figure(2)
hold on
grid on
plot(pos_L_plot(2:end-1),flip(L(2:end-1)),pos_L_plot(2:end-1),L1+L0x0./pos_L_plot(2:end-1))
title('Inductance L(x)=L0+L1x0/x')
legend('Experimental','Fitted')
xlabel('Position [m]')
ylabel('Inductance [H]')
hold off

plot_pos = [0:0.0001:0.0185];
figure(3)
hold on
grid on
plot(pos_L(2:end-1),L(2:end-1),'o')
plot(plot_pos(20:end),L1+L0x0./plot_pos(20:end))
title('Inductance L(x) (Upper Magnet)')
legend('Experimental points','Hyperbolic fitting')
xlabel('Distance from the magnet [m]')
ylabel('Inductance [H]')
hold off

%% RL Identification (Bottom Magnet)

fileInductance = {'Induttanza2_v1/Induttanza1', 'Induttanza2_v1/Induttanza2', 'Induttanza2_v1/Induttanza3', 'Induttanza2_v1/Induttanza4', 'Induttanza2_v1/Induttanza5', 'Induttanza2_v1/Induttanza6', 'Induttanza2_v1/Induttanza7', 'Induttanza2_v1/Induttanza8', 'Induttanza2_v1/Induttanza9', 'Induttanza2_v1/Induttanza10', 'Induttanza2_v1/Induttanza11', 'Induttanza2_v1/Induttanza12'};

Res = 4;  % Resistance
V = 11;   % Voltage for the experiment (1 duty cycle)

for i=1:12
    experim = load(fileInductance{i});

    position = EMsDist-BallDiam-mean(experim.MLS2EMExpData.signals(1).values(1:300));
    curr_ss = mean(experim.MLS2EMExpData.signals(3).values(770:800,2));
    curr_traj_fit = fit(experim.MLS2EMExpData.time(702:720),experim.MLS2EMExpData.signals(3).values(702:720,2),'smoothingspline');
    curr_traj_vect = curr_traj_fit(7.01:0.001:7.19);

    for k=1:length(curr_traj_vect)
        if curr_traj_vect(k)>=curr_ss*0.632
            tau = k*0.001;
            break
        end
    end

    L2(i) = tau*Res;
    pos_L2(i) = position;
end

L2 = L2([2:8, 10:end-1]);
pos_L2 = pos_L2([2:8, 10:end-1]);

figure(100)
plot(pos_L2,L2)
grid on
title('Linf(x)');
xlabel('Position [m]');
ylabel('Inductance [H]');

L12 = 0.119;
L0x02 = 0.00003;

figure(102)
hold on
grid on
plot(pos_L2,flip(L2),'o')
plot(pos_L2,L12+L0x02./pos_L2)
plot(plot_pos(20:end),L12+L0x02./plot_pos(20:end),'r--')
title('Inductance L(x) (Bottom Magnet)')
legend('Experimental points','Hyperbolic fitting','Non-working region')
xlabel('Distance from the magnet [m]')
ylabel('Inductance [H]')
ylim([0.115 0.15])
hold off


%% PWM to V

PWM = 0:0.1:1;
Volt = [0.0043, 0.61, 1.75, 2.9, 4.06, 5.21, 6.36, 7.51, 8.7, 9.86, 11]; % Obtained with the multimeter

PWM_V_fitted = fit(PWM',Volt','poly1');
PWM_V = PWM_V_fitted(0:0.1:1);

m = PWM_V_fitted.p1;
q = PWM_V_fitted.p2;

figure(4)
hold on
plot(PWM_V_fitted,PWM,m*PWM+q)
plot(PWM,m*PWM)
title('Relationship PWM - Voltage')
grid on
legend('Experimental points','m*x+q','m*x')
xlabel('PWM Duty Cycle')
ylabel('Voltage [V]')

%% PWM to current

PWM = [0.12 0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34 0.35 0.37 0.39 0.41 0.43 0.45 0.46 0.50 0.51 0.54 0.56 0.59 0.6 0.61 0.64 0.66 0.7 0.71];
i_PWM = [0.214 0.341 0.401 0.429 0.505 0.561 0.624 0.667 0.713 0.768 0.809 0.852 0.907 0.953 1.049 1.069 1.133 1.148 1.226 1.254 1.35 1.442 1.51 1.54 1.56 1.66 1.7 1.826 1.872];

PWM_i = fit(PWM',i_PWM','poly1');

%% Magnetic force constant Ksup

% i_dL = [0.214 0.341 0.401 0.429 0.505 0.561 0.624 0.667 0.713 0.768 0.809 0.852 0.907 0.953 1.049 1.069 1.133 1.148 1.226 1.254 1.35 1.442 1.51 1.54 1.56 1.66 1.7 1.826 1.872];
% pos_dL = [0.00016 0.00084 0.00132 0.00196 0.00275 0.00373 0.00442 0.00534 0.00617 0.00692 0.00763 0.00841 0.00902 0.00976 0.01042 0.01108 0.01171 0.01228 0.01295 0.01371 0.01445 0.01504 0.01567 0.01522 0.01611 0.01677 0.0174 0.01788 0.0183];
% 
% for i=1:length(pos_dL)
%     K1(i) = (BallMass*g*pos_dL(i)^2)/i_dL(i)^2;
% end

Kem_1=[6.72059206644115e-05	6.92145376989519e-05	6.92332150285897e-05	6.80659146785482e-05	6.99479714473845e-05	7.05722629122651e-05	6.80128124860478e-05	7.01221702001213e-05	7.05016878442255e-05	7.05515419323440e-05	6.98455307642310e-05	7.13055533313152e-05	6.96274521174565e-05	6.90147402074580e-05	6.99419931300404e-05	6.95805738707539e-05	6.93269258852514e-05	6.90843621353604e-05	7.00006512636940e-05	6.91149414076029e-05	6.97519014539535e-05	6.83573053782350e-05	6.89933518537361e-05	6.79953062557338e-05	6.77206336827605e-05	6.82595333358458e-05	6.6494e-05	6.41548968339966e-05	6.08172215535155e-05	6.05299131311581e-05	6.01818568879911e-05	5.87123846492744e-05	5.77159494392489e-05	5.52986546480238e-05	5.22809928128838e-05	4.60062308129917e-05	4.12885327629886e-05	3.57707157204477e-05	2.88156521526937e-05	2.23928764595895e-05	1.69237291476886e-05];
pos_1=[0.0194279333333333	0.0190477333333333	0.0188606000000000	0.0181366333333333	0.0178862666666667	0.0176649333333333	0.0172163000000000	0.0169433666666667	0.0166807000000000	0.0164495000000000	0.0160222000000000	0.0153453666666667	0.0148285333333333	0.0144294333333333	0.0141447000000000	0.0138273666666667	0.0133592000000000	0.0130921666666667	0.0128335333333333	0.0124948500000000	0.0121420000000000	0.0118023333333333	0.0114197500000000	0.0110494800000000	0.0106465400000000	0.0103677000000000	0.0100000000000000	0.00959026666666667	0.00925280000000000	0.00890150000000000	0.00854740000000000	0.00817620000000000	0.00779310000000000	0.00711953333333333	0.00644370000000000	0.00575746666666667	0.00497293333333333	0.00423913333333333	0.00351336666666667	0.00277350000000000	0.00202256666666667];

figure(5)
hold on
plot(pos_1,interp1(pos_1,Kem_1,pos_1))
title('Magnetic force constanst (Uppper Magnet)')
xlabel('Distance from the magnet [m]')
ylabel('Ksup')
grid on
hold off

%% Magnetic force constant Kinf

pos_2 = [0.0176 0.01529 0.0125 0.0084 0.0026];
curr_Kinf_sup = [2.041 2.041 2.041 1.27 0.76];
curr_Kinf_inf = [0.35 1.45 2.38 1.71 1.11];

xd = EMsDist-BallDiam;

for i=1:length(pos_2)
    Ksup = interp1(pos_1,Kem_1,pos_2(i));
    Kem_2(i) = (Ksup*curr_Kinf_sup(i)^2/pos_2(i)^2-BallMass*g)*(xd-pos_2(i))^2/curr_Kinf_inf(i)^2;
end

pos_2 = [0.0176 0.0166 0.016 0.01529 0.01429 0.0130 0.0125 0.01 0.0084 0.006 0.0026];
Kem_2 = Kem_1(1:length(pos_2));
for i=1:length(Kem_2)
    Kem_2(i)=Kem_2(i)*1.1;
end

figure(50)
hold on
plot(xd-pos_2, Kem_2)
title('Magnetic force constanst (Bottom Magnet)')
xlabel('Distance from the magnet [m]')
ylabel('Kinf')
ylim([1e-5,8e-5])
xlim([0,0.034])
grid on
hold off


%% State-space representation (SISO System)

% Continuous-time

xeq = 0.01
% K1f = Kem_1(round(xeq*length(Kem_1)/0.02));
K1f = interp1(pos_1,Kem_1,xeq);
Lx0 = interp1(pos_L,L,xeq);
ieq = sqrt((BallMass*g)/K1f)*xeq
ueq = (Res*ieq)/m

a21 = (2*K1f*ieq^2)/(BallMass*xeq^3);
a23 = -(2*K1f*ieq)/(BallMass*xeq^2);
% a31 = -((L0x0/xeq^2)*(m*ueq-Res*ieq))/(L1+L0x0/xeq)^2;
a31 = 0;
a33 = -Res/Lx0;
b3 = m/Lx0;

A = [ 0    1    0;
     a21   0   a23;
     a31   0   a33];
B = [0; 0; b3];
C = [1  0  0];
D = 0;

SISOLinearizedMaglev = ss(A,B,C,D)
[SISOnumLinearized,SISOdenLinearized] = ss2tf(A,B,C,D);
SISOLinearizedMaglevTF = tf(SISOnumLinearized,SISOdenLinearized)

% Discrete-time

SISOLinearizedMaglevDiscrete = c2d(SISOLinearizedMaglev,dt,'tustin')
F = SISOLinearizedMaglevDiscrete.A;
G = SISOLinearizedMaglevDiscrete.B;
H = SISOLinearizedMaglevDiscrete.C;

[SISOnumLinearizedDiscr,SISOdenLinearizedDiscr] = ss2tf(F,G,H,0);
SISOLinearizedMaglevTFDiscr = tf(SISOnumLinearizedDiscr,SISOdenLinearizedDiscr,dt)

%% Pole-Placement

% State-feedback controller (continuous-time)
% p = [-50, -55, -40];
p = [-8.01*12, -4.02*12, -4];
K_PP = place(A,B,p);

% State-feedback controller with integrators (continuous-time)
Abar = A;
Abar(4,1) = -1;
Abar(4,4) = 0;
Bbar = B;
Bbar(4,1) = 0;
% pbar = [-50, -55, -40, -30];
% pbar = [-40, -40.5, -35, -10];
pbar = [-4*12, -4.01*12, -4.02*12, -4];
K_PP_integrators = place(Abar,Bbar,pbar);

Cbar = C;
Cbar(1,4) = 0;
PPSysBar = ss(Abar,Bbar,Cbar,D);
% PPSysBarCL = ss((Abar-Bbar*K_PP_integrators),Bbar,Cbar,D);
% pzmap(PPSysBarCL)


% State observer (continuous-time)
% p_obs = [-150, -80.1, -250.2];
p_obs = pbar(1:3);
L_PP_obs = place(A',C',p_obs)';



% State-feedback controller (discrete-time)
pDisc=exp(p.*dt);
K_PP_Disc = place(F,G,pDisc);

% State-feedback controller with integrators (dicrete-time)
PPSysBar_disc = c2d(PPSysBar,dt,'tustin');
Fbar = PPSysBar_disc.A;
Gbar = PPSysBar_disc.B;
pbar_disc = exp(pbar.*dt);
K_PP_integrators_disc = place(Fbar,Gbar,pbar_disc);

% State observer (discrete-time)
% p_obs_disc = exp(p_obs.*dt);
p_obs_disc = exp((10*pbar(1:3)).*dt);
L_PP_obs_disc = place(F',(H*F)',p_obs_disc);

%% LQR

x_d = [xeq; 0; ieq];   % Desired state

Q = [30 0  0;         % Penalize postion error
     0   0.001  0;     % Penalize velocity error
     0   0  10];       % Penalize current error

R = 5.5;              % Penalize control effort

N = 0;
 
[K_LQR,S_LQR,CLP_LQR] = lqr(SISOLinearizedMaglev,Q,R);                  % Continuous-time LQR
[K_LQR_Discr,S_LQR_Discr,CLP_LQR_Discr] = dlqr(F,G,Q,R,N);              % Discrete-time LQR

Q_int = [50   0    0   0;
     0   0.001  0   0;
     0     0   10   0;
     0     0    0   100];

R_int = 1.5;

[K_LQR_int,S_LQR_int,CLP_LQR_int] = lqi(SISOLinearizedMaglev,Q_int,R_int,N);                            % Continuous-time LQR with integrator
[K_LQR_Discr_int,S_LQR_Discr_int,CLP_LQR_Discr_int] = lqi(SISOLinearizedMaglevDiscrete,Q_int,R_int,N);  % Discrete-time LQR with integrator

[K_LQR_2,S_LQR_2,CLP_LQR_2] = lqr(PPSysBar,Q_int,R_int);                        % Continuous-time LQR enlarged
[K_LQR_2_disc,S_LQR_2_disc,CLP_LQR_2_disc] = dlqr(Fbar,Gbar,Q_int,R_int);       % Discrete-time LQR enlarged

K_lqr_r = (C*(eye(3)-(A+B*K_LQR))^(-1)*B)^(-1);


%% Luenberger's observer

% m = 10;
% xeq = 0.01
% % K1f = Kem_1(round(xeq*length(Kem_1)/0.02));
% K1f = interp1(pos_1,Kem_1,xeq);
% Lx0 = interp1(pos_L,L,xeq);
% ieq = sqrt((BallMass*g)/K1f)*xeq
% ueq = (Res*ieq)/m
% 
% a21 = (2*K1f*ieq^2)/(BallMass*xeq^3);
% a23 = -(2*K1f*ieq)/(BallMass*xeq^2);
% % a31 = -((L0x0/xeq^2)*(m*ueq-Res*ieq))/(L1+L0x0/xeq)^2;
% a31 = 0;
% a33 = -Res/Lx0;
% b3 = m/Lx0;
% 
% A = [ 0    1    0;
%      a21   0   a23;
%      a31   0   a33];
% B = [0; 0; b3];
% C = [1  0  0];
% D = 0;
% 
% SISOLinearizedMaglev = ss(A,B,C,D)
% [SISOnumLinearized,SISOdenLinearized] = ss2tf(A,B,C,D);
% SISOLinearizedMaglevTF = tf(SISOnumLinearized,SISOdenLinearized)
% 
% SISOLinearizedMaglevDiscrete = c2d(SISOLinearizedMaglev,dt,'tustin')
% F = SISOLinearizedMaglevDiscrete.A;
% G = SISOLinearizedMaglevDiscrete.B;
% H = SISOLinearizedMaglevDiscrete.C;
% 
% [SISOnumLinearizedDiscr,SISOdenLinearizedDiscr] = ss2tf(F,G,H,0);
% SISOLinearizedMaglevTFDiscr = tf(SISOnumLinearizedDiscr,SISOdenLinearizedDiscr,dt)
% 
% % Observer gain
% p_obs_disc = exp((10*pbar(1:3)).*dt);
% L_PP_obs_disc = place(F',(H*F)',p_obs_disc);

%% State-space representation (Double Magnet)

% % Continuous-time
% 
% xeq = 0.01
% xd = EMsDist-BallDiam;
% K1f = interp1(pos_1,Kem_1,xeq);
% K2f = interp1(xd-pos_2,Kem_2,xd-xeq);
% Lx0_sup = interp1(pos_L,L,xeq);
% Lx0_inf = interp1(pos_L2,L2,xd-xeq);
% 
% o = (K1f*m^2)/(BallMass*Res^2*xeq^2)-(K2f*m^2)/(BallMass*Res^2*(xd-xeq)^2);
% ueq = sqrt(g/o);
% ieq_sup = m*ueq/Res;
% ieq_inf = ieq_sup;
% 
% a21 = (2*K1f*ieq_sup^2)/(BallMass*xeq^3)+(2*K1f*ieq_inf^2)/(BallMass*(xd-xeq)^3);
% a23 = -(2*K1f*ieq_sup)/(BallMass*xeq^2);
% a24 = -(2*K1f*ieq_inf)/(BallMass*(xd-xeq)^2);
% a33 = -Res/Lx0_sup;
% a44 = -Res/Lx0_inf;
% b3 = m/Lx0_sup;
% b4 = m/Lx0_inf;
% 
% 
% A = [ 0    1    0    0;
%      a21   0   a23  a24;
%       0    0   a33   0;
%       0    0    0   a44];
% B = [0;
%      0;
%      b3;
%      b4];
% C = [1  0  0  0];
% D = [0];
% 
% SISOLinearizedMaglev2 = ss(A,B,C,D)
% % [MISOnumLinearized,MISOdenLinearized] = ss2tf(A,B,C,D);
% SISOLinearizedMaglevTF2 = tf(SISOLinearizedMaglev2)
% 
% % Discrete-time
% 
% SISOLinearizedMaglevDiscrete2 = c2d(SISOLinearizedMaglev2,dt,'tustin')
% F = SISOLinearizedMaglevDiscrete2.A;
% G = SISOLinearizedMaglevDiscrete2.B;
% H = SISOLinearizedMaglevDiscrete2.C;
% 
% % [MISOnumLinearizedDiscr,MISOdenLinearizedDiscr] = ss2tf(F,G,H,[0 0]);
% SISOLinearizedMaglevTFDiscr2 = tf(SISOLinearizedMaglevDiscrete2)

%% Double Magnet Controllers

% % Pole-Placement
% p2 = [-3, -3.*12, -3.01*12, -3.02*12];
% K_PP = place(A,B,p2);
% pDisc2=exp(p2.*dt);
% K_PP_Disc = place(F,G,pDisc2);
% 
% Abar = A;
% Abar(5,1) = -1;
% Abar(5,5) = 0;
% Bbar = B;
% Bbar(5,1) = 0;
% Cbar = C;
% Cbar(1,5) = 0;
% 
% pbar = p2;
% pbar(5) = pbar(4)*1.01;
% K_PP_integrators = place(Abar,Bbar,pbar);
% PPSysBar = ss(Abar,Bbar,Cbar,D);
% PPSysBar_disc = c2d(PPSysBar,dt,'tustin');
% Fbar = PPSysBar_disc.A;
% Gbar = PPSysBar_disc.B;
% pbar_disc = exp(pbar.*dt);
% K_PP_integrators_disc = place(Fbar,Gbar,pbar_disc);
% 
% % LQR
% 
% Q_int = diag([50; 0; 0; 0; 100]);
% R_int = 1.5;
% [K_LQR_int,S_LQR_int,CLP_LQR_int] = lqi(SISOLinearizedMaglev2,Q_int,R_int,N);  

%% Cose

s = tf('s');

Lslqr = K_LQR_int*inv(s*eye(4)-Abar)*Bbar;
Tslqr = Lslqr/(1+Lslqr);

Ls = K_LQR*inv(s*eye(3)-A)*B;
Ts = Ls/(1+Ls);

Lspp = K_PP_integrators*inv(s*eye(4)-Abar)*Bbar;
Tspp = Lspp/(1+Lspp);

% Non attendibile in quanto sistema instabile
figure
hold on
margin(Tspp)
margin(Tslqr)
title('Complementary Sensitivity Function')
legend('PP','LQR')
hold off

% Banda reale
bandwidth(Tspp)
bandwidth(Tslqr)
