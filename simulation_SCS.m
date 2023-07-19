clear all;
close all;

%% Part 1: Defining parameters and constants
%start and end time of pulses for different fibers
tON_AB = 0; tOFF_AB = 0.01; % sec
tON_AD = 0.02; tOFF_AD = 0.035; % sec 
tON_C = 0.1; tOFF_C = 0.175; % sec

AB_max = 2; AD_max = 0.5; C_max = 1.5;  % maximum amplitude of different fiber stimulations for pain
dt=0.001; % step interval to be used for ODE
gABW = 0.6; gABI=0.6; gADW=0.3; gCE=0.6; gCW=0.4; gEW=0.4; gIE=0.05; gIW=0.15;     % weight of synaptic connections between Excitatory, Inhibitory and Projection neurons 
Te = 0.02; Ti = 0.02; Tw = 0.001; %ms    % intrinsic time scale for different neuron populations

%% Part 2: Generating input to dorsal horn circuit
fAB = zeros(1,0.2/dt);   % pre-allocation of firing neuron variables
fAD = zeros(1,0.2/dt);
fC = zeros(1,0.2/dt);
t1 = 0.001:0.001:0.2;
for i=1:1:0.2/dt
    t=i*dt;
    fAB(i) = AB_max*heaviside(t-tON_AB)*heaviside(tOFF_AB-t);
    fAD(i) = AD_max*heaviside(t-tON_AD)*heaviside(tOFF_AD-t);
    fC(i) = C_max*heaviside(t-tON_C)*heaviside(tOFF_C-t);
end

figure(1);
plot(t1,fAB, 'LineWidth', 1.5); hold on
plot(t1,fAD, 'LineWidth', 1.5); hold on
plot(t1,fC, 'LineWidth', 1.5); hold on
title('Sensory afferent input to dorsal horn circuit')
%title('Simulated model input to the dorsal horn circuit from the afferent fibers after pre-processing in the dorsal root ganglion')
xlabel('Time (s)')
ylabel('Sensory Afferent Input')
legend('Abeta', 'Adelta', 'C')
hold off

%% Part 3: Defining firing rate response functions
N=3;
E_max = 60; I_max = 80; W_max = 50; %Hz
Winf = zeros(1,N*1000);            % pre-allocating  variables
Einf = zeros(1,N*1000);
Iinf = zeros(1,N*1000);
t=0.001:0.001:N;
for i=2:1:N*1000
    c=i*dt;
    Winf(i) = W_max*0.5*(1+tanh((c-1.2)/0.6));
    Einf(i) = E_max*0.5*(1+tanh((c-1.3)/0.35));
    Iinf(i) = I_max*0.5*(1+tanh((c-1)/0.8));
end
figure(2);
plot(t,Winf, 'LineWidth', 1.5); hold on
plot(t,Einf, 'LineWidth', 1.5); hold on
plot(t,Iinf, 'LineWidth', 1.5);
ylabel('Firing rate (Hz)')
xlabel('Input')
title('Input-Output Curves')
legend('Winf', 'Einf', 'Iinf')
hold off


%% Part 4: Simulating pain response in different nerve populations
T=0.4;   %s
fW = zeros(1,T/dt);
fE = zeros(1,T/dt);
fI = zeros(1,T/dt);

t1 = dt:dt:T;
tmp1=zeros(1,T/dt);
tmp2=zeros(1,T/dt);
tmp3=zeros(1,T/dt);
C_max = 0.5; AD_max = 2.5;
for i=1:1:T*1000-1
    t=i*dt;
    fAB = AB_max*heaviside(t-tON_AB)*heaviside(tOFF_AB-t);
    fAD = AD_max*heaviside(t-tON_AD)*heaviside(tOFF_AD-t);
    fC = C_max*heaviside(t-tON_C)*heaviside(tOFF_C-t);

    c1 = gABW*fAB + gADW*fAD + gCW*fC + gEW*fE(i) - gIW*fI(i);
    Winf = W_max*0.5*(1+tanh((c1-0.1)/1));    %  b = 1.2828
    dfW = dt*(Winf-fW(i))/Tw;
    
    c2 = gCE*fC - gIE*fI(i);
    Einf = E_max*0.5*(1+tanh((c2-0.4267)/0.3));  % b = 0.4267
    dfE = dt*(Einf-fE(i))/Te;
    
    c3 = gABI*fAB;
    Iinf = I_max*0.5*(1+tanh((c3-0.6)/0.45));   %  b = 0.6
    dfI = dt*(Iinf-fI(i))/Ti;

    fW(i+1)=fW(i)+dfW;
    fE(i+1)=fE(i)+dfE;
    fI(i+1)=fI(i)+dfI;
    tmp1(i)=c1;
    tmp2(i)=c2;
    tmp3(i)=c3;
end

fig=figure(3);
subplot(3,1,1)
plot(t1,fW); %hold on
title('Projection neurons (W)')
subplot(3,1,2)
plot(t1,fE); %hold on
title('Excitatory neurons (E)')
subplot(3,1,3)
plot(t1,fI);
title('Inhibitory neurons (I)')

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Firing Rate (Hz)');
xlabel(han,'Time (s)');
%  title(han,'Delayed Inhibition');


%% Part 5: Simulating pain inhibition phenomenon

clear;
AB_max = 2; 
dt=0.001;
T=0.25;   %s, total time of simulation
t1 = dt:dt:T;
N=5; 

fig = figure(4);
for x=1:1:N
    fAB = [AB_max*ones(1,10) zeros(1,25*x) AB_max*ones(1,10) zeros(1,400-20-25*x)];
    subplot(N,1,x); plot(t1, run_sim(fAB, T)); xline((10+25*x)*dt,'--r');
end

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Firing Rate (Hz)');
xlabel(han,'Time (s)');
title(han,'Delayed Inhibition in Projection population (W)');


function fW = run_sim(fAB, T)

    tON_AD = 0.02; tOFF_AD = 0.035;
    tON_C = 0.1; tOFF_C = 0.175;
%     AD_max = 0.5; C_max = 1.5;

    gABW = 0.6; gABI=0.6; gADW=0.3; gCE=0.6; gCW=0.4; gEW=0.4; gIE=0.05; gIW=0.15;
    Te = 0.02; Ti = 0.02; Tw = 0.001; %ms
    E_max = 60; %Hz
    I_max = 80; %Hz
    W_max = 50; %Hz

    dt=0.001;
    fW = zeros(1,T/dt);
    fE = zeros(1,T/dt);
    fI = zeros(1,T/dt);
    C_max = 0.5; AD_max = 2.5;
    for i=1:1:T*1000-1
        t=i*dt;
    %     fAB = AB_max*heaviside(t-tON_AB)*heaviside(tOFF_AB-t);
        fAD = AD_max*heaviside(t-tON_AD)*heaviside(tOFF_AD-t);
        fC = C_max*heaviside(t-tON_C)*heaviside(tOFF_C-t);
    
        c1 = gABW*fAB(i) + gADW*fAD + gCW*fC + gEW*fE(i) - gIW*fI(i);
        Winf = W_max*0.5*(1+tanh((c1-0.5)/0.4));
        dfW = dt*(Winf-fW(i))/Tw;
     
        c2 = gCE*fC - gIE*fI(i);
        Einf = E_max*0.5*(1+tanh((c2-1.2)/0.30));
        dfE = dt*(Einf-fE(i))/Te;
        
        c3 = gABI*fAB(i);
        Iinf = I_max*0.5*(1+tanh((c3-1)/0.45));
        dfI = dt*(Iinf-fI(i))/Ti;
    
        fW(i+1)=fW(i)+dfW;
        fE(i+1)=fE(i)+dfE;
        fI(i+1)=fI(i)+dfI;
    end
end




