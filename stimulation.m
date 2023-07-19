clear;


gABW = 0.6;
gABI=0.6;
gADW=0.3;
gCE=0.6;
gCW=0.4;
gEW=0.4;
gIE=0.05;
gIW=0.15;
AB_max = 2;
AD_max = 0.5;
C_max = 1.5;

Te = 20; %ms
Ti = 20; %ms
Tw = 1; %ms

T = 0.2; %sec , total time of simulation
dt = 0.001; % 1ms
 
tON_AB = 0.0011;
tOFF_AB = 0.01;
tON_AD = 0.02;
tOFF_AD = 0.03;
tON_C = 0.1;
tOFF_C = 0.175;

fW = zeros(1,T/dt);
fE = zeros(1,T/dt);
fI = zeros(1,T/dt);

for t=0.002:dt:T
    i=int8(t*1000);
    fAB = AB_max*heaviside(t-tON_AB)*heaviside(tOFF_AB-t);
    fAD = AD_max*heaviside(t-tON_AD)*heaviside(tOFF_AD-t);
    fC = C_max*heaviside(t-tON_C)*heaviside(tOFF_C-t);
    c1=gABW*fAB+gADW*fAD+gCW*fC+gEW*fE(i)-gIW*fI(i);
    dfW = dt*(firing_rate(c1,80)-fW(i))/Tw;
    c2=gCE*fC-gIE*fI(i);
    dfE = dt*(firing_rate(c2,60)-fE(i))/Te;
    c3 = gABI*fAB;
    dfI = dt*(firing_rate(c3,50)-fI(i))/Ti;
    
    
    fW(i)=fW(i-1)+dfW;
    fE(i)=fE(i-1)+dfE;
    fI(i)=fI(i-1)+dfI;

end

% t1 = 0.001:0.001:T;
% plot(t1,fW); hold on;
% plot(t1, fE); hold on;
% plot(t1, fI); hold off

%% 
E_max = 60; %Hz
I_max = 80; %Hz
W_max = 50; %Hz
rw = W_max/2;  % 25
ew = E_max/2;  % 30
iw = I_max/2;  % 40
N=3;
t1 = 0.001:0.001:N;
Winf = zeros(1,N*1000);
Einf = zeros(1,N*1000);
Iinf = zeros(1,N*1000);
dt=0.001;
for i=1:1:N*1000
    c = i*dt;
    Winf(i) = rw*(1+tanh(c-1.2));
    Einf(i) = ew*(1+tanh(c-1.3));
    Iinf(i) = iw*(1+tanh(c-1));
end
plot(t1,Winf); hold on
plot(t1,Einf); hold on
plot(t1,Iinf);
legend('Winf', 'Einf', 'Iinf')
title('Response functions') 

%% Plotting Fig3 (input to model)

fAB = zeros(1,0.2/dt);
fAD = zeros(1,0.2/dt);
fC = zeros(1,0.2/dt);
t1 = 0.001:0.001:0.2;
for t=0.001:0.001:0.2
%     t*1000, heaviside(t-tON_AB)
    fAB(int8(t*1000)) = AB_max*heaviside(t-tON_AB)*heaviside(tOFF_AB-t);
    fAD(int8(t*1000)) = AD_max*heaviside(t-tON_AD)*heaviside(tOFF_AD-t);
    fC(int8(t*1000)) = C_max*heaviside(t-tON_C)*heaviside(tOFF_C-t);
end
figure;
plot(t1,fAB); hold on
plot(t1,fAD); hold on
plot(t1,fC);
title('Simulated model input to the dorsal horn circuit from the afferent fibers after pre-processing in the dorsal root ganglion')
xlabel('Time (s)')
ylabel('Sensory Afferent Input')
legend('Abeta', 'Adelta', 'C')
hold off

% Fig 5
T=0.4;
fW = zeros(1,T/dt-1);
fE = zeros(1,T/dt-1);
fI = zeros(1,T/dt-1);

fAB = [fAB zeros(1,199)]; fAD = [fAD zeros(1,199)]; fC = [fC zeros(1,199)];
for t=0.002:dt:T
    i=int8(t*1000);

    c1=gABW*fAB(i)+gADW*fAD(i)+gCW*fC(i)+gEW*fE(i)-gIW*fI(i);
    dfW = dt*(firing_rate(c1,80)-fW(i))/Tw;
    c2=gCE*fC(i)-gIE*fI(i);
    dfE = dt*(firing_rate(c2,60)-fE(i))/Te;
    c3 = gABI*fAB(i);
    dfI = dt*(firing_rate(c3,50)-fI(i))/Ti;
    fW(i)=fW(i-1)+dfW;
    fE(i)=fE(i-1)+dfE; 
    fI(i)=fI(i-1)+dfI;

end
% figure;
% t1 = 0.002:0.001:T;
% plot(t1,fW); hold on;
% plot(t1, fE); hold on;
% plot(t1, fI); hold off

%%

c=normrnd(1.25,0.75,8200);
ad=normrnd(0.12,0.083,900);
ab=normrnd(0.024,0.013,900);
% c = 30 + (40).*rand(820,1);
% ab = 0.5 + (2-0.5).*rand(90,1);
% ad = 5 + (30-5).*rand(90,1);
c(c < 0) = [];
ad(ad < 0) = [];
ab(ab < 0) = [];
c=c(1:820);
ad = ad(1:90*4);
ab=ab(1:90*4);
c=1./c;
ad=1./ad;
ab=1./ab;
% ca = cat(2,c,ad,ab);
ca = cat(2,c,ad,ab);
ca(ca > 300) = [];
histogram(ca)



%% 

function val = firing_rate(c,max)
    r = max/2;
    val = r*(1+tanh((c-r)/r));
end
