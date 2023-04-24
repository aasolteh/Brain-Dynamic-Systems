clc; clear;

%% The Morris-Lecar Model

% Initializations
mscm2   = 10;
mv      = 0.001;
muFcm2  = 0.01;

gL      = 8;
EL      = -80;
gNa     = 20;
ENa     = 60;
gK      = 10;
EK      = -90;
V2n     = -25;
V2m     = -20;
kn      = 5;
km      = 15;
Cm      = 1;
tau     = 1;
I       = 0;

%% Part I
% Start PPlane
% pplane8;

% Equations & Finding Equilbrium Points
syms V n

dV = (I - gL * (V - EL) - gK * n * (V-EK) - gNa * 1/(1+exp((-20-V)/km)) * (V-ENa))/Cm;
dn = 1/tau * (1/(1+exp((V2n-V)/kn)) - n);

morisLecarEqns = [dV == 0, dn == 0];
S    = vpasolve(morisLecarEqns, [V n]);
Veq  = double(S.V);
neq  = double(S.n);

% Plotting Nullclines
nullcline1 = @(v)((I - gL .* (v - EL) - gNa .* 1./(1+exp((-20-v)./km)) .* (v-ENa))./ (Cm .* (gK .* (v-EK))));
nullcline2 = @(v)1./(1+exp((V2n-v)./kn));
v          = linspace(-100, -20, 100);
plot(v,nullcline1(v),'b'); 
hold on
plot(v,nullcline2(v),'r');
title('Nullcline of differential equations')

% Type of Equilbrium Points
J       = jacobian([dV, dn], [V, n]);
eValues = eig(subs(J, [V, n], [Veq, neq]));
eValue1 = double(eValues(1));
eValue2 = double(eValues(2));

%% Part II

stepamplitudes = 0.25:0.2:10;
for i = 1:length(stepamplitudes)
    [t,y] = ode45(@(t,y) func(t,y,stepamplitudes(i)),[0 100],[-66; 0]);
    subplot(7, 7, i)
    plot(t, y(:,1))
end
threshold = stepamplitudes(23);
clear t y 

Ib = threshold - 1;
Ia = threshold + 1;

[tb,yb] = ode45(@(t,y) func(t,y,Ib),[0 100],[-66; 0]);
[ta,ya] = ode45(@(t,y) func(t,y,Ia),[0 100],[-66; 0]);

figure 
subplot(1, 2, 1)
plot(tb, yb(:,1))
xlabel('t')
ylabel('v')
title('Voltage Plot for current below threshold')
subplot(1, 2, 2)
plot(ta, ya(:,1))
xlabel('t')
ylabel('v')
title('Voltage Plot for current above threshold')

%% Part III
% Find Equilibrium Points
syms V n

dVb = (Ib - gL * (V - EL) - gK * n * (V-EK) - gNa * 1/(1+exp((-20-V)/km)) * (V-ENa))/Cm;
dnb = 1/tau * (1/(1+exp((V2n-V)/kn)) - n);

morisLecarEqnsb = [dVb == 0, dnb == 0];
Sb    = vpasolve(morisLecarEqnsb, [V n]);
Veqb  = double(Sb.V);
neqb  = double(Sb.n);

dVa = (Ia - gL * (V - EL) - gK * n * (V-EK) - gNa * 1/(1+exp((-20-V)/km)) * (V-ENa))/Cm;
dna = 1/tau * (1/(1+exp((V2n-V)/kn)) - n);

morisLecarEqnsa = [dVa == 0, dna == 0];
Sa    = vpasolve(morisLecarEqnsb, [V n]);
Veqa  = double(Sa.V);
neqa  = double(Sa.n);

% Type of Equilbrium Points
Jb       = jacobian([dVb, dnb], [V, n]);
eValuesb = eig(subs(Jb, [V, n], [Veqb, neqb]));
eValue1b = double(eValuesb(1));
eValue2b = double(eValuesb(2));

Ja       = jacobian([dVa, dna], [V, n]);
eValuesa = eig(subs(Ja, [V, n], [Veqa, neqa]));
eValue1a = double(eValuesa(1));
eValue2a = double(eValuesa(2));

% Plotting Nullclines

nullcline1a = @(v)((Ia - gL .* (v - EL) - gNa .* 1./(1+exp((-20-v)./km)) .* (v-ENa))./ (Cm .* (gK .* (v-EK))));
nullcline2a = @(v)1./(1+exp((V2n-v)./kn));

nullcline1b = @(v)((Ib - gL .* (v - EL) - gNa .* 1./(1+exp((-20-v)./km)) .* (v-ENa))./ (Cm .* (gK .* (v-EK))));
nullcline2b = @(v)1./(1+exp((V2n-v)./kn));

v          = linspace(-100, -20, 100);
subplot(1, 2, 1)
plot(v,nullcline1b(v),'b'); 
hold on
plot(v,nullcline2b(v),'r');
title('Nullcline of differential equations for below threshold')

subplot(1, 2, 2)
plot(v,nullcline1a(v),'b'); 
hold on
plot(v,nullcline2a(v),'r');
title('Nullcline of differential equations for above threshold')

%% Part IV
stepamplitudes = 0.01:1:100;
for i = 1:length(stepamplitudes)
    [t,y] = ode45(@(t,y) func2(t,y,stepamplitudes(i)),[0 100],[-66; 0]);
    subplot(10, 10, i)
    plot(t, y(:,1))
end

impulseThreshold = stepamplitudes(28);

clear t y 

Ib = impulseThreshold - 1;
Ia = impulseThreshold + 1;

[tb,yb] = ode45(@(t,y) func2(t,y,Ib),[0 100],[-66; 0]);
[ta,ya] = ode45(@(t,y) func2(t,y,Ia),[0 100],[-66; 0]);

figure 

subplot(1, 2, 1)
plot(tb, yb(:,1))
xlabel('t')
ylabel('v')
title('Voltage Plot for current below threshold for impulse input')
subplot(1, 2, 2)
plot(ta, ya(:,1))
xlabel('t')
ylabel('v')
title('Voltage Plot for current above threshold for impulse input')

%% LIF model

% Initializations
Vth         = -55 * 10^(-3);
Vreset      = -75 * 10^(-3);
tau         = 10 * 10^(-3);
gL          = 10 * 10^(-9);
Vinit       = -75 * 10^(-3);
EL          = -75 * 10^(-3);
tref        = 2 * 10^(-3);
T           = 400 * 10^(-3);
dt          = 0.1 * 10^(-3);

t           = 0:dt:T;
Vm          = nan(1, size(t, 2));
Vm(1)       = Vinit;

%% Part I
f = 10:20:1000;
oscillationAmplitude = nan(1, size(f, 2));

for j=1:length(f)
    I = sin(2*pi*t*f(j)) * 0.1 * 10^(-9);
    i = 1;
    while i < size(t, 2)
        if Vm(i) > Vth
            Vm(i+1:i+tref/dt) = Vreset;
            i = i + tref/dt;
        else
            Vm(i+1) = Vm(i) + dt * (-(Vm(i) - EL) + I(i) / gL) / tau;
            i = i + 1;
        end
    end
    oscillationAmplitude(j) = abs(max(Vm) - min(Vm)) / 2;
end

plot(f, oscillationAmplitude)
xlabel('f (Hz)')
ylabel('Amplitude of Oscillation')
title('Effect of frequency on amplitude of oscillation')

%% Part II
I = 0.1 * 10^(-9) * mystepfunc(t);

i = 1;
while i < size(t, 2)
    if Vm(i) > Vth
        Vm(i+1:i+tref/dt) = Vreset;
        i = i + tref/dt;
    else
        Vm(i+1) = Vm(i) + dt * (-(Vm(i) - EL) + I(i) / gL) / tau;
        i = i + 1;
    end
end

Fs = 1/T;
fft_i = fft(I);
fft_v = fft(Vm);
L = size(t, 2);

P2_i = abs(fft_i/L);
P1_i = P2_i(1:L/2+1);
P1_i(2:end-1) = 2*P1_i(2:end-1);

P2_v = abs(fft_v/L);
P1_v = P2_v(1:L/2+1);
P1_v(2:end-1) = 2*P1_v(2:end-1);

H = P1_v ./ P1_i;

f = Fs*(0:(L/2))/L;
plot(f,H) 
title('Frequency response of LIF neuron')
xlabel('f (Hz)')
ylabel('|H(f)|')


%% Part III
I = 0.15 * 10^(-9) * mystepfunc(t);

i = 1;
while i < size(t, 2)
    if Vm(i) > Vth
        Vm(i+1:i+tref/dt) = Vreset;
        i = i + tref/dt;
    else
        Vm(i+1) = Vm(i) + dt * (-(Vm(i) - EL) + I(i) / gL) / tau;
        i = i + 1;
    end
end

Vthreshold = Vth * mystepfunc(t);

plot(t, Vm)
hold on
plot(t, Vthreshold, '--')
title('Membrane potential vs threshold potential')
xlabel('t')
ylabel('v')
ylim([-0.075 -0.05])
legend('Membrane Potential', 'Threshold')

%% Part IV
tref = 0.003;
Iamps = 0.25:0.5:100;
numSpikes = nan(1, size(Iamps, 2));

for j = 1:length(Iamps)
    i = 1;
    I = Iamps(j) * 10^(-9) * mystepfunc(t);
    while i < size(t, 2)
        if Vm(i) > Vth
            Vm(i+1:i+tref/dt) = Vreset;
            i = i + tref/dt;
        else
            Vm(i+1) = Vm(i) + dt * (-(Vm(i) - EL) + I(i) / gL) / tau;
            i = i + 1;
        end
    end
    numSpikes(j) = sum(Vm == max(Vm)) / T;
end

plot(Iamps, numSpikes)
xlabel('I (nA)')
ylabel('Firing Rate')
title('Effect of current amplitude on firing rate')

%% The Hodgkin–Huxley Model
% Initializations 
V0  = -65;
n0  = 0.3;
m0  = 0.6;
h0  = 0.05;
ENa = 55;
EK  = -77;
EL  = -65;
gNa = 40;
gL  = 0.3;
gK  = 35;
C   = 1;
I   = 0;

% Equations Initialization
syms V n m h

dV = I - gK * n^4 * (V - EK) - gNa * m^3 * h * (V - ENa) - gL * (V - EL) / C;
dn = alphaN(V) * (1-n) - betaN(V) * n;
dm = alphaM(V) * (1-m) - betaM(V) * m;
dh = alphaH(V) * (1-h) - betaH(V) * h;

%% Part I
I1     = 20;
width1 = 0.2;
I2     = 0;
width2 = 0.2;
[t,y] = ode45(@(t,y) func3(t,y,I1,width1,I2,width2),[0 100],[V0; n0; m0; h0]);

plot(t, y(:, 1))
xlabel('time (ms)')
ylabel('Voltage (mv)')
title('Action potential over time')

Iamps = 50:50:1000;
figure
for j = 1:length(Iamps)
    [t,y] = ode45(@(t,y) func3(t,y,Iamps(j),width1,I2,width2),[0 100],[V0; n0; m0; h0]);
    subplot(4, 5, j)
    plot(t, y(:, 1))
    xlabel('time (ms)')
    ylabel('Voltage (mv)')
    title(['Action potential over time for I = ' , num2str(Iamps(j))])
end

width1 = 0.3;
figure
for j = 1:length(Iamps)
    [t,y] = ode45(@(t,y) func3(t,y,Iamps(j),width1,I2,width2),[0 100],[V0; n0; m0; h0]);
    subplot(4, 5, j)
    plot(t, y(:, 1))
    xlabel('time (ms)')
    ylabel('Voltage (mv)')
    title(['Action potential over time for I = ' , num2str(Iamps(j))])
end

width1 = 0.4;
figure
for j = 1:length(Iamps)
    [t,y] = ode45(@(t,y) func3(t,y,Iamps(j),width1,I2,width2),[0 100],[V0; n0; m0; h0]);
    subplot(4, 5, j)
    plot(t, y(:, 1))
    xlabel('time (ms)')
    ylabel('Voltage (mv)')
    title(['Action potential over time for I = ' , num2str(Iamps(j))])
end

width1 = 0.5;
figure
for j = 1:length(Iamps)
    [t,y] = ode45(@(t,y) func3(t,y,Iamps(j),width1,I2,width2),[0 100],[V0; n0; m0; h0]);
    subplot(4, 5, j)
    plot(t, y(:, 1))
    xlabel('time (ms)')
    ylabel('Voltage (mv)')
    title(['Action potential over time for I = ' , num2str(Iamps(j))])
end

width1 = 0.6;
figure
for j = 1:length(Iamps)
    [t,y] = ode45(@(t,y) func3(t,y,Iamps(j),width1,I2,width2),[0 100],[V0; n0; m0; h0]);
    subplot(4, 5, j)
    plot(t, y(:, 1))
    xlabel('time (ms)')
    ylabel('Voltage (mv)')
    title(['Action potential over time for I = ' , num2str(Iamps(j))])
end

width1 = 0.15;
figure
for j = 1:length(Iamps)
    [t,y] = ode45(@(t,y) func3(t,y,Iamps(j),width1,I2,width2),[0 100],[V0; n0; m0; h0]);
    subplot(4, 5, j)
    plot(t, y(:, 1))
    xlabel('time (ms)')
    ylabel('Voltage (mv)')
    title(['Action potential over time for I = ' , num2str(Iamps(j))])
end

Iexcitation = [1000, 750, 550, 400, 300, 250];
widths = [0.15, 0.2, 0.3, 0.4, 0.5, 0.6];
plot(widths, Iexcitation);
ylabel('I(mu A)');
title('Effect of excitation width on spiking current')

%% Part II
I1     = 20;
width1 = 0.2;
[t,y] = ode45(@(t,y) func3(t,y,I1,width1,I2,width2),[0 100],[V0; n0; m0; h0]);

I1     = 800;
[t1,y1] = ode45(@(t1,y1) func3(t1,y1,I1,width1,I2,width2),[0 100],[V0; n0; m0; h0]);

subplot(1, 2, 1)
plot(t, y(:, 2))
xlabel('time (ms)')
ylabel('n')
title('n over time under excitation')

subplot(1, 2, 2)
plot(t, y(:, 2))
xlabel('time (ms)')
ylabel('n')
title('n over time over excitation')

figure
subplot(1, 2, 1)
plot(t, y(:, 3))
xlabel('time (ms)')
ylabel('m')
title('m over time under excitation')

subplot(1, 2, 2)
plot(t1, y1(:, 3))
xlabel('time (ms)')
ylabel('m')
title('m over time over excitation')

figure
subplot(1, 2, 1)
plot(t, y(:, 4))
xlabel('time (ms)')
ylabel('h')
title('h over time under excitation')

subplot(1, 2, 2)
plot(t1, y1(:, 4))
xlabel('time (ms)')
ylabel('h')
title('h over time over excitation')

%% Part III
IK1  = gK * y(:, 2).^4 .* (y(:, 1) - EK);
IK2  = gK * y1(:, 2).^4 .* (y1(:, 1) - EK);
INa1 = gNa * y(:, 3).^3 .* y(:, 4) .* (y(:, 1) - ENa);
INa2 = gNa * y1(:, 3).^3 .* y1(:, 4) .* (y1(:, 1) - ENa);
IL1  = gL * (y(:, 1) - EL);
IL2  = gL * (y1(:, 1) - EL);

subplot(1, 2, 1)
plot(t, IK1)
xlabel('time (ms)')
ylabel('IK (mu A)')
title('IK over time under excitation')

subplot(1, 2, 2)
plot(t1, IK2)
xlabel('time (ms)')
ylabel('IK (mu A)')
title('IK over time over excitation')

figure
subplot(1, 2, 1)
plot(t, INa1)
xlabel('time (ms)')
ylabel('INa (mu A)')
title('INa over time under excitation')

subplot(1, 2, 2)
plot(t1, INa2)
xlabel('time (ms)')
ylabel('INa (mu A)')
title('INa over time over excitation')

figure
subplot(1, 2, 1)
plot(t, IL1)
xlabel('time (ms)')
ylabel('IL (mu A)')
title('IL over time under excitation')

subplot(1, 2, 2)
plot(t1, IL2)
xlabel('time (ms)')
ylabel('IL (mu A)')
title('IL over time over excitation')

%% The FitzHugh-Nagumo Model
pplane8;


%% Functions Initializations for Hodgkin-Huxley Model
function aM = alphaM(V)
aM = 0.182 * (V + 35) ./ (1 - exp(-(V + 35) / 9));
end

function bM = betaM(V)
bM = -0.124 * (V + 35) ./ (1 - exp((V + 35) / 9));
end

function aH = alphaH(V)
aH = 0.25 * exp(-(V + 90) / 12);
end

function bH = betaH(V)
bH = 0.25 * exp((V + 62) / 6) ./ exp((V + 90) / 12);
end

function aN = alphaN(V)
aN = 0.02 * (V - 25) ./ (1 - exp(-(V - 25) / 9));
end

function bN = betaN(V)
bN = -0.002 * (V - 25) ./ (1 - exp((V - 25) / 9));
end