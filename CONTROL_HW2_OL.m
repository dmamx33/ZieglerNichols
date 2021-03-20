close all
clear all
clc
%%%% TF CREATION
zeta1 = 1;
month = 12;
[num, den]=ord2(month,zeta1);
A1=tf(num, den);
  %%%%
zeta2 = 1.3;
wn2 = 23; 
[num, den]=ord2(wn2,zeta2);
A2=tf(num, den);
Ksys= 9600;
  %%%%
G = series(Ksys,series(A1,A2));
%G=tf([0 0 9.48],[1 10.666 90]);TEST TF
[num den] = tfdata(G,'v');
%%%% T VECTOR
TS= 0.00001;
Tfinal = 10;
t=0:TS:Tfinal;
%%%% RESPONSES
derivator = tf([1 0],[0 1]);
y = step(G,t);
dy = step(series(G,derivator),t);
%%%% IDENTIFY y(L) = mL +b
[MAX imax]=max(dy);
Tmax=t(imax);
s=sprintf('Maximum value of %f at %f sec',MAX,Tmax);
disp(s)
%%%% SYMBOLIC TF CREATION
syms s dens nums J(s) P(s) m LL b yy n f
nums=poly2sym(num,s);
dens=poly2sym(den,s);
J(s)=nums/dens;
pretty(J(s))
%%%% IDENTIFY K with Final value theorem
P(s)=s*J(s);
P(s)=subs(J(s),s,0);
K = double(P(s));
s=sprintf('Value of K=%f',K);
disp(s)
%%%% IDENTIFY b 
f = (m*LL) - b;
b=solve(f==yy,b);
b=subs(b,yy,sym(y(imax)));
b=subs(b,m,sym(MAX));
b=subs(b,LL,sym(Tmax));
%pretty(b);
%%%%% CREATION ZN LINE
znl= (MAX*t) - double(b);
%%%%% IDENTIFY L & T
Lind = find(abs(znl)==min(abs(znl)));
L=t(Lind)
aux = znl - K;
Tind=find(abs(aux)==min(abs(aux)));
T=t(Tind)-t(Lind)
%%%% PLOTTING
figure
hold on
plot(t,y,'LineWidth',2)
plot(t,dy,'LineWidth',2)
%plot(t,MAX*t-0.025)
plot(t,znl)
line([Tmax Tmax],[-10 10],'Color','black','LineStyle',':')%vertical
line([0 10],[K K],'Color','[0.4 0.4 0.4]','LineStyle','--')%horizontal
line([T T],[-10 10],'Color','black','LineStyle','--')
line([t(Lind) t(Lind)],[-10 10],'Color','black','LineStyle','--')
%yline(MAX)%For MATLAB R2018b
title('Ziegler Nichols Open Loop method')
xlabel('Time')
ylabel('Amplitude')
legend('Step response Y(s)','Derivate of Y(s)','Ziegler Nichols line')
axis([0 1.5 -0.1 0.5])
grid
hold off
%%%%ZN CONSTANTS
  %%%ONLY P CONTROLLER
    %P = T/L;
  %%%ONLY PI CONTROLLER
    %P=0.9*(T/L)
    %tau_i=L/0.3;
  %%ONLY PID CONTROLLER
     P=1.2*(T/L)
     tau_i=2*L
     tau_d=0.5*L
%%% PID DEF
PID = tf([P*tau_d*tau_i P*tau_i P],[0 tau_i 0])
%%%% CONTROL SYSTEM SIMULATION
ycomp=step(feedback(series(PID,G),1),t);
%%%% Modifications
     P2=P*4.2;
     tau_i2=tau_i*1.2;
     tau_d2=tau_d*2;
     PID2 = tf([P2*tau_d2*tau_i2 P2*tau_i2 P2],[0 tau_i2 0])
     ycomp2=step(feedback(series(PID2,G),1),t);
%%%% PLOTTING
figure
hold on
plot(t,ycomp,'LineWidth',2)
plot(t,ycomp2,'LineWidth',2)
plot(t,y,'LineWidth',2)
plot(t,ones(length(t),1),'--k')
legend('Output with ZN PID controller','Output with PID modified','Output Open loop','Input')
title('Ziegler Nichols Open Loop method')
xlabel('Time')
ylabel('Amplitude')
axis([0 2 0 1.1])
grid
hold off
