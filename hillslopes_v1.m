%%Rachel Havranek
%Hillslopes code

%% Initialize

%horizontal array: this is analagous to my z array in the last fcts I did
dx = 20;
L=1000;
xmax = L; % length L
xmin=-1*L;
x=xmin:dx:xmax;
N=length(x);

%time array:
tmax=100000; % years
dt=1000; %time step, years
t=0:dt:tmax; %creates an array of time steps (time)
%Variables:
rhor=2750; %kg/m^3, typical granite
rhos=2650; %kg/m^3, from google search
wdot=1; %m/1,000,000 years
wdot=wdot/1000; %m/1000 years
kappa=0.003; %m^2/yr from: Roering et al., 2001
k=0.002;
edot=4; %mm/yr -- this is from Whipple, Beyond Bedrock (Bob's 2nd author)
edot=edot/1000; %m/yr
H=2000; %m, max height of the hill 

%% Create an initialization space
Q=zeros(N,1); %an array of empty fluxes which we'll imprint onto
z=zeros(N,1); %an array of empty heights which we'll imprint onto
dzdx=((-rhor*wdot)/k)*x; 
zmax=H; 

imax=length(t); %loop goes every time step
nplots=100; %so we only see 100 plots, and not every time
tplot=tmax/nplots; %the amount of time between each plot

%% equations from our notes
%erosion - z(x) - dzdx - Q - dQdx - dHdt - H
%dHdt=(rhor/rhos)*wdot+kappa*second derivative of dz/dx

%% RUN
for i=1:imax
 
 z=(H-edot*t)+((rhor*wdot)/(2*k))*(L.^2-x.^2); % This equation takes H (the max height of the hill), subtracts the amount of erosion 
 
 %a problem of my code is that it errodes at the same amount
 %everywhere.. not just at the edges...
 
 dzdx(1:N-1)= diff(z)/dx; %make a dz/dx array - slope of the hill
 
 %Calculate the slope flux
 Q=-k*dzdx; %the flux based on the slope of the hill 

%rate of change of flux 
 dQdx=diff(Q)/dx; %change in slope of the hill - second derivative (concave down)
  %dQdx=wdot*rhor; %for steady state only

%change in max height:
dHdt=(rhor/rhos)*wdot-(1/rhos)*dQdx;

%update H:
H(2:N-1)= H(2:N-1)-(rhor/rhos)*wdot-(1/rhos)*dQdx;
end 

if(rem(t(i),tplot)==0)  
figure(1)
     plot(x,z)
     hold on
     xlabel('distance','fontname','arial','fontsize',21)
     ylabel('height','fontname','arial','fontsize',21)
     set(gca,'fontsize',18,'fontname','arial')
     pause(0.1)
     hold off
end

