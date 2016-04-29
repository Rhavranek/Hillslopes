%% Rachel Havranek, Feb 9, 2016
% Version 2 the Hillslopes model

%% Initialize
%horizontal array: this is analagous to my z array in the last fcts I did
dx = 2;
L=100;
xmax = L; % length L
xmin=-1*L;
x=xmin:dx:xmax;
N=length(x);

%topography
zbmax=100;
S0=0.2;
zb=zbmax-S0*abs(x);

%time array:
tmax=10e5; % years
dt=100; %time step, years
t=0:dt:tmax; %creates an array of time steps (time)

%Variables:
rhor=2750; %kg/m^3, typical granite
rhos=1300; %kg/m^3, from google search
wdot0=5e-5; 
kappa=0.003; %m^2/yr from: Roering et al., 2001
k=kappa*rhos;
edot=2e-3; %m/yr 
H0=1;
H=H0*ones(size(x));
Hstar=0.3; %m

z=zb+H;
 
imax=length(t); %loop goes every time step
nplots=100; %so we only see 100 plots, and not every time
tplot=tmax/nplots; %the amount of time between each plot

%% RUN

for i=1:imax
%weathering of bedrock:
wdot=wdot0*exp(-H/Hstar);

%change in the height of the bedrock
dzdx= diff(z)/dx; %make a dz/dx array - slope of the hill
 
 %Calculate the slope flux
Q=-k*dzdx; %the flux based on the slope of the hill

%rate of change of flux 
dQdx=diff(Q)/dx; %change in slope of the hill - second derivative (concave down)
 
dHdt=zeros(size(x));
dHdt(2:end-1)= ((rhor/rhos)*wdot(2:end-1))-((1/rhos)*dQdx);

H(2:end-1)=H(2:end-1)+(dHdt(2:end-1)*dt);
H=max(0,H);
zb(2:end-1)=zb(2:end-1)-(wdot(2:end-1)*dt);

H(1)=0;
H(end)=0;
zb(1)=zb(1)-(edot*dt);
zb(end)=zb(end)-(edot*dt);
%define H(1)/end and zb (1,end)

z=zb+H;
if(rem(t(i),tplot)==0)  
figure(1)
     plot(x,z,'k', 'linewidth', 2)
     hold on
     plot (x,zb,'r', 'linewidth',3)
     axis ([-L,L,70, 110])
      X=[x,fliplr(x)]; %this makes a vector from x, and flip x left to right
       Y=[zb,fliplr(z)];%makes a vector b/n zb and z
    fill(X,Y,'k'); % this colors  - fill the vector X through Y 
       Z=[zb,fliplr(ones(1,length(x)))];
       fill(X,Z,'m');
    %ht=text(95,-1000,['  ',num2str(t(i)/1000), ' kyears  '],'fontsize',18);
     xlabel('distance','fontname','arial','fontsize',21)
     ylabel('height','fontname','arial','fontsize',21)
     set(gca,'fontsize',18,'fontname','arial')
     pause(0.1)
     hold off
end
end


