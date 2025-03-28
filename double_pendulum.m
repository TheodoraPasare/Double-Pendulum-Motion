clear; close all; clc;
RT=1; % selector real time(1) / slow motion(0)
g=9.80665; % m/s^2; gravitational acceleration
% Physical parameters of  the system:
L1=2.2; L2=1.3; % m; lengths of the rods
m1=1.2; m2=1.7; % kg; masses of the bodies
% initial conditions - any angles between -180 and +180
theta10=100; theta20=-10; % degrees; initial angles
theta10=theta10*pi/180; theta20=theta20*pi/180; % conversion to radians
OM10=-30; OM20=20; % degrees/s; initial angular velocities
OM10=OM10*pi/180; OM20=OM20*pi/180; % conversion to rad/s
% Defining characteristic durations:
omega1=sqrt(g/L1); omega2=sqrt(g/L2); % pulsations of the components
T1=2*pi/omega1; T2=2*pi/omega2; % the periods of the components
T=max(T1,T2); % the characteristic time of the pendulum's motion
ti=0; tf=5*T; N=200000; t=linspace(ti,tf,N); dt=t(2)-t(1); % discreet time
% Starting values:
theta1=zeros(1,N); theta2=theta1; % angles
OM1=zeros(1,N); OM2=OM1; % angular speeds
theta1(1)=theta10; theta2(1)=theta20; % starting angles step 1
theta1(2)=theta10+OM10*dt; theta2(2)=theta20+OM20*dt; % starting angles step 1
OM1(1)=OM10; OM2(1)=OM20; % starting angular speeds
%Helpful notations:
miu=1+m1/m2; % dimensionless coefficient
r=L2/L1; % dimensionless coefficient
a11=miu; a22=r; % coefficients of the main diagonal (constants) 
for i=2:N-1 % reccurent cycle
    aux=theta2(i)-theta1(i);
    a21=cos(aux); a12=a21*r; % coefficients of the secondary diagonal (variables) 
    % For the current angular velocities, we use the backward difference derivatives:
    OM1(i)=(theta1(i)-theta1(i-1))/dt; % velocity of body 1 at step i 
    OM2(i)=(theta2(i)-theta2(i-1))/dt; % velocity of body 2 at step i
    b1=r*OM2(i)^2*sin(aux)-g/L1*miu*sin(theta1(i)); % 'free' term 1
    b2=-OM1(i)^2*sin(aux)-g/L1*sin(theta2(i)); % 'free' term 2
    A=[a11,a12;a21,a22]; B=[b1;b2]; % system matrix and the column of 'free' terms
    E=A\B; % solution of the linear system in matrix form 
    eps1=E(1); eps2=E(2); % current angular accelerations
    % Second order reccurences:
    theta1(i+1)=2*theta1(i)-theta1(i-1)+dt^2*eps1; % body 1
    theta2(i+1)=2*theta2(i)-theta2(i-1)+dt^2*eps2; % body 2
end;
OM1(N)=(theta1(N)-theta1(N-1))/dt; % angular velocity of body 1 at step N
OM2(N)=(theta2(N)-theta2(N-1))/dt; % angular velocity of body 2 at step N
toc; % displays calculation time for the numerical solutions
% Cartesian coordinates of the bodies:
x1=L1*sin(theta1); x2=x1+L2*sin(theta2); % horizontal coordinates
y1=-L1*cos(theta1); y2=y1-L2*cos(theta2); % vertical coordinates
% Kinetic energy, potential energy, total energy
T=1/2*(m1*L1^2*OM1.^2+m2*(L1^2*OM1.^2+L2^2*OM2.^2+2*L1*L2*OM1.*OM2.*cos(theta2-theta1)));
U=-g*((m1+m2)*L1*cos(theta1)+m2*L2*cos(theta2));
H=T+U; % the hamiltonian of the system - total energy
figure(1);
Lmax=L1+L2; % the graphic frame
coef=30; % controls the graphic dimensions of bodies
rg1=coef*m1^(1/3); rg2=coef*m2^(1/3); % radius of graphics
tic; simt=0; % starts the timer and initializes the simulation time
while simt<=tf % graphic cycle
  j=abs(t-simt)==min(abs(t-simt)); % search for the closest t in the discretization
  plot([0 x1(j) x2(j)],[0 y1(j) y2(j)],'-g','LineWidth',3); hold on; % draws the rods
  xlabel('x/m'); ylabel('y/m');
  plot(0,0,'.k','MarkerSize',10); % suspension joint
  plot(x1(j),y1(j),'.r','MarkerSize',rg1); % body 1
  plot(x1(j),y1(j),'.k','MarkerSize',10); % joint of body 1
  plot(x2(j),y2(j),'.b','MarkerSize',rg2); % body 2
  axis([-Lmax Lmax -Lmax Lmax]); axis square; % graphic frame
  text(3/5*Lmax,3/5*Lmax,['E = ',num2str(round(H(j))),' J']);
  if RT==1 % real time(1) / slow motion(0)
    simt=toc; % updates the simulation time with the system clock
    text(3/5*Lmax,4/5*Lmax,['t = ',num2str(round(t(j))),' s']);
  else
    simt=simt+1e-2; % increments time with a centisecond
    text(3/5*Lmax,4/5*Lmax,['t=',num2str(round(t(j)*100)),' cs']);
  end
  pause(1e-6); hold off
end
