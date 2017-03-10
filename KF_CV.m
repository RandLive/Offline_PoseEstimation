PO = POM;
datain=PO';

% calculating the velocity of the moving catheter

Vx_state=zeros(1,size(datain,2));
Vy=zeros(1,size(datain,2));
Vz=zeros(1,size(datain,2));
Vphi=zeros(1,size(datain,2));
Vtheta=zeros(1,size(datain,2));

deltaT = 1;

for i=2:size(datain,2)
     Vx_state(i)=(datain(1,i)-datain(1,i-1))/deltaT;
     Vy(i)=(datain(2,i)-datain(2,i-1))/deltaT;
     Vz(i)=(datain(3,i)-datain(3,i-1))/deltaT;
     Vphi(i)=(datain(4,i)-datain(4,i-1))/deltaT;
     Vtheta(i)=(datain(5,i)-datain(5,i-1))/deltaT;
end

datain=[datain;Vx_state;Vy;Vz;Vphi;Vtheta];

M=10;                                                                                  % M is the observation chanel number
N=size(datain,2);                                                                % N is the detection times 
F=[eye(M/2),deltaT*eye(M/2);zeros(M/2),eye(M/2)];   % state transition matrix_state 
B=0;                                                                                      % control matrix_state
H=[eye(M/2),zeros(M/2)];                                                 % obsevation matrix_state/measurement matrix_state

%initialize the noise covariance matrix_state

P=100*eye(M);                                                                   %state variance matrix_state/error covariance matrix_state
Q=1e-12*[(deltaT^3)/3*eye(M/2),(deltaT^2)/2*eye(M/2);(deltaT^2)/2*eye(M/2),deltaT*eye(M/2)];     %Process covariance matrix_state
% R=cov(datain(1:3,:)');                                                     %observation covariance matrix_state
R=zeros(M/2);

for tt=1:M/2
    R(tt,tt)=var(datain(tt,:));
end

%x_state,x_state_,u,u_,z initialize
u=0; %u is the control value, u_ is the updated control value for the target

x_state=zeros(M,N);
x_state_=zeros(M,N);
x_state(:,1)=datain(:,1);

for i=2:N
     x_state_(:,i)=F*x_state(:,i-1)+B*u;
     P_=F*P*F'+Q;
     Kt=P_*H'/(H*P_*H'+R);
     x_state(:,i)=x_state_(:,i)+Kt*(datain(1:M/2,i)-H*x_state_(:,i));
     xx_state(:,i)=x_state(1:5,i);
     P=(eye(M)-Kt*H)*P_;
end

xx=x_state(1:M/2,N)';
yy=datain(1:M/2,N)';