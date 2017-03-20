global mu0 mur_sende mur_empfang R_sende N_sende Ia

%% Add noise
U_g = Ug;

%% Relevent Constants
mu0=4*pi*10^-7;			% Permeability of vacumm
mur_sende=1;			% Relative permeability of transmitter coils
mur_empfang=898.4;      % Relative permeability of the sensor coil
Gain_tr=10.9;           % The gain of the pre amplifier for the transmitter coil array
s=n_sample;           % Sample length

%% Cross-section area and gain for the sensor coil
% R_sende=0.02;
% N_sende=207;
% R_empfang=0.0002832;                                        % Radius of the sensor
% N_empfang=500;                                              % Turns of the sensor
A_Empfang=pi*R_empfang^2;                                   % Cross-section area of the transmitter coil
% Verstaerkung=567*3.295*3/4*6739.7/6000/Gain_tr;			% The gain of the amplifier for the sensor coil, 14.8 it the resistance of the circuit
Verstaerkung=1.573943202562500e+03/Gain_tr;

%% Calculation of current flow in the transmitter coils
t = linspace(0,1/100,s);
% Select sinewave?
if flag == 0
    I = U_g/(sqrt(14.8^2+(2*pi*f*1.165*10^-3)^2))*Gain_tr;			% 14.8 is the circuit resistance and 1.165e-3 is the inductance
    Ia = sqrt(2)*rms(I(1:1000));									% The amplitude of sinewave current
    Omega = 2*pi*f;												    % Angular speed for circuilar wave
    % Select ramp signal?
elseif flag == 1
    T = max(t);
    I = U_g/(sqrt(14.8^2+(2*pi*(1/max(t))*1.165*10^-3)^2))*Gain_tr;
    Ia = I(max(size(I)));										 	% The maximum value of the ramp current
    Omega = 1/T;													% Angular speed - general caculation: d(Phi)/d(t). Here, we have 1 rad over total time T in ramp. The gradient remains.
    % Select DC signal?
elseif flag == 2
    T = max(t);
    I = U_g/(sqrt(14.8^2+(2*pi*(1/max(t))*1.165*10^-3)^2))*Gain_tr;
    Ia = mean(I);													% The mean value of the DC current
    Omega = 1/T;
end

%% Orientations of each transmitter coil
oo=[o1;o2;o3;o4;o5;o6;o7;o8];

%% Distances between sensor and enmitter
r=[r1;r2;r3;r4;r5;r6;r7;r8];
rr1=rs-r1;rr2=rs-r2;rr3=rs-r3;rr4=rs-r4;
rr5=rs-r5;rr6=rs-r6;rr7=rs-r7;rr8=rs-r8;
rr=[rr1;rr2;rr3;rr4;rr5;rr6;rr7;rr8];


%% Magnetic dipole calculation
m=zeros(8,3,1000);
K = pi*N_sende*I*R_sende^2;
for i=1:1000
    thetae=oo(1,3);phie=oo(1,2);
    m(1,:,i)=[sind(thetae)*cosd(phie) sind(thetae)*sind(phie) cosd(thetae)]*K(i);    
    thetae=oo(2,3);phie=oo(2,2);
    m(2,:,i)=[sind(thetae)*cosd(phie) sind(thetae)*sind(phie) cosd(thetae)]*K(i);    
    thetae=oo(3,3);phie=oo(3,2);
    m(3,:,i)=[sind(thetae)*cosd(phie) sind(thetae)*sind(phie) cosd(thetae)]*K(i);    
    thetae=oo(4,3);phie=oo(4,2);
    m(4,:,i)=[sind(thetae)*cosd(phie) sind(thetae)*sind(phie) cosd(thetae)]*K(i);    
    thetae=oo(5,3);phie=oo(5,2);
    m(5,:,i)=[sind(thetae)*cosd(phie) sind(thetae)*sind(phie) cosd(thetae)]*K(i);    
    thetae=oo(6,3);phie=oo(6,2);
    m(6,:,i)=[sind(thetae)*cosd(phie) sind(thetae)*sind(phie) cosd(thetae)]*K(i);    
    thetae=oo(7,3);phie=oo(7,2);
    m(7,:,i)=[sind(thetae)*cosd(phie) sind(thetae)*sind(phie) cosd(thetae)]*K(i);    
    thetae=oo(8,3);phie=oo(8,2);
    m(8,:,i)=[sind(thetae)*cosd(phie) sind(thetae)*sind(phie) cosd(thetae)]*K(i);
end


%% Constant values for voltage (induced in sensor coil) calculation
C1 = mu0/(4*pi)*mur_sende;                                                   % Constant of B-field
C2 = Omega;																   % Value of angular frequency
% C2=2*pi*f;
C3 = -mur_empfang*N_empfang*A_Empfang;                                       % Constant for calculating the voltage

%% B-filed Calculation
B = zeros(s,8,3);
U2 = zeros(s,8);
U1 = zeros(1,8);
U1_nonoise = zeros(1,8);

if flag == 0
    PHI = zeros(s-1,8);
    dPHI = zeros(s-2,8);
    dt = zeros(8,s-2);
    V = zeros(s-2,8);
    V_Visual = zeros(s-2,8);
else
    PHI = zeros(s,8);
    dPHI = zeros(s-1,8);
    dt = zeros(8,s-1);
    V = zeros(s-1,8);
    V_Visual = zeros(s-1,8);
end

C = C1*C2*C3*Verstaerkung;
Ca = -mur_sende*mur_empfang*N_empfang*A_Empfang*Verstaerkung;
n_sensor = [sind(os(3))*cosd(os(2)); sind(os(3))*sind(os(2)); cosd(os(3))];

for j=1:1000
    for i=1:8
        B(j,i,:)=Bfeld_Dipole(m(i,:,j),rr(i,:));                           % This B1 is the calculated Magnetic flux density alone X,Y,Z axis at the position where sensor coil placed
    end
end

if flag == 0
    B = [B(250*100/f+1: 1000,:,:); B(1:250*100/f,:,:)];
end
    
Bx=B(:,:,1);By=B(:,:,2);Bz=B(:,:,3);

for j_1=1:8
    PHI(:,j_1)=Ca*(n_sensor(1).*Bx(:,j_1)+n_sensor(2).*By(:,j_1)+n_sensor(3).*Bz(:,j_1));
    dPHI(:,j_1)=diff(PHI(:,j_1));
    dt=diff(t)';
    V(:,j_1)=dPHI(:,j_1)./dt(j_1);
    V_Visual(:,j_1)=V(:,j_1)+wgn(max(size(V)),1,nl,'dBW','real');
end

U_norm = zeros(1,8);

%% Calculation of the Voltage Induced in the sensor coil.
for k=1:8
    
    U_noise = wgn(1,1,nl,'dBW','real');
    
    if flag == 0                                                           % Sine Wave
        if mean(V(20:30,k))<0
            U1(k) = mean(abs(V(25:924,k)))*(pi/2);
        else
            U1(k) = -mean(abs(V(25:924,k)))*(pi/2); 
        end
              
    elseif flag == 1                                                       % Ramp Wave, take signals between 9 and 10 ms.
        U1(k) = mean(V(900:999,k));
    else                                                                   % DC waveform
        U1(k) = mean(V(1:999,k));
    end
        
    U_norm(k) = U1(k)+ U_noise;
   
    
end