clear;clc;

global Anz_Punkte Anz_Spulen ...
    m C r funktionscount zeit U  ...
    max_abstand ...
    StartX StartY StartZ StartPhi StartTheta

N_Channel=8;
f=1e3;

cd DataSave
% % % -80dBW
load Um_NoNoise
load POR
Um=Um_NoNoise;
POR=POR;

% % -40dBW
% load Um_Noise_40dBW
% load POR_40dBW
% Um=Um_Noise_40dBW;
% POR=POR_40dBW;

cd ..

StartX = 0.1;
StartY = 0.1;
StartZ = 0.1;
StartPhi = 0;
StartTheta = 0;

POR=POR; %#ok<*ASGSL>
U1=Um;
omega=2*pi*f;


Anz_Punkte=1;                         % Anzahl der Testpunkte ?Number of testing point?
Anz_Spulen=8;                           % Anzahl der Spulen    (Number of coils)                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% required

mu0=4*pi*10^-7;                         % ?                   (Permeability)
mur_sende=1;                            % relative Permeabilität Sendespule  (Permeability of emitting coils)
mur_empfang=898.4;                      % relative Permeabilität Empfangsspule  (Permeability of sensing coils)

N_sende=207;                            % Windungszahl Sendespule  (number of windings - emitting coil)
N_empfang=500;                          % Windungszahl Empfangsspule (number of windings - sensing coil)
R_sende=0.015;                           % Radius Sendespule (Radius of emitting coils)
R_empfang=283.2*10^-6;                  % Radius Empfangsspule (Radius of sensing coils)
I=0.6602;
% I = 0.660156294774090;


% Frequenz Sendespule (frequency of signal on sensing coil)

A_Empfang=pi*R_empfang^2;               % Normalenfläche Empfangsspule (normal area of the sensing coil) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% direct value assign

C1=mu0/(4*pi)*mur_sende;                % Konstande für B-Feld (constant of B field)
C2=2*pi*f;                              % Konstante für Ableitung nach t (constant for w)
C3=-mur_empfang*N_empfang*A_Empfang;      % Konstante für Spannungsberechnung (constant for voltage calculation)
Verstaerkung= 144.3985;
C=Verstaerkung*C1*C2*C3;                % (Constant C)       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% direct value assign
K=pi*N_sende*I*R_sende^2;               % Vorfaktor magnetisches Dipol (factor of magnetic dipol)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% direct value assign
funktionscount=2000;                    % Abbruchbedingung für fsolve (stop option for fsolve functon)
zeit=0.5;                               % Zeitbegrenzung für fsolve (time limitation for fsolve)
max_abstand=2;
m=zeros(Anz_Spulen,3);
% r=zeros(Anz_Spulen,3);
faktor= [1 1 1 1 1 1 1 1]';            % (factor for calibration calculation)


%% 1P1O
% oo_ = [1 0 0;...
%       1 0 0;...
%       1 0 0;...
%       1 0 0;...
%       1 0 0;...
%       1 0 0;...
%       1 0 0;...
%       1 0 0];

%% 1PMO
oo_ = [1 0 0;...
    1 20 30;...
    1 40 60;...
    1 60 90;...
    1 80 120;...
    1 100 150;...
    1 120 180;...
    1 140 210];

% VOI
%
% oo_ = [1 0 0;...
%       1 0 0;...
%       1 0 0;...
%       1 0 0;...
%       1 0 0;...
%       1 0 0;...
%       1 0 0;...
%       1 0 0];

%% with or without optimization
% oo = oo_+coord_shift(:,1:3);
oo=oo_;

% magnetisches Dipol
phie=oo(1,2);    thetae=oo(1,3);        m(1,:)=faktor(1)*K*[sind(thetae)*cosd(phie)      sind(thetae)*sind(phie)        cosd(thetae)];
phie=oo(2,2);    thetae=oo(2,3);        m(2,:)=faktor(2)*K*[sind(thetae)*cosd(phie)      sind(thetae)*sind(phie)        cosd(thetae)];
phie=oo(3,2);    thetae=oo(3,3);    	m(3,:)=faktor(3)*K*[sind(thetae)*cosd(phie)      sind(thetae)*sind(phie)        cosd(thetae)];
phie=oo(4,2);    thetae=oo(4,3);   	    m(4,:)=faktor(4)*K*[sind(thetae)*cosd(phie)      sind(thetae)*sind(phie)        cosd(thetae)];
phie=oo(5,2);    thetae=oo(5,3);        m(5,:)=faktor(5)*K*[sind(thetae)*cosd(phie)      sind(thetae)*sind(phie)        cosd(thetae)];
phie=oo(6,2);	 thetae=oo(6,3);    	m(6,:)=faktor(6)*K*[sind(thetae)*cosd(phie)      sind(thetae)*sind(phie)        cosd(thetae)];
phie=oo(7,2);    thetae=oo(7,3);     	m(7,:)=faktor(7)*K*[sind(thetae)*cosd(phie)      sind(thetae)*sind(phie)        cosd(thetae)];
phie=oo(8,2);    thetae=oo(8,3);        m(8,:)=faktor(8)*K*[sind(thetae)*cosd(phie)      sind(thetae)*sind(phie)        cosd(thetae)];

% % Ortskoordinaten der Spulen (position coordinate of the emitting coils)
%% 1P1O, 1PMO
r_ = [-0.06 0.05 0;...
    -0.02 0.05 0;...
    0.02 0.05 0;...
    0.06 0.05 0;...
    -0.06 -0.05 0;...
    -0.02 -0.05 0;...
    0.02 -0.05 0;...
    0.06 -0.05 0];

%% VOI
% r_ =[-0.28  -0.28 0.58;...
%     -0.28  -0.28 0.02;...
%     0.28   0.28 0.58;...
%     0.28   0.28 0.02;...
%     -0.28   0.28 0.58;...
%     -0.28   0.28 0.02;...
%     0.28  -0.28 0.58;...
%     0.28  -0.28 0.02];

% r_ =[-0.3  -0.3 0.6;...
%     -0.3  -0.3 0.0;...
%     0.3   0.3 0.6;...
%     0.3   0.3 0.06;...
%     -0.3   0.3 0.6;...
%     -0.3   0.3 0.0;...
%     0.3  -0.3 0.6;...
%     0.3  -0.3 0.0];

%% with or without optimization
% r=r_+coord_shift(:,4:6);
r=r_;



length_U = max(size(U1));
POM = zeros(length_U,5);
POM_KF = zeros(length_U,5);

mode = input(' 1-no update \n 2-update with previous value \n 3-update with filtered value \n 4-update with ground truth \n Enter the working mode:');

StartX = 0;
StartY = 0;
StartZ = 0.1;
StartPhi = 45;
StartTheta = 45;

for i=1:length_U
    
    switch mode
        case 1
            StartX = 0;
            StartY = 0;
            StartZ = 0.1;
            StartPhi = 0;
            StartTheta = 0;
            U=U1(i,:);
            [x,y,z,phi,theta]=POS(U,StartX,StartY,StartZ,StartPhi,StartTheta);
            POM(i,:)=[x,y,z,phi,theta];
            disp(i);
            pause(0);
            if i>=15
                POM_KF(i,:) = fun_KF_CV(POM(i-14:i,:))';
            else
                POM_KF(i,:) = POM(i,:);
            end
        case 2
            [x,y,z,phi,theta]=POS(U,StartX,StartY,StartZ,StartPhi,StartTheta);
            POM(i,:)=[x,y,z,phi,theta];
            disp(i);
            pause(0);
            StartX = x;
            StartY = y;
            StartZ = z;
            StartPhi = phi;
            StartTheta = theta;
            U=U1(i,:);
            if i>=15
                POM_KF(i,:) = fun_KF_CV(POM(i-14:i,:))';
            else
                POM_KF(i,:) = POM(i,:);
            end
            
        case 3
            U=U1(i,:);
            [x,y,z,phi,theta]=POS(U,StartX,StartY,StartZ,StartPhi,StartTheta);
            POM(i,:)=[x,y,z,phi,theta];
            if i>=15
                POM_KF(i,:) = fun_KF_CV(POM(i-14:i,:))';
            else
                POM_KF(i,:) = POM(i,:);
            end
            disp(i);
            pause(0);
            StartX = POM_KF(i,1);
            StartY = POM_KF(i,2);
            StartZ = POM_KF(i,3);
            StartPhi = POM_KF(i,4);
            StartTheta = POM_KF(i,5);
            if i>=15
                POM_KF(i,:) = fun_KF_CV(POM(i-14:i,:))';
            else
                POM_KF(i,:) = POM(i,:);
            end
            
        case 4
            StartX = POR(i,1);
            StartY = POR(i,2);
            StartZ = POR(i,3);
            StartPhi = POR(i,4);
            StartTheta = POR(i,5);
            U=U1(i,:);
            [x,y,z,phi,theta]=POS(U,StartX,StartY,StartZ,StartPhi,StartTheta);
            POM(i,:)=[x,y,z,phi,theta];
            disp(i);
            pause(0);
            if i>=15
                POM_KF(i,:) = fun_KF_CV(POM(i-14:i,:))';
            else
                POM_KF(i,:) = POM(i,:);
            end
    end
    
end

POM = [1000*POM(:,1:3), POM(:,4:5)];
POM_KF = [1000*POM_KF(:,1:3), POM_KF(:,4:5)];
POR = [1000*POR(:,1:3), POR(:,4:5)];

scatter3(POM(:,1),POM(:,2),POM(:,3),'bo');hold on;scatter3(POR(:,1),POR(:,2),POR(:,3),'r.');legend('Estimated Pose','Ground Truth');axis([0 120 0 250 350 700]);
set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',20);xlabel('X Position (mm)');ylabel('Y Position (mm)');zlabel('Z Position (mm)');
