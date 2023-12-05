B_B = 2.50;
D_B = 0.40;
LOA = 8.4;
rho_W = 997;
H_B = 1.7;
g = 9.81;
M_d = 150;
M_e = 10;
r_d = 0.230;
r_e = 0.090;
e_d = 0.150;
l_e = 0.500;
A1 = 1.5;
c1 = 0.428;
c2 = 0.006;
c3 = 0.149;
w = l_e / sqrt(B_B / (2 * g));
C_g = 5000;
B_g = 2500;
Omega = (8500*pi/30);

%primeira equação

k1 = -12*sqrt(2)*g*B_B*sqrt(B_B/g)*c2/(LOA*B_B^2 + 12*c3*B_B^2 + LOA*H_B^2);

k2 = -6*Omega*M_d*r_d^2/(B_B*D_B*(LOA*B_B^2+12*c3*B_B^2+LOA*H_B^2)*rho_W);

k3 = 6*Omega*M_d*r_d^2/(M_d*e_d^2 - M_e*e_d^2 + M_e*l_e^2 + 3*M_d*r_d^2 + 3*M_e*r_e^2);

k4 = -12*B_g/(M_d*e_d^2 - M_e*e_d^2 + M_e*l_e^2 + 3*M_d*r_d^2 + 3*M_e*r_e^2);

k5 = (12*g*LOA*pi*B_B*rho_W*D_B^2 - 2*g*LOA*rho_W*pi*B_B^3)/(2*pi*B_B*D_B*(LOA*B_B^2+12*c3*B_B^2+LOA*H_B^2)*rho_W);

k6 = -12*C_g/(M_d*e_d^2- M_e*e_d^2 + l_e^2*M_e + 3*M_d*r_d^2 + 3*M_e*r_e^2);



 k7 = (A1*LOA*w*rho_W*B_B^3 - 6*A1*LOA*w*B_B*rho_W*D_B^2)/(2*pi*B_B*D_B*(LOA*B_B^2+12*c3*B_B^2+LOA*H_B^2)*rho_W);
% 
% k2 = 6*A1*LOA*w*B_B*D_B^2*rho_W;
% 
% k3 = 2*g*LOA*pi*B_B^3*rho_W;
% 
% k4 = 12*g*LOA*pi*B_B*D_B^2*rho_W;
% 
% k5 = 2*pi*B_B*D_B*(LOA*B_B^2+12*B_B^3*c3+LOA*H_B^2)*rho_W; %primeiro que aparece dividindo
% 
% k6 = -6*Omega*M_d*r_d^2/(B_B*D_B*(LOA*B_B^2+12*B_B^2*c3+LOA*H_B^2)*rho_W);

 % segundo termo que aparece dividindo

% k8 = -12*sqrt(2)*g*B_B*sqrt(B_B/g)*c2/(LOA*B_B^2 + 12*B_B^2 + 12*B_B^2*c3 + LOA*H_B^2);



K1 = [k1,k2;
      k3,k4
    ];

K2 = [k5,0;
      0,k6

    ];

K3 = [k7;
      0;
    ];



X0 = [0;
      0
];

%segunda equação

% k10 = 12*C_g/(e_d^2*M_d - e_d^2*M_e + l_e^2*M_e + 3*M_d*r_d^2 + 3*M_e*r_e^2);
%  
% k11 = 6*Omega*M_d*r_d^2/(e_d^2*M_d - e_d^2*M_e + l_e^2*M_e + 3*M_d*r_d^2 + 3*M_e*r_e^2); % termo comum da divisão
%  
% k12 = 12*B_g/(e_d^2*M_d - e_d^2*M_e + l_e^2*M_e + 3*M_d*r_d^2 + 3*M_e*r_e^2);
% 
% k13 = 6*Omega*M_d*r_d^2;




