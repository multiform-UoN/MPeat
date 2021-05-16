function [ Kl,Cl,Ml ] = element( local_data )

  Kl = zeros(4,4);
  Cl = zeros(4,4);
  Ml = zeros(4,4);
  
%  Constitutive matrix - Hooke's Law and Darcy's Law
%
  Emod    = local_data.mater(1,1);  %Young's Modulus
  k       = local_data.mater(1,2);  %Hydraulic conductivity
  time    = local_data.mater(1,3);  %Time  
  sw      = local_data.mater(1,4);  %Degree of saturation 
  ss      = local_data.mater(1,5);  %Spesific storage or pore pressure coefficient
  t_final = local_data.mater(1,6);  %Maximum time simmulated
  
  constant_parameter = 0.45;
  area      = exp(-(time/(constant_parameter*t_final)));
  L         = abs(local_data.coords(2,1)-local_data.coords(1,1));
  
%-------------------------------- Evalaute constitutive matrices
C22 = (k*area/L)*[1  -1; -1  1];  
B11 = (sw*area*Emod/L)*[1  -1; -1  1];  
B22 = (ss*area*L)*[1/3  1/6; 1/6  1/3]; 
B12 = (area/L)*[L/2  L-(L/2); -L/2  -L+(L/2)];
B21 = B12';
%-------------------------------- Assemble local stiffness matrix 
    Kl = [ zeros(2,2)               zeros(2,2);
           zeros(2,2)                     C22 ];

%-------------------------------- Reorder stiffness matrix for degree of
%freedom

% ----- Reorder columns  
            tempcol = Kl(:,2); Kl(:,2) = Kl(:,3); Kl(:,3) = tempcol;
            
% ----- Reorder rows
            temprow = Kl(2,:); Kl(2,:) = Kl(3,:); Kl(3,:) = temprow;
            
%-------------------------------- Assemble local capacitance matrix 
    Cl = [ B11                              B12 ;
           B21                              B22  ];
           
%-------------------------------- Reorder capacitance matrix for dof

% ----- Reorder columns  
            tempcol = Cl(:,2); Cl(:,2) = Cl(:,3); Cl(:,3) = tempcol;
% ----- Reorder rows
            temprow = Cl(2,:); Cl(2,:) = Cl(3,:); Cl(3,:) = temprow;


