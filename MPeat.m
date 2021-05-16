%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% MPeat Program 1D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
%MPeat is the long-term peatland development model that includes mechanical, 
%ecological, and hydrological processes via poroelasticity theory, which coupling 
%between fluid flow and solid deformation. 

%First version 2021

%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% Input Parameter Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
acrotelm_decay = 5e-2; %Acrotelm decay rate(yr^{-1})
catotelm_decay = 8e-5; %Catotelm decay rate(yr^{-1})
hydraulic_param = 15; %Hydraulic conductivity parameter (-)
modulus_param_1 = 2e5; %Young’s modulus parameter 1 (Pa) 
modulus_param_2 = 0.1; %Young’s modulus parameter 2 (-)
time_max = 6000; %Maximum time simmulated (year)
annual_timestep = 10; %Annual timestep for ecological and hydrological submodel
timestep = 1/annual_timestep;
annual_timestep_mechanic = 10; %Annual timestep for mechanical submodel
mechanic_timestep = 1/annual_timestep_mechanic; 
gravity = 9.8; %Gravitational acceleration (m s^{-2})
hydraulic_initial = 1e-2/3.171e-8; %Hydraulic conductivity initial value (m yr{-1})
sw_un = 0.4; %Degree of saturation of water in the acrotelm (-)
sw_sat = 1; %Degree of saturation of water in the catotelm (-)
weight_water = 9800; %Specific weight of water (N m^{-3})
water_retention_1 = 0.5; %Water retention empirical constant 1 (-)
water_retention_2 = 0.4; %Water retention empirical constant 2 (-)
mechanic_time_average = 50; %Average mechanical time (yr)
turning_point = 0; %Turning point for acrotelm and catotelm transition (-)
shrub_proportion = 0.61; %Shrub proportion (-)
sedge_proportion = 0.09; %Sedge or herb proportion (-)
spagnum_proportion = 0.3; %Sphagnum proportion (-)
shrub_wet_constant = 0.4; %Shrub wet condition constant (-)
sedge_wet_constant = 0.4; %Sedge or herb wet condition constant (-)
spagnum_wet_constant = 20; %Sphagnum wet condition constant (-)
limit_density = 175; %Maximum bulk density (kg m^{-3})
limit_porosity = 0.15; %Minimum active porosity (-)
annual_time = 1;
layer_index = 0;
time_counter = 0;
wt_sum = 0;
carbon_content = 0.4; %Carbon content (-)
peatland_radius = 500; %Peatland radius (m)

%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% Initial Array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
layer_initial_mass = zeros(1,time_max);
layer_mass = zeros(1,time_max);
layer_remaining_mass = zeros(1,time_max);
layer_thickness = zeros(1,time_max);
layer_elevation = zeros(1,time_max);
wet_proportion = zeros(1,time_max);
hydraulic = zeros(1,time_max);
layer_transmissivity = zeros(1,time_max);
peatland_height = zeros(1,time_max);
water_table_height = zeros(1,time_max);
water_table_depth = zeros(1,time_max);
density = zeros(1,time_max);
porosity = zeros(1,time_max);
modulus = zeros(1,time_max);

%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% Initial Value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
layer_mass(1) = (9.3^2)*0.001; % Peat mass (kg)
density(1) = 50; %Initial bulk density (kg m^{-3})
layer_initial_mass(1) =layer_mass(1); %Initial mass (kg)
layer_remaining_mass(1) = 1; %Proportion of remaining mass (-)
wet_proportion(1) = 1; %Wet proportion in the acrotelm or catotelm (-)
layer_thickness(1) = layer_mass(1)/density(1); %Layer thickness (m)
layer_elevation(1) = layer_thickness(1); %Layer elevation (m)
hydraulic(1) = hydraulic_initial; %Hydraulic conductivity (m yr^{-1})
layer_transmissivity(1) = layer_thickness(1)*hydraulic(1); %Layer transmisivity (m^2 yr{-1})
bog_height = layer_elevation(1); %Peatland height (m)
wt_height = bog_height; %Water table height (m)
bog_transmissivity = layer_transmissivity(1);
displacement(1) = 0; %Vertical displcement (m)
porosity(1) = 0.8; %Active porosity initial value (-)


%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% Main Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
for ma = 1 : time_max/mechanic_time_average;
    mechanic_time(ma) = ma*mechanic_time_average;
end

while (1);
    
while (1);
time_counter = time_counter + 1;
if(time_counter > annual_timestep);
  break;
end;
bog_height = 0;
turning_point = 0;

%Calculate the changes of layer mass  with constant decay rate for each zone (acrotelm or catotelm). 
for layer_index = 1: annual_time;
    layer_mass (layer_index) = (layer_mass (layer_index) * wet_proportion (layer_index)*(exp (- catotelm_decay * timestep))) + (layer_mass (layer_index) * (1 - wet_proportion (layer_index)) * (exp (- acrotelm_decay * timestep)));
    layer_remaining_mass(layer_index) = layer_mass(layer_index) ./ layer_initial_mass(layer_index);
%Calculate the changes of bulk density.
    density(layer_index) = density(layer_index)*((layer_elevation(layer_index)/(layer_elevation(layer_index)-abs(displacement(layer_index))*(1+water_table_depth(layer_index))))); 
    if  density(layer_index) > limit_density;
        density(layer_index) = limit_density;
    else
        density(layer_index) = density(layer_index);
    end
    layer_thickness(layer_index) = layer_mass(layer_index) ./ density(layer_index);
%Calculate the changes of  Young's modulus. 
    modulus(layer_index) = modulus_param_1*(1+(layer_remaining_mass(layer_index)^modulus_param_2));
%Calculate the changes of  active porosity.  
    porosity(layer_index) = porosity(layer_index)*(((layer_elevation(layer_index)-abs(displacement(layer_index))*(1+water_table_depth(layer_index)))/layer_elevation(layer_index)));
    if  porosity(layer_index) < limit_porosity;
        porosity(layer_index) = limit_porosity;
    else
        porosity(layer_index) = porosity(layer_index);
    end
    if  layer_index == 1;
        layer_elevation(layer_index) = layer_thickness(layer_index);
    else
        layer_elevation(layer_index) = layer_thickness(layer_index) + layer_elevation(layer_index-1);
    end  
end;

bog_height = layer_elevation(annual_time);
peatland_height(annual_time) = bog_height;  

%Calculate the changes of water table height.
wt_height = wt_height +(timestep * net_rainfall(annual_time) / porosity(annual_time))-((timestep * 2 * bog_transmissivity * wt_height/(peatland_radius ^ 2)) / porosity(annual_time));
bog_transmissivity = 0.0;
    if      wt_height > bog_height ;
            wt_height = bog_height; 
    end;
water_table_height(annual_time) = wt_height;
%Calculate the changes of water table depth.
wt_depth = bog_height - wt_height;
water_table_depth(annual_time) = wt_depth;
wt_sum = wt_sum + wt_depth;

%Calculate the wet proportion.
for layer_index = 1: annual_time;
    if  layer_index == 1;
        if  wt_height >= layer_elevation(layer_index);
            wet_proportion(1) = 1;
        else
            wet_proportion(1) = wt_height ./ layer_thickness(1);
        end;
    elseif  wt_height >= layer_elevation(layer_index);
            wet_proportion(layer_index) = 1;
    elseif  wt_height <= layer_elevation(layer_index - 1) ;
            wet_proportion(layer_index) = 0;
    else;
            wet_proportion(layer_index) = (wt_height - layer_elevation(layer_index - 1))./ layer_thickness(layer_index);
    end;
%Calculate the changes of  hydraulic conductivity and transmisivity. 
    hydraulic(layer_index) = hydraulic_initial*((porosity(layer_index)/porosity(1))^hydraulic_param);
    layer_transmissivity(layer_index) = hydraulic(layer_index) * layer_thickness(layer_index)* wet_proportion(layer_index);
    bog_transmissivity = bog_transmissivity + layer_transmissivity(layer_index);
end; 

for lb = 1 : annual_time;
    if  layer_elevation(lb) <= water_table_height(annual_time);
        turning_point = turning_point+1;    
    end
end

%Calculate the layer height for peat.
for lc = 1: annual_time;
    if  lc == 1
        layer_height(lc) = 0;
    elseif layer_elevation(lc) <= water_table_height(annual_time);    
        layer_height(lc) = (lc)*(water_table_height(annual_time)/(turning_point));
    else
        layer_height(lc) = layer_elevation(lc);
    end
end
%Calculate the layer depth for peat.
for ld = 1 : annual_time;
    layer_depth(ld) = layer_height(annual_time)-layer_height(ld);
end
end;

time_counter = 0;
annual_time = annual_time + 1;

if  annual_time > time_max;
  break;
end;

wt_depth_av = wt_sum ./(annual_timestep);
wt_sum = 0.0;
%Calculate the new layer addition based on water table position and
%temperature.
if  wt_depth_av > 0.668
    layer_mass(annual_time) = 0.0000001;
else;
    layer_mass(annual_time) =(((9.3 +(133 * wt_depth_av) -(0.022*((wt_depth_av*100)^2)))^2)*0.001)*(0.1575*temp(annual_time)+0.0091);
end;

%New layer initial properties
porosity(annual_time) = porosity(1);
density(annual_time) = density(1);
layer_initial_mass(annual_time) = layer_mass(annual_time);
layer_remaining_mass(annual_time) = 1;
wet_proportion(annual_time) = 0;
layer_thickness(annual_time) = layer_mass(annual_time) / density(annual_time);
layer_elevation(annual_time) = bog_height + layer_thickness(annual_time);

%Load sources from new layer addition and plant weight
layer_addition_weight(annual_time) = layer_initial_mass(annual_time)*gravity;
plant_weight(annual_time) = (((((10^(((log10(layer_initial_mass(annual_time)))+0.409)/0.985))*shrub_proportion*(1+shrub_wet_constant))+((10^(((log10(layer_initial_mass(annual_time)))+0.001)/1))*sedge_proportion*(1+sedge_wet_constant))))+((0.144*spagnum_proportion)*(1+spagnum_wet_constant)))*gravity;

%Mechanical submodel main calculation 
mb = ismember(annual_time,mechanic_time);
if mb == 1 ; 
%Pore pressure coefficient
for ab = 1:annual_time-1
%Saturated specific storage or pore pressure coefficient
    sat_ss(ab) = 0.014;
%Unsaturated specific storage or pore pressure coefficient
    unsat_ss(ab) = porosity(ab)/(((weight_water*(1-water_retention_1))/(water_retention_1*water_retention_2))*(sw_un^(-1/water_retention_1))*((1-sw_un^(1/water_retention_1))^water_retention_1));
end
% Finite element input 
input.constant      = [annual_time annual_time-1 2 2 annual_time-1 1];
nodes               = annual_time; %Number of nodes
numel               = annual_time -1; %Number of element
nodaldof            = 2; %Number of degree of freedom (displacement and pore pressure)
nodesperelement     = 2; %Number of nodes per element
numat               = annual_time-1; %Number of material for nonlinear poroelasticity
ndimensions         = 1; %Number of dimension
numeq               = nodes*nodaldof;
totalelementdof     = nodaldof*nodesperelement;

%Peat physical properties data and maximum time iteration
for aa = 1:annual_time-1;
    first_column_ml(aa) = modulus(aa); 
    second_column_ml(aa) = hydraulic(aa); 
    third_column_ml(aa) = aa;
        if  wet_proportion(aa) == 1;
            fourth_column_ml(aa) = sw_sat; 
            fifth_column_ml(aa) = sat_ss(aa); 
        else
            fourth_column_ml(aa) = sw_un; 
            fifth_column_ml(aa) = unsat_ss(aa); 
        end
            sixth_column_ml(aa) = time_max; 
end
%[Young's modulus, Hydraulic conductivity, Active porosity, Degree of
%saturation, Specific storage, Maximum time simmulated]
material =[first_column_ml' second_column_ml'  third_column_ml' fourth_column_ml'  fifth_column_ml' sixth_column_ml']; 

Kg  = zeros(numeq,numeq);
Cg  = zeros(numeq,numeq);
Mg  = zeros(numeq,numeq);
f   = zeros(numeq,1);
u   = zeros(numeq,1);

%Node coordinate [node_number x1]
for bb = 1:annual_time;    
    first_column_ne(bb) = bb;
    second_column_ne(bb) = layer_elevation(bb);
    node_data = [first_column_ne' second_column_ne'];
end

%Element coordinate [element_number physical_properties node_1 node_2]
for cc = 1:annual_time-1;
    first_column_et(cc) = cc;
    second_column_et(cc) = cc; 
    third_column_et(cc) = cc;
    fourth_column_et(cc) = cc+1;
    element_data = [first_column_et' second_column_et' third_column_et' fourth_column_et'];
end

%Constrain or boundary condition [node_number vertical_diplacement pore_water_pressure]
% 1 fixed displacement or pore water pressure
% 0 non-fixed displacement or pore water pressure
for dd = 1:annual_time;
    first_column_cn(dd) = dd;
    if(dd==1);
        second_column_cn(dd) = 1;
    else
        second_column_cn(dd) = 0;
    end
    if(dd == annual_time)
        third_column_cn(dd) = 1;
    else
        third_column_cn(dd) = 0;
    end
    constrain_data = [first_column_cn' second_column_cn' third_column_cn'];
end

%Load data [node_number load_solid load_fluid ]
%Magnitude of prescribed 
%load_solid : load or displacement
%load_fluid : pore_water_pressure or flux

  for ee = 1:annual_time;
            first_column_ld(ee) = ee;
            if(ee==1);
            second_column_ld(ee) = 0;
            elseif(ee == annual_time);
            second_column_ld(ee) = -(plant_weight(ee)+layer_addition_weight(ee));
            elseif(ee == annual_time-1);
            second_column_ld(ee) = 0;
            else
            second_column_ld(ee) = 0;
            end       
            if (ee == annual_time)
            third_column_ld(ee) = 0;
            else
            third_column_ld(ee) = net_rainfall(ee);    
            end
            load_data = [first_column_ld' second_column_ld' third_column_ld'];          
 end 


%Global stifness matrix 
for iel=1:numel;
    [local_data] = localize (iel, node_data, element_data, material, input);
            ElementNo         = iel;
            MaterialArray     = local_data.mater;
            CoordinateArray   = local_data.coords;
            dofAddressArray   = local_data.dofs;
%Create local stiffness matrix             
    MaterialType = element_data(iel,2);   
    [Kl,Cl,Ml] = element(local_data);
%Assemble global stiffness matrix from local stiffness matrix   
    [Kg] = assemble (Kg, Kl, local_data, totalelementdof);
    [Cg] = assemble (Cg, Cl, local_data, totalelementdof); 
    [Mg] = assemble (Mg, Ml, local_data, totalelementdof);
    
end
%Assemble load and constrain  vectors for solid and fluid
    [f]     = vectorstretch (load_data, nodes, 1+nodaldof);  
    [con]   = vectorstretch (constrain_data,  nodes, 1+nodaldof) ; 
         u  = f.* (con);                                        
         f  = f.* (1-con);                                      
         f  = f - Kg*u;                                         

        initialvalue = 0; 
        Ktemp = Kg; 
        Ctemp = Cg; 
        Mtemp = Mg; 
        ftemp = f;       
        ftemp = f;  
        utemp = ones(numeq,1)*initialvalue;
%Reorder the mtrix       
    for ieq = numeq: -1: 1;
        if (con(ieq,1) == 1) Ktemp(:,ieq)=[]; Ktemp(ieq,:)=[]; 
                             Ctemp(:,ieq)=[]; Ctemp(ieq,:)=[]; 
                             Mtemp(:,ieq)=[]; Mtemp(ieq,:)=[]; 
                             ftemp(ieq,:)=[]; utemp(ieq,:)=[]; 
        end
    end
%Solve system equations    
     for mt = 1:annual_timestep_mechanic;                
                utemp = (Ktemp + (1/mechanic_timestep)*Ctemp)\(ftemp + (1/mechanic_timestep)*(Ctemp*utemp));
                ftemp = ftemp*0;
                warning('off')
     end
                                          
%Resubstitute for unknowns from boundary condition           
    icount=0;                                                   
    for ieq = 1:numeq;
        if (con(ieq,1) == 0) icount = icount + 1; u(ieq,1) = utemp(icount,1);
        end
    end
    u;
%Reorder the solid displacement and pore water pressure for mechanical time average to speed up the calculation.    
   [numrows_u numcols_u] = size(u);   
    for k=1:numrows_u/2;
        displacement(k)=u((2*k)-1);
        excess_pore_press(k)=u(2*k);
    end  
else
    if  annual_time < mechanic_time_average
        dissy_trans = zeros(1,annual_time);
        press_trans = zeros(1,annual_time);
    else
        dissy_trans = displacement;
        press_trans = excess_pore_press; 
    end   
    [numrows_diss numcols_diss] = size(dissy_trans); 
    for    k=1:annual_time;
      if   k < numcols_diss
           displacement(k)= dissy_trans(k);
           excess_pore_press(k)= press_trans(k);
      else
           displacement(k) = dissy_trans(numcols_diss);
           excess_pore_press(k) = press_trans(numcols_diss);
      end
    end
end
%End of the mechanical submodel 

%Update peatland height 
bog_height = layer_elevation(annual_time);

end
%Cumulative carbon calculation
cumulative_carbon = peatland_height.*mean(density).*carbon_content;
