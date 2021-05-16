function [local_data] = localize (iel, node_data, element_data, material, input)
% Localize global data at the local element level
% Define global system variables
  nodes           = input.constant(1);
  numel           = input.constant(2);
  nodaldof        = input.constant(3);
  nodesperelement = input.constant(4);
  numat           = input.constant(5);
  ndimensions     = input.constant(6);


% ----------------------------------------------Zero arrays 
  local_data.mater   = zeros(1,6);
  local_data.coords  = zeros(nodesperelement,ndimensions);
  local_data.dofs    = zeros(1,nodesperelement*nodaldof);
% ----------------------------------------------Add local data to arrays
    material_number = element_data(iel,2);  
%----------------------------------------------- Material data
    local_data.mater   = material(material_number, 1:6);                      
%----------------------------------------------- Coordinate data
        for inodes=1:nodesperelement                                              
            local_data.coords(inodes,1:ndimensions) = node_data(element_data(iel,2+inodes),2:ndimensions+1); 
        end
%----------------------------------------------- degree of freedom  
        for inodes=1:nodesperelement
            globalnode = element_data(iel,2+inodes);
            startdof   = (globalnode-1)*nodaldof ;                              
            for idof=1:nodaldof
            local_data.dofs(1,(inodes-1)*nodaldof+idof) = startdof+idof;  
            end
        end