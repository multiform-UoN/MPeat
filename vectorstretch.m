function [ newvector ] = vectorstretch( oldvector, rows, cols )
% Remove nodal data from (first column) load vectors and to 
% transform them in to a single load vector for multiple dof systems.

    for i=1:rows;
        for j=2:cols;
        newvector((i-1)*(cols-1)+j-1 ,1) = oldvector(i,j);
        end
    end
