function [matchMatrix, satFuel] = maneuverAssignment_2(costMatrix, currentFuel)
    % in maneuverAssignment_2 we assign the satellite with the smallest
    % amount of fuel to the cheapest trajectory. For the rest we solve
    % assignment problem with maximin optimization

    F = currentFuel; % Remaining level of fuel in formation satellite
    C = costMatrix;
    penalty = sum(abs(F))+sum(abs(C),'all'); %penalty constant
    
    satFuel = F;
    %% assigning a satellite with smallest fuel level to the cheapest trajectory 
    [~, sat_min]  = min(F);
    [~, traj_min] = min(C(sat_min,:));
    
    for i = 1:length(currentFuel)
        if i ~= traj_min
            C(sat_min,i) = penalty;
        end
        if i ~= sat_min
            C(i,traj_min) = penalty;
        end
    end
    
    while (min(F) > 0)
        A = matchpairs(C, penalty); %solving assignment problem with the objective of minimizing overall fuel spent for maneuver

        for i = 1:length(F) %updating the remaining fuel
            F(A(i,1)) = currentFuel(A(i,1)) - C(A(i,1),A(i,2));
        end

        min_fuel_left = min(F); %finding the sat with the lowest fuel level

        if min_fuel_left > 0
            for i = 1:size(C,1) %updating costs so as to forbid maneuvers we do not want to consider
            for j = 1:size(C,2)
                C(i, j) = max(C(i, j), (currentFuel(i) - C(i,j) <= min_fuel_left) * penalty);
            end
            end
            
            satFuel = F;
            matchMatrix = sortrows(A,1);
            
        end
    end
end