function [matchMatrix, satFuel] = maneuverAssignment(costMatrix, currentFuel)

    F = currentFuel; % Remaining level of fuel in formation satellite
    C = costMatrix;
    penalty = sum(abs(F))+sum(abs(C),'all'); %penalty constant

    satFuel = F;
    mean_fuel_consumption = [];
    Fmin = [];
    assignment = [];
    
    while (min(F) > 0)
        A = matchpairs(C, penalty); %solving assignment problem with the objective of minimizing overall fuel spent for maneuver

        for i = 1:length(F) %updating the remaining fuel
            F(A(i,1)) = currentFuel(A(i,1)) - C(A(i,1),A(i,2));
        end

        [min_fuel_left, index_min_fuel_left] = min(F); %finding the sat with the lowest fuel level

        if min_fuel_left > 0
            for i = 1:size(C,1) %updating costs so as to forbid maneuvers we do not want to consider
            for j = 1:size(C,2)
                C(i, j) = max(C(i, j), (currentFuel(i) - C(i,j) <= min_fuel_left) * penalty);
            end
            end
            
            satFuel = F;
            matchMatrix = sortrows(A,1);
            mean_fuel_consumption = [mean_fuel_consumption; mean(currentFuel - F)];
            Fmin = [Fmin; min_fuel_left];
            assignment = [assignment; matchMatrix(index_min_fuel_left,:)];

            maximin.mean_fuel_consumption = mean_fuel_consumption;
            maximin.Fmin = Fmin;
            maximin.assignment = assignment;
            save('C:\SatelliteFormationFlying\data\maximin_optimization', 'maximin');
        end
    end
end