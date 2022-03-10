function [osculatingElements] = calculateOsculatingElements(rv, consts)

% Calculates osculating elements of the orbit
% 
% Inputs:
%     r
%       : 1x3 radius-vector of current position, in [m]
%     v
%       : 1x3 velocity of current position, in [m/s]
%                  
% Outputs:
%     osculatingElements
%       : structure, which contains values of osculating elements, namely
%         semimajorAxis
%           : semimajor axis of orbit, in [m]
%         eccentricity
%           : eccentricity of orbit
%         periapsisArgument
%           : periapsis argument of orbit, in degrees
%         inclination
%           : inclination of orbit, in degrees
%         ascending node longitude
%           : ascending node longitude of orbit, in degrees
%         periapsis
%           : periapsis of orbit, in [m]
%         apogee
%           : apogee of orbit, in [m]
%         period
%           : period of orbit, in [s]
%         
%
% Notes: for mathematical details see
%   1) Book Balk, Demin, Kunitsyn, "Sbornik zadach no nebesnoj mehanike i kosmodinamike", p. 81-82
%   2) Reports
%
    r = rv(1:3);
    v = rv(4:6);
    EarthGravity = consts.muEarth;
    rad2deg = 180/pi;
    tolerance = 1e-7;
  
  sigma = cross(r, v);
  parameter = vecnorm(sigma)^2 / EarthGravity;
  altitude = vecnorm(r);
  energyConstant = vecnorm(v)^2 - 2 * EarthGravity / altitude;
  eccentricity2 = 1 + energyConstant * parameter / EarthGravity;
  
  if eccentricity2 >= 1
    warning('Eccentricity is greater than 1: the orbit is not elliptic!');
  end
  
  if (eccentricity2 < 0)
    eccentricity = 0;
  else
    eccentricity = sqrt(eccentricity2);
  end
  
  semimajorAxis  = parameter / (1 - eccentricity2);
  periapsis = parameter / (1 + eccentricity);
  apogee = parameter / (1 - eccentricity);
  period = 2 * pi * sqrt(semimajorAxis^3 / EarthGravity);
  
  eSigma = sigma / norm(sigma);
  cosInclination = eSigma(3);
  if (cosInclination >= 1 - tolerance)
    inclination = 0;
    ascendingNodeLongitude = 0;
    cosInclination = 1;
    sinInclination = 0;
    cosAscendingNodeLongitude = 1;
    sinAscendingNodeLongitude = 0;        
  else
    inclination = acos(cosInclination);
    sinInclination = sin(inclination);
    
    cosAscendingNodeLongitude = -eSigma(2) / sinInclination;
    if (cosAscendingNodeLongitude >= 1)
       ascendingNodeLongitude = 0;
       cosAscendingNodeLongitude = 1;
    elseif (cosAscendingNodeLongitude <= -1)
       ascendingNodeLongitude = pi;
       cosAscendingNodeLongitude = -1;
    else
      ascendingNodeLongitude = acos(cosAscendingNodeLongitude);
    end
    if (eSigma(1) <= tolerance) 
      ascendingNodeLongitude = 2 * pi - ascendingNodeLongitude;
    end
    sinAscendingNodeLongitude = sin(ascendingNodeLongitude);
  end

  if (eccentricity <= tolerance)
    periapsisArgument = 0; 
  else
    cosAnomaly = (parameter / altitude - 1) / eccentricity;
    sinAnomaly =  dot(r, v) * parameter / (eccentricity * norm(sigma) * altitude);
    if cosAnomaly < -1
      anomaly = pi;
    elseif cosAnomaly > 1
      anomaly = 0;
    else
      anomaly = acos(cosAnomaly);
    end
    if (sinAnomaly <= 0) 
      anomaly = 2 * pi - anomaly;
    end
    cosAnomalyPlusPeriapsisArgument = (r(1) * cosAscendingNodeLongitude + r(2) * sinAscendingNodeLongitude) / altitude;
    anomalyPlusPeriapsisArgument = acos(cosAnomalyPlusPeriapsisArgument);
    sinAnomalyPlusPeriapsisArgument = (cosInclination * (-sinAscendingNodeLongitude * r(1) + ...
      cosAscendingNodeLongitude * r(2)) + sinInclination * r(3)) / altitude;
    if (sinAnomalyPlusPeriapsisArgument <= 0) 
      anomalyPlusPeriapsisArgument = 2 * pi - anomalyPlusPeriapsisArgument;
    end
    periapsisArgument = anomalyPlusPeriapsisArgument - anomaly;
    if (periapsisArgument <= 0) 
      periapsisArgument = periapsisArgument + 2 * pi;
    end
  end
  
  [~,meanAnomaly] = newtonnu (eccentricity, anomaly);
  %% conversion from radians to degrees
  
%   osculatingElements.semimajorAxis = semimajorAxis;
%   osculatingElements.eccentricity = eccentricity;
%   osculatingElements.periapsisArgument = periapsisArgument;
%   osculatingElements.inclination = inclination;
%   osculatingElements.ascendingNodeLongitude = ascendingNodeLongitude;
%   osculatingElements.periapsis = periapsis;
%   osculatingElements.apogee = apogee;
%   osculatingElements.period = period;
% 
  osculatingElements = zeros(8,1);
  osculatingElements(1) = semimajorAxis;
  osculatingElements(2) = eccentricity;
  osculatingElements(3) = inclination;
  osculatingElements(4) = ascendingNodeLongitude;
  osculatingElements(5) = periapsisArgument;
  osculatingElements(6) = anomaly;
  osculatingElements(7) = meanAnomaly;
  osculatingElements(8) = anomalyPlusPeriapsisArgument;
    
end