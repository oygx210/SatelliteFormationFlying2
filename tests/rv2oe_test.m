function oe = rv2oe_test(rv, consts)

% Calculates osculating elements of the orbit
% Inputs: state vector [m, m/s]

%Output
% [a; e; i; RAAN; AOP; nu; M; u; TMA]

% Notes: for mathematical details see
%   1) Book Balk, Demin, Kunitsyn, "Sbornik zadach no nebesnoj mehanike i kosmodinamike", p. 81-82

r = rv(1:3);
v = rv(4:6);

EarthGravity = consts.muEarth;
tolerance = 1e-5;

sigma = cross(r, v);
parameter = norm(sigma)^2 / EarthGravity;
r_norm = vecnorm(r);
v_norm = vecnorm(v);
energyConstant = v_norm^2 - 2 * EarthGravity / r_norm;
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

    cosAnomaly = (parameter / r_norm - 1) / eccentricity;
    sinAnomaly =  dot(r, v) * parameter / (eccentricity * norm(sigma) * r_norm);
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
    cosAnomalyPlusPeriapsisArgument = (r(1) * cosAscendingNodeLongitude + r(2) * sinAscendingNodeLongitude) / r_norm;
    anomalyPlusPeriapsisArgument = acos(cosAnomalyPlusPeriapsisArgument);
    sinAnomalyPlusPeriapsisArgument = (cosInclination * (-sinAscendingNodeLongitude * r(1) + ...
      cosAscendingNodeLongitude * r(2)) + sinInclination * r(3)) / r_norm;
    if (sinAnomalyPlusPeriapsisArgument <= 0) 
      anomalyPlusPeriapsisArgument = 2 * pi - anomalyPlusPeriapsisArgument;
    end
    periapsisArgument = anomalyPlusPeriapsisArgument - anomaly;
    if (periapsisArgument <= 0) 
      periapsisArgument = periapsisArgument + 2 * pi;
    end
    
    [~, MeanAnomaly] = newtonnu (eccentricity, anomaly);
    ArgumentofLatitude = anomalyPlusPeriapsisArgument;
    
    TrueMeanAnomaly = mod(periapsisArgument + MeanAnomaly,2*pi);
    
    oe = [semimajorAxis;
          eccentricity;
          inclination;
          ascendingNodeLongitude;
          periapsisArgument;
          anomaly;
          MeanAnomaly;
          ArgumentofLatitude;
          TrueMeanAnomaly];
end

