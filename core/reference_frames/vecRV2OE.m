%                  function rv2oe   
%  this function finds the orbital elements given the geocentric
%    equatorial position and velocity vectors.
%
%
% 
%  inputs          description                    range / units
%    rv          - ijk position;velocity vector   m; m/s
%    consts      - standard set of constants for orbital mechanics simulation
%
%  outputs       :
%    a           - semimajor axis                 km
%    ecc         - eccentricity
%    incl        - inclination                    0.0  to pi rad
%    omega       - longitude of ascending node    0.0  to 2pi rad
%    argp        - argument of perigee            0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%    M           - mean anomaly                   0.0  to 2pi rad
%    u           - argument of latitude           0.0  to 2pi rad
%
%  locals        :
%    h = cross(r,v)        - angular momentum h vector      km2 / s
%    ebar        - eccentricity     e vector
%    nbar        - line of nodes    n vector
%    c1          - v**2 - u/r
%    rdotv       - r dot v
%    hk          - hk unit vector
%    sme         - specfic mechanical energy      km2 / s2
%
%  coupling      :
%    mag         - magnitude of a vector
%    angl        - find the angl between two vectors
%    newtonnu    - find the mean anomaly
%
%  references    :
%    vallado       2007, 121, alg 9, ex 2-5

function oe = vecRV2OE(rv, consts)
        b = size(rv,2);

        r = rv(1:3,:);
        v = rv(4:6,:);

        r_norm = vecnorm(r);
        v_norm = vecnorm(v);
        % ------------------  find h n and e vectors   ----------------
        h = cross(r,v);
        h_norm = vecnorm(h);

        n = [-h(2,:) ; h(1,:); zeros(1,b)];
        n_norm = vecnorm(n);
        
        e = (v_norm.^2 - consts.muEarth./r_norm).*r ./consts.muEarth - dot(r,v).*v ./ consts.muEarth;
        ecc = vecnorm(e);
        if ecc > 1
            disp('hyperbolic trajectory');
        end
        % ------------------  find semi-major axis   ----------------
        sme = (v_norm .* v_norm *0.5) - (consts.muEarth ./ r_norm);
        a = -consts.muEarth  ./ (2 * sme);
        % -----------------  find inclination   -------------------
        hk = h(3,:)./h_norm;
        incl= acos(hk);
        % ----------  find longitude of ascending node ------------
        omega= atan2(n(2,:), n(1,:));
        % ---------------- find argument of perigee ---------------
        argp = acos(dot(n,e)./n_norm./ecc);
        check1 = (e(3,:) < 0);
        check2 = (e(3,:) > 0);
        argp= 2*pi.*check1 - argp.*check1 + check2.*argp; 
        % ------------  find true anomaly at epoch    -------------
        nu =  real(acos( dot(e,r)./ecc./r_norm));
        check1 = (dot(r,v) < 0 );
        check2 = (dot(r,v) > 0 );
        nu= 2*pi.*check1 - nu.*check1 + nu.*check2;
 
        M = mean_anomaly(ecc,nu);     

        % ------------  find argument of latitude at epoch    -------------
        u = acos(dot(n,r)./n_norm./r_norm);
        check1 = (r(3,:) < 0);
        check2 = (r(3,:) > 0);
        u = 2*pi.*check1 - u.*check1 + u.*check2;

        oe = [a ; ecc ; incl ; omega ; argp ; nu ; M ; u];

