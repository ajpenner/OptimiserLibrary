function P = BivariateNormalCDF(b1,b2,rho)
% returns the bivariate normal probability
% BVN(b1,b2,rho) = int_-inf ^b1 int_-inf ^b2 1/(2pi sqrt(1-rho^2) * exp[-(x^2-2
% rho x y + y^2)/(2(1-rho))] dx dy
%
% Based on algorithms by Drezner & Wesolowsky (1990), 
% as reported by Alan Genz (2004).
% Philipp Hennig, January 2012

%persistent x1 x2 x3 x4 x5 XGL w1 w2 w3 w4 w5 WGL% Gauss-Legendre points and weights
persistent XGL WGL sq2 twopi% Gauss-Legendre points and weights
if isempty(WGL)
%     x1 = 0;
%     x2 = +1/3 * sqrt(5 - 2*sqrt(10/7));  % +0.5385
%     x3 = -1/3 * sqrt(5 - 2*sqrt(10/7));  % -0.5385
%     x4 = +1/3 * sqrt(5 + 2*sqrt(10/7));  % +0.9062
%     x5 = -1/3 * sqrt(5 + 2*sqrt(10/7));  % -0.9062
%     XGL= [x1; x2; x3; x4; x5];
%     w1 = 128 / 225;
%     w2 = (322 + 13 * sqrt(70)) / 900;
%     w3 = w2;
%     w4 = (322 - 13 * sqrt(70)) / 900;
%     w5 = w4;
%     WGL= [w1; w2; w3; w4; w5];
      x1 = 0.095012509837637;
      x2 = 0.281603550779259;
      x3 = 0.458016777657227;
      x4 = 0.617876244402644;
      x5 = 0.755404408355003;
      x6 = 0.865631202387832;
      x7 = 0.944575023073233;
      x8 = 0.989400934991650;
      XGL = [x1;-x1;x2;-x2;x3;-x3;x4;-x4;x5;-x5;x6;-x6;x7;-x7;x8;-x8];
      w1 = 0.189450610455069;
      w2 = 0.182603415044924;
      w3 = 0.169156519395003;
      w4 = 0.149595988816577;
      w5 = 0.124628971255534;
      w6 = 0.095158511682493;
      w7 = 0.062253523938648;
      w8 = 0.027152459411754;
      WGL = [w1;w1;w2;w2;w3;w3;w4;w4;w5;w5;w6;w6;w7;w7;w8;w8];
      sq2 = sqrt(2); twopi = 2*pi;
end

if rho == 0 % trivial case
    % P = GaussCDF(b1) * GaussCDF(b2);
    P = 0.25 * (1 + erf(b1 / sqrt(2))) * (1 + erf(b2 / sqrt(2)));
    return;
end

switch abs(rho) <= 0.99
    case 1 % substitute integration over r=[0,rho] with r=sin(theta)
        %t1  = GaussCDF(b1) * GaussCDF(b2); % first term
        t1  = 0.25 * (1 + erf(b1 / sq2)) * (1 + erf(b2 / sq2));
        upl = asin(rho); % upper integration limit
        xev = upl/2 .* (XGL + 1); % evaluation points for Gauss-Legendre approximation
        t2  = upl/2 .* WGL' * exp(-(b1.*b1 + b2.*b2 - 2 .* b1 .* b2 .* sin(xev))...
            ./(2.*cos(xev).*cos(xev)));
        P   = t1 + t2 / twopi;
    case 0 % different integration region, to avoid singularity in integrand
        h = -b1; k = -b2; s = sign(rho);  % classic co-ordinate system
        if s == 1; 
            Lhks = 0.5 * (1 + erf(-max(h,k) / sq2)); % GaussCDF(-max(h,k));
        else
            %Lhks = max(0,GaussCDF(-h)-GaussCDF(k)); 
            Lhks = max(0,0.5* (erf(-h / sqrt(2)) - erf(k/sq2)));
        end;
        a   = sqrt(1-rho*rho); b = abs(h-s*k); c = (4-s*h*k)/8; % some constants
        A   = s / (2*pi) * exp(- s*h*k/2);
        B   = a * (1-c*(b*b-a*a)/3) * exp(-b*b/(2*a*a))...
               - b*(1-c*b*b/3) * sqrt(2*pi)*0.5 * (1 + erf(-b/a / sq2));
        xev = a / 2 .* (XGL+1);
        xev2= xev.*xev;
        FVALS = exp(-((h-s*k)*(h-s*k)) ./ (2 .* xev2)) ...
            .* (exp(-(s*h*k)./(1+ sqrt(1-xev2))) ./ sqrt(1-xev2) ...
            - exp(-s*h*k/2) .* (1 + (4-s*h*k).* xev2 / 8));
        C   = a/2 .* WGL' * FVALS;
        P = Lhks - A * B - (s/twopi) * C;
end
