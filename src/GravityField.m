classdef GravityField < handle
    properties
        Cnm_true %
        Snm_true %
        
        mu_true  %(m^3/s^2)
        radius % (m)
    end
    
    %% Constructor 
    methods
        function [self] = GravityField(radius, mu, Cnm, Snm)
            self.radius = radius;
            self.mu_true = mu;
            self.Cnm_true = Cnm;
            self.Snm_true = Snm;
        end
    end
    
    %% Public Methods
    methods (Access = public)
        % Calculate Acceleration:
        function [accel] = acceleration(self, r, E, Cnm, Snm, mu)
            if nargin == 3
                mu  = self.mu_true;
                Cnm = self.Cnm_true;
                Snm = self.Snm_true;
            elseif nargin == 5
                mu = self.mu_true;
            end

            % Body-fixed position 
            r_bf = E * r;
            
            % Calculate nmax values:
            n_max = size(Cnm,1)-1;
            m_max = n_max-1;

            % Auxiliary quantities
            d = norm(r_bf);                     % distance
            latgc = asin(r_bf(3)/d);
            lon = atan2(r_bf(2),r_bf(1));

            [pnm, dpnm] = self.Legendre(n_max, m_max, latgc);

            dUdr = 0;
            dUdlatgc = 0;
            dUdlon = 0;
            q3 = 0; q2 = q3; q1 = q2;
            for n = 0:n_max
                b1 = (-mu/d^2)*(self.radius/d)^n*(n+1);
                b2 =  (mu/d)*(self.radius/d)^n;
                b3 =  (mu/d)*(self.radius/d)^n;
                for m = 0:n
                    q1 = q1 + pnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
                    q2 = q2 + dpnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
                    q3 = q3 + m*pnm(n+1,m+1)*(Snm(n+1,m+1)*cos(m*lon)-Cnm(n+1,m+1)*sin(m*lon));
                end
                dUdr     = dUdr     + q1*b1;
                dUdlatgc = dUdlatgc + q2*b2;
                dUdlon   = dUdlon   + q3*b3;
                q3 = 0; q2 = q3; q1 = q2;
            end

            % Body-fixed acceleration
            r2xy = r_bf(1)^2+r_bf(2)^2;

            ax = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(1)-(1/r2xy*dUdlon)*r_bf(2);
            ay = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(2)+(1/r2xy*dUdlon)*r_bf(1);
            az =  1/d*dUdr*r_bf(3)+sqrt(r2xy)/d^2*dUdlatgc;

            a_bf = [ax ay az]';

            % Inertial acceleration 
            accel = E'*a_bf;
        end
    end
    
    %% Public Visualization Methods:
    methods (Access = public)
        function [self] = draw(self,Cnm,Snm,mu)
            % Check if alternative gravity field terms have been provided:
            if nargin == 1
                Cnm_use = self.Cnm_true;
                Snm_use = self.Snm_true;
                mu_use  = self.mu_true;
            else
                Cnm_use = Cnm;
                Snm_use = Snm;
                mu_use  = mu;
            end
            
            % Setup the latitude and longitude:
            latitude  = linspace(-pi/2,pi/2,180);
            longitude = linspace(0, 2*pi,360);
            
            % Evaluate the field for specified el/az:
            accel = zeros(numel(latitude),numel(longitude),3);
            for ii = 1:length(latitude)
                for jj = 1:length(longitude)
                    [x,y,z] = sph2cart(longitude(jj),latitude(ii),self.radius);
                    r = [x,y,z]';
                    accel(ii,jj,:) = self.acceleration(r, eye(3), Cnm_use,Snm_use,mu_use);
                end
            end
            
            % Draw the field:
            [~,accel_norm_field] = normw(accel);
            surf(longitude,latitude,accel_norm_field,'EdgeColor','none')
            colorbar
            axis equal
            view([0 90])
        end
    end
    
    %% Private Methods:
    methods (Access = private)
        function [pnm, dpnm] = Legendre(~,n,m,fi)
            pnm = zeros(n+1,m+1);
            dpnm = zeros(n+1,m+1);

            pnm(1,1)=1;
            dpnm(1,1)=0;
            pnm(2,2)=sqrt(3)*cos(fi);
            dpnm(2,2)=-sqrt(3)*sin(fi);
            % diagonal coefficients
            for i=2:n    
                pnm(i+1,i+1)= sqrt((2*i+1)/(2*i))*cos(fi)*pnm(i,i);
            end
            for i=2:n
                dpnm(i+1,i+1)= sqrt((2*i+1)/(2*i))*((cos(fi)*dpnm(i,i))- ...
                              (sin(fi)*pnm(i,i)));
            end
            % horizontal first step coefficients
            for i=1:n
                pnm(i+1,i)= sqrt(2*i+1)*sin(fi)*pnm(i,i);
            end
            for i=1:n
                dpnm(i+1,i)= sqrt(2*i+1)*((cos(fi)*pnm(i,i))+(sin(fi)*dpnm(i,i)));
            end
            % horizontal second step coefficients
            j=0;
            k=2;
            while(1)
                for i=k:n        
                    pnm(i+1,j+1)=sqrt((2*i+1)/((i-j)*(i+j)))*((sqrt(2*i-1)*sin(fi)*pnm(i,j+1))...
                        -(sqrt(((i+j-1)*(i-j-1))/(2*i-3))*pnm(i-1,j+1)));
                end
                j = j+1;
                k = k+1;
                if (j>m)
                    break
                end
            end
            j = 0;
            k = 2;
            while(1)
                for i=k:n        
                    dpnm(i+1,j+1)=sqrt((2*i+1)/((i-j)*(i+j)))*((sqrt(2*i-1)*sin(fi)*dpnm(i,j+1))...
                         +(sqrt(2*i-1)*cos(fi)*pnm(i,j+1))-(sqrt(((i+j-1)*(i-j-1))/(2*i-3))*dpnm(i-1,j+1)));
                end
                j = j+1;
                k = k+1;
                if (j>m)
                    break
                end
            end
        end
    end
end