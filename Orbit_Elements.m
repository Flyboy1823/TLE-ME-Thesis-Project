function Orbit_Elements(a, e, i, w, Omega, M, TF)
   
    hold on
    Re=6378;  %earth's radius [km]
    u=398600; %gravitational perameter, [km^3/sec^2]
    J2=.001082; %Non spherical coefficient
    
    %%  Calculate Orbital Charachteristics
    h=cross(r0, v0); %specific relative angular momentum [km/s^2]
    e=norm((cross(v0, h)/u)-(r0/norm(r0))); %eccentricity
    p=((norm(h))^2)/u;  %orbit's parameter
    ra=p/(1-e); %apogee [km]
    rp=p/(1+e); %perigee [km]
    a=(ra+rp)/2; %semi major axis [km]
    
    %% Calculation of Oribt for Earth as a Sphere
    counter=0; %used to count each itiration
    delta=100; %itiration step size
    t=0;  %initial time
    r_temp=r0;  %initial position
    v_temp=v0;  %initial velocity
    while t<=TF
        counter=counter+1;
        time(counter)=t;
        r(counter, :)=r_temp;   %next position iteratio
        v(counter, :)=v_temp;   %next velocit iteration
        epsi(counter)=norm(v(counter, :))^2/2 - u/norm(r(counter, :)); %specfic orbital energy
        H(counter)=norm(cross(r(counter, :), v(counter, :))); %magnitude of specific relative angular momentum
        R(counter)=norm(r(counter, :)); %Radius of motion (magnitiude)
        V(counter)=norm(v(counter, :)); %Velocity of motion (magnitude)
        
        if dot(r(counter, :), v(counter, :))>0   %to calculate the true anomoly it is dependent on the side of the graph its on
            phi(counter)=acos(H(counter)/(R(counter)*V(counter)))*180/pi;
        else
            phi(counter)=-acos(H(counter)/(R(counter)*V(counter)))*180/pi;
        end
        %Runga Kutta m4th order method for both r and v
        k1=delta*Derivative(r(counter, :), v(counter, :)); 
        k2=delta*Derivative(r(counter, :)+k1(1, :)/2, v(counter, :)+k1(2, :)/2);
        k3=delta*Derivative(r(counter, :)+k2(1, :)/2, v(counter, :)+k2(2, :)/2);
        k4=delta*Derivative(r(counter, :)+k3(1, :), v(counter, :)+k3(2, :));
        t=t+delta; %increase itteration of time
        r_temp=r(counter, :)+(k1(1,:)+2*k2(1,:)+2*k3(1,:)+k4(1,:))/6; %set next iteration of r to runga kutta approxiamtion
        v_temp=v(counter, :)+(k1(2,:)+2*k2(2,:)+2*k3(2,:)+k4(2,:))/6;  %set next iteration of v to runga kutta approxiamtion
    end
    
    %Plot Sphere with Radius Re
    earth_sphere('km');
    axis equal
    hold on
    % [x, y, z] =sphere;
    % surf(x*Re, y*Re, z*Re);
    % axis equal
    
    hold on
    
    %Plot orbit
    plot3(r(:, 1), r(:, 2), r(:, 3));
    title('Oribt of Satelite Around Earth');
    xlabel('x in ECI [km]');
    ylabel('y in ECI [km]');
    zlabel('z in ECI [km]');
    
    hold off
    
    
    
    %% Derivative Function
    function M = Derivative(x, y)
    Re=6378;  %earth's radius [km]
    u=398600; %gravitational perameter, [km^3/sec^2]
    J2=.001082; %Non spherical coefficient
    M(1, :)=y;   %x dot, y dot and v dot are already given in the velocity vector
    M(2, :)=-u*(x/(norm(x)^3));  %accoridng to equation r double dot = -u* r hat / r magnitude ^3
    end
    
    %% Derivative Function or earth not as a sphere
    function M = Derivative2(x, y)
    Re=6378;  %earth's radius [km]
    u=398600; %gravitational perameter, [km^3/sec^2]
    J2=.001082; %Non spherical coefficient
    M(1, :)=y;   %x dot, y dot and v dot are already given in the velocity vector
    M(2, 1)=(-u*(x(1)/(norm(x)^3)))*(1-(3/2)*J2*((Re/norm(x))^2)*(5*((x(3)^2)/(norm(x)^2))-1));  %accoridng to equation given in the homeowork sheet for non spherical earth
    M(2, 2)=(-u*(x(2)/(norm(x)^3)))*(1-(3/2)*J2*((Re/norm(x))^2)*(5*((x(3)^2)/(norm(x)^2))-1));  %accoridng to equation given in the homeowork sheet for non spherical earth
    M(2, 3)=(-u*(x(3)/(norm(x)^3)))*(1-(3/2)*J2*((Re/norm(x))^2)*(5*((x(3)^2)/(norm(x)^2))-3));  %accoridng to equation given in the homeowork sheet for non spherical earth
    end
end
    
    
    
    
    
    
    
    
    
    











