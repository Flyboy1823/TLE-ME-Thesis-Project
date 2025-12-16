function [r_t , v_t] = Orbit_Prop_Test(r0, v0, prop_time)
   
    Re=6378;  %earth's radius [km]
    u=398600; %gravitational perameter, [km^3/sec^2]
    J2=.001082; %Non spherical coefficient
    
    %% Calculation of Oribt for Earth as a Sphere
    counter=0; %used to count each itiration
    t=0;  %initial time
    r = [r0'/1000];
    v = [v0'/1000];

    r_temp=r0'/1000;  %initial position
    v_temp=v0'/1000;  %initial velocity
    for counter=2:length(prop_time)
        delta=seconds(prop_time(counter)-prop_time(counter-1));

        %Runga Kutta m4th order method for both r and v
        k1=delta*Derivative2(r(counter-1, :), v(counter-1, :)); 
        k2=delta*Derivative2(r(counter-1, :)+k1(1, :)/2, v(counter-1, :)+k1(2, :)/2);
        k3=delta*Derivative2(r(counter-1, :)+k2(1, :)/2, v(counter-1, :)+k2(2, :)/2);
        k4=delta*Derivative2(r(counter-1, :)+k3(1, :), v(counter-1, :)+k3(2, :));
        r(counter,:)=r(counter-1, :)+(k1(1,:)+2*k2(1,:)+2*k3(1,:)+k4(1,:))/6; %set next iteration of r to runga kutta approxiamtion
        v(counter,:)=v(counter-1, :)+(k1(2,:)+2*k2(2,:)+2*k3(2,:)+k4(2,:))/6;  %set next iteration of v to runga kutta approxiamtion
    end

    r_t=r'*1000;
    v_t=v'*1000;

    
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
    
  