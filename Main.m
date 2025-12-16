clc;close all;

%% Calculating Orbital Parameters of Desired Satellites
% all angles in radians and distance in [Km]
[a,e,i,Om,om,n,M,Epoch,rp] = prs_TLE ('46807'); 

figure
plot(Epoch,a)
xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
ylabel('$a\ [Km]$','interpreter','latex','FontSize',12')
title('$Semi\ Major\ Axis\ Vs.\ Time$','interpreter','latex','FontSize',12')
subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
grid on

figure
plot(Epoch,n)
xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
ylabel('$n\ [rad/sec]$','interpreter','latex','FontSize',12')
title('$Mean \ Motion\ Vs.\ Time$','interpreter','latex','FontSize',12')
subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
grid on


figure
plot(Epoch,rad2deg(om))
xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
ylabel('$\omega\ [{^\circ}]$','interpreter','latex','FontSize',12')
title('$Argument \ of\ Perigee\ Vs.\ Time$','interpreter','latex','FontSize',12')
subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
grid on

figure
plot(Epoch,rad2deg(Om))
xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
ylabel('$\Omega\ [{^\circ}]$','interpreter','latex','FontSize',12')
title('$RAAN\ Vs.\ Time$','interpreter','latex','FontSize',12')
subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
grid on

figure
plot(Epoch,e)
xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
ylabel('$Eccentricity$','interpreter','latex','FontSize',12')
title('$Eccentricity \ Vs.\ Time$','interpreter','latex','FontSize',12')
subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
grid on

figure
plot(Epoch,rad2deg(i))
xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
ylabel('$i\ [{^\circ}]$','interpreter','latex','FontSize',12')
title('$Inclination \ Vs.\ Time$','interpreter','latex','FontSize',12')
subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
grid on

figure
plot(Epoch,rad2deg(M))
xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
ylabel('$M\ [{^\circ}]$','interpreter','latex','FontSize',12')
title('$Mean Anomaly \ Vs.\ Time$','interpreter','latex','FontSize',12')
subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
grid on

%% Plot 3d Earth
%%%% Make into a function
%[r_satelllite v_satellite]=
% for counter=1:length(Epoch)
%     if i(counter)==0 & e(counter)==0
% 
%     else if i(counter)==0 & e(counter)~=0
% 
%     else if i(counter)~=0 & e(counter)==0

%     else
% 
% 
% 
% 
% 
% 
% 
% end





%% [a_NS47,e_NS47,i_NS47,Om_NS47,om_NS47,n_NS47,...
%     M_NS47,Epoch_NS47,rp_NS47] = prs_TLE ('26360'); %NAVSTAR 47 - Active, latest addition to the GPS fleet) 
% 
% % NAVSTAR 47 Analysis
% 
% figure
% plot(Epoch_NS47(300:end),a_NS47(300:end))
% xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
% ylabel('$a\ [Km]$','interpreter','latex','FontSize',12')
% title('$Semi\ Major\ Axis\ Vs.\ Time$','interpreter','latex','FontSize',12')
% subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
% grid on
% 
% figure
% plot(Epoch_NS47(300:end),n_NS47(300:end))
% xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
% ylabel('$n\ [rad/sec]$','interpreter','latex','FontSize',12')
% title('$Mean \ Motion\ Vs.\ Time$','interpreter','latex','FontSize',12')
% subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
% grid on
% 
% 
% figure
% plot(Epoch_NS47(300:end),rad2deg(om_NS47(300:end)))
% xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
% ylabel('$\omega\ [{^\circ}]$','interpreter','latex','FontSize',12')
% title('$Argument \ of\ Perigee\ Vs.\ Time$','interpreter','latex','FontSize',12')
% subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
% grid on
% 
% figure
% plot(Epoch_NS47(300:end),rad2deg(Om_NS47(300:end)))
% xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
% ylabel('$\Omega\ [{^\circ}]$','interpreter','latex','FontSize',12')
% title('$RAAN\ Vs.\ Time$','interpreter','latex','FontSize',12')
% subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
% grid on
% 
% figure
% plot(Epoch_NS47(300:end),e_NS47(300:end))
% xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
% ylabel('$Eccentricity$','interpreter','latex','FontSize',12')
% title('$Eccentricity \ Vs.\ Time$','interpreter','latex','FontSize',12')
% subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
% grid on
% 
% figure
% plot(Epoch_NS47(300:end),rad2deg(i_NS47(300:end)))
% xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
% ylabel('$i\ [{^\circ}]$','interpreter','latex','FontSize',12')
% title('$Inclination \ Vs.\ Time$','interpreter','latex','FontSize',12')
% subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
% grid on
% 
% figure
% plot(Epoch_NS47(300:end),rad2deg(M_NS47(300:end)))
% xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
% ylabel('$M\ [{^\circ}]$','interpreter','latex','FontSize',12')
% title('$Mean Anomaly \ Vs.\ Time$','interpreter','latex','FontSize',12')
% subtitle('$NAVSTAR\ 47$' ,'interpreter','latex','FontSize',9)
% grid on
% 
% 

%% NAVSTAR 43 Analysis
% 
% [a,e,i,Om,om,n...
%     ,M,Epoch,rp] = prs_TLE ('24876'); %NAVSTAR 43 
% 
% 
% figure
% plot(Epoch(300:end),a(300:end))
% xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
% ylabel('$a\ [Km]$','interpreter','latex','FontSize',12')
% title('$Semi\ Major\ Axis\ Vs.\ Time$','interpreter','latex','FontSize',12')
% subtitle('$NAVSTAR\ 43$' ,'interpreter','latex','FontSize',9)
% grid on
% 
% 
% figure
% plot(Epoch(300:end),rad2deg(om(300:end)))
% xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
% ylabel('$\omega\ [{^\circ}]$','interpreter','latex','FontSize',12')
% title('$Argument \ of\ Perigee\ Vs.\ Time$','interpreter','latex','FontSize',12')
% subtitle('$NAVSTAR\ 43$' ,'interpreter','latex','FontSize',9)
% grid on
% 
% figure
% plot(Epoch(300:end),rad2deg(Om(300:end)))
% xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
% ylabel('$\Omega\ [{^\circ}]$','interpreter','latex','FontSize',12')
% title('$RAAN\ Vs.\ Time$','interpreter','latex','FontSize',12')
% subtitle('$NAVSTAR\ 43$' ,'interpreter','latex','FontSize',9)
% grid on
% 
% figure
% plot(Epoch(300:length(e)),e(300:end))
% xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
% ylabel('$Eccentricity$','interpreter','latex','FontSize',12')
% title('$Eccentricity \ Vs.\ Time$','interpreter','latex','FontSize',12')
% subtitle('$NAVSTAR\ 43$' ,'interpreter','latex','FontSize',9)
% grid on
% 
% figure
% plot(Epoch(300:end),rad2deg(i(300:end)))
% xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
% ylabel('$i\ [{^\circ}]$','interpreter','latex','FontSize',12')
% title('$Inclination \ Vs.\ Time$','interpreter','latex','FontSize',12')
% subtitle('$NAVSTAR\ 43$' ,'interpreter','latex','FontSize',9)
% grid on
% 
% figure
% plot(Epoch(300:end),rad2deg(M(300:end)))
% xlabel('$Epoch\ [UTC]$','interpreter','latex','FontSize',12')
% ylabel('$M\ [{^\circ}]$','interpreter','latex','FontSize',12')
% title('$Mean Anomaly \ Vs.\ Time$','interpreter','latex','FontSize',12')
% subtitle('$NAVSTAR\ 43$' ,'interpreter','latex','FontSize',9)
% grid on


%% Useful code for interpolation and Epoch Merging (not used in the report)
% %% Merging the Epoch vectors
% 
% merged_Epoch = zeros(1,length(Epoch_NS47)+length(Epoch_NS55));
% cnt_47 = 1;
% cnt_55 = 1;
% 
% for i=1:length(merged_Epoch)
%     if cnt_47>length(Epoch_NS47)
%         merged_Epoch(i) = Epoch_NS55(cnt_55);
%         cnt_55 = cnt_55+1;
%     elseif cnt_55>length(Epoch_NS55)
%         merged_Epoch(i) = Epoch_NS47(cnt_47);
%         cnt_47 = cnt_47+1;
%     elseif Epoch_NS47(cnt_47)<= Epoch_NS55(cnt_55)
%         merged_Epoch(i) = Epoch_NS47(cnt_47);
%         cnt_47 = cnt_47+1;
%     elseif Epoch_NS55(cnt_55)<= Epoch_NS47(cnt_47)
%         merged_Epoch(i) = Epoch_NS55(cnt_55);
%         cnt_55 = cnt_55+1;
%     end
% end




% [x_NS47,y_NS47,z_NS47] = orb2cart(a_NS47,e_NS47,i_NS47,Om_NS47,om_NS47,M_NS47);
% [x_NS55,y_NS55,z_NS55] = orb2cart(a_NS55,e_NS55,i_NS55,Om_NS55,om_NS55,M_NS55);

% N = round(length(merged_Epoch)/length(a_NS47));
% M = round(length(merged_Epoch)/length(a_NS55));
% 
% a_NS47_interp = interp(a_NS47,N);e_NS47_interp = interp(e_NS47,N);i_NS47_interp = interp(i_NS47,N);...
%     Om_NS47_interp = interp(Om_NS47,N); om_NS47_interp = interp(om_NS47,N);M_NS47_interp = interp(M_NS47,N);
% 
% [x_NS47_interp,y_NS47_interp,z_NS47_interp] = orb2cart(a_NS47_interp,e_NS47_interp,i_NS47_interp,Om_NS47_interp,om_NS47_interp,M_NS47_interp);
% 
% a_NS55_interp = interp(a_NS55,M);e_NS55_interp = interp(e_NS55,M);i_NS55_interp = interp(i_NS55,M);...
%     Om_NS55_interp = interp(Om_NS55,M); om_NS55_interp = interp(om_NS55,M);M_NS55_interp = interp(M_NS55,M);
% [x_NS55_interp,y_NS55_interp,z_NS55_interp] = orb2cart(a_NS55_interp,e_NS55_interp,i_NS55_interp,Om_NS55_interp,om_NS55_interp,M_NS55_interp);


%interpolation integrity test
% figure()
% plot(Epoch_NS47,a_NS47)
% hold on
% plot(merged_Epoch,a_NS47_interp(1:length(merged_Epoch)))
% legend('initial','interpolated')
% grid on

% interp_x_NS47 = interp(x_NS47,N);interp_y_NS47 = interp(y_NS47,N);...
%     interp_z_NS47 = interp(z_NS47,N);
% interp_x_NS55 = interp(x_NS55,M);interp_y_NS55 = interp(y_NS55,M);...
%     interp_z_NS55 = interp(z_NS55,M);
% interp_x_NS47_b = interp_x_NS47(1:length(interp_x_NS55));interp_y_NS47_b = interp_y_NS47(1:length(interp_y_NS55));...
%     interp_z_NS47_b = interp_z_NS47(1:length(interp_z_NS55));
% 
% rel_x = x_NS47_interp(1:length(x_NS55_interp))-x_NS55_interp;
% rel_y = y_NS47_interp(1:length(y_NS55_interp))-y_NS55_interp;
% rel_z = z_NS47_interp(1:length(z_NS55_interp))-z_NS55_interp;
% 
% rel_d = sqrt((rel_x).^2+(rel_y).^2+(rel_z).^2);






