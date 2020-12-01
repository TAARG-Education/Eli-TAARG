function [alpha_e] = Angolo_Attacco_Effettivo(lambda,r_segn,beta,dbeta,mu,psi,theta)
% ANGOLO_ATTACCO_EFFETTIVO Summary of this function goes here
%   Detailed explanation goes here
u_P = zeros(length(r_segn),length(psi));
u_T = zeros(length(r_segn),length(psi));
u_R = zeros(length(r_segn),length(psi));
phi = zeros(length(r_segn),length(psi));
alpha_e = zeros(length(r_segn),length(psi));
for i = 1:length(psi)
    u_P(:,i) = lambda + r_segn'.*dbeta(i) + beta(i).*mu.*cos(psi(i));
    u_T(:,i) = r_segn' + mu.*sin(psi(i));
    u_R(:,i) = mu*cos(psi(i));
    phi(:,i) = u_P(:,i)./u_T(:,i); 
    alpha_e(:,i) = theta' - phi(:,i);
end

alpha_stall_1 = convang(20,'deg','rad');
alpha_stall_2 = convang(-10,'deg','rad');
stall = zeros(length(r_segn),length(psi));
non_stall = NaN(length(r_segn),length(psi));
for iii = 1:length(r_segn)
    for jjj = 1:length(psi)
        if alpha_e(iii,jjj) - alpha_stall_1 > 0 || alpha_e(iii,jjj) - alpha_stall_2 < 0
            stall(iii,jjj) = alpha_e(iii,jjj);
        else
            non_stall(iii,jjj) = alpha_e(iii,jjj);
        end
    end
end

figure(1)
set(figure(1),'Color','w');
clf
fill(x2,y2,'k'); hold on; axis equal; axis off; grid off
text(0.6,0.9,['\mu = ',num2str(mu)],'Color','r','FontSize',16);
contourf(x,y,convang(non_stall','rad','deg'),25,'ShowText','on'); drawnow; axis equal;
end