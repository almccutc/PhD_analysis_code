clear all
clc

load('F:\4Vworkspace.mat','Utotal','Wtotal')

%load('F:\6Vworkspace.mat','Utotal','Wtotal')

%% 
clc
u_o = Utotal;  v_o = Wtotal;

u_mean = nanmean(u_o,3); v_mean = nanmean(v_o,3);
u_f = u_o-u_mean; v_f = v_o-v_mean;
%obj.u_rms = sqrt(nanmean((obj.u_f.^2),3)); obj.v_rms = sqrt(nanmean((obj.v_f.^2),3));
tke = 0.5*(2*(u_f.^2) + (v_f.^2));


tke_mean = nanmean(tke,3);
tke_mean = tke_mean(40:205,30:150,:);
%tke_mean = tke_mean(20:100,15:75,:);

tkeMean = mean(tke_mean(:))

figure (2)
            imagesc(tke_mean)
            colorbar
            caxis([-0.02,0.02])
            title('TKE (m^2/s^2)')