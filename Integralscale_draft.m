
%% autocorrelation calculation
load('samplePIV_data_AM','u_original','v_original','calxy')

u_o = u_original(:,32:217,:).*100; v_o = v_original(:,32:217,:).*100; %calibrated to cm/s
obj.calibration = calxy;
u_o = permute(u_o,[3,1,2]);
v_o = permute(v_o,[3,1,2]);

subwindow=16;    %pixels between subwindow centers 
[~, Ny, Nx]=size(u_o); 

x_c=(Nx+1)/2; % determines the centerline position
y_c=(Ny+1)/2; % determines the centerline position
rad_x=[0,(subwindow:subwindow*2:Nx*subwindow).*obj.calibration]; %calibrate to cm
rad_y=[0,(subwindow:subwindow*2:Ny*subwindow).*obj.calibration]; %calibrate to cm
heights=([-Ny/2:1:-1 1:1:Ny/2]).*obj.calibration.*subwindow; %calibrate to cm
widths=([-Nx/2:1:-1 1:1:Nx/2]).*obj.calibration.*subwindow; %calibrate to cm

%preallocate matrices for speed
obj.a_u_11_1 = [ones([Ny 1]) NaN([Ny Nx/2])]; obj.a_v_33_1 = [ones([Ny 1]) NaN([Ny Nx/2])];
obj.a_u_11_3 = [ones([Nx 1]) NaN([Nx Ny/2])]; obj.a_v_33_3 = [ones([Nx 1]) NaN([Nx Ny/2])];

 
%spatial, horizontal autocorrelation calculation, for the case of an even # of vertical and horizontal subwindows          
for row=1:Ny % calculated at every height for each subwindow 
    for radius=0.5:1:x_c-0.5

        %11,1 - longitudinal, Horizontal velocity, horizontal separation
        obj.a_u_11_1(row,radius+1.5)=mean(u_o(:,row,x_c-radius).*u_o(:,row,x_c+radius),'omitnan')...
            ./sqrt(mean(u_o(:,row,x_c-radius).^2,'omitnan').*mean(u_o(:,row,x_c+radius).^2,'omitnan'));
        
        %33,1 - transverse, Vertical velocity, horizontal separation
        obj.a_v_33_1(row,radius+1.5)=mean(v_o(:,row,x_c-radius).*v_o(:,row,x_c+radius),'omitnan')...
            ./sqrt(mean(v_o(:,row,x_c-radius).^2,'omitnan').*mean(v_o(:,row,x_c+radius).^2,'omitnan'));
    end
end

%spatial, vertical autocorrelation calculation, for the case of an even # of vertical and horizontal subwindows          
for column=1:Nx % calculated at every width for each subwindow 
    for radius=0.5:1:y_c-0.5

        %11,3 - longitudinal, Horizontal velocity, vertical separation
        obj.a_u_11_3(column, radius+1.5)=mean(u_o(:,y_c-radius,column).*u_o(:,y_c+radius,column),'omitnan')...
            ./sqrt(mean(u_o(:,y_c-radius,column).^2,'omitnan').*mean(u_o(:,y_c+radius,column).^2,'omitnan'));
        
        %33,3 - transverse, Vertical velocity, vertical separation
        obj.a_v_33_3(column,radius+1.5)=mean(v_o(:,y_c-radius,column).*v_o(:,y_c+radius,column),'omitnan')...
            ./sqrt(mean(v_o(:,y_c-radius,column).^2,'omitnan').*mean(v_o(:,y_c+radius,column).^2,'omitnan'));
    end
end

%calculate L at all heights
for row=1:Ny
    
    % Length scale calculation
    newExpFuncG=@(rr,ll) exp(-rr./ll); 
    %newExpFuncG_traverse=@(rr,ll) exp(-rr./ll).*(1-(rr./(2.*ll))); 

    L0=6; % Starting guess for L

    g_11_1=@(ll)sum((obj.a_u_11_1(row,:)-newExpFuncG(rad_x,ll)).^2,2);
    g_33_1=@(ll)sum((obj.a_v_33_1(row,:)-newExpFuncG(rad_x,ll)).^2,2);

    % Minimize to find Lstar
    L_11_1(row,1)=fminunc(g_11_1,L0); 
    L_33_1(row,1)=fminunc(g_33_1,L0); 

    exp_11_1(row,:)=newExpFuncG(rad_x,L_11_1(row,1)); 
    exp_33_1(row,:)=newExpFuncG(rad_x,L_33_1(row,1)); 
end    

%calculate L at all widths
for column=1:Nx
    
    % Length scale calculation
    newExpFuncG=@(rr,ll) exp(-rr./ll); 
    %newExpFuncG_traverse=@(rr,ll) exp(-rr./ll).*(1-(rr./(2.*ll))); 

    L0=6; % Starting guess for L

    g_11_3=@(ll)sum((obj.a_u_11_3(column,:)-newExpFuncG(rad_y,ll)).^2,2);
    g_33_3=@(ll)sum((obj.a_v_33_3(column,:)-newExpFuncG(rad_y,ll)).^2,2);

    % Minimize to find Lstar
    L_11_3(column,1)=fminunc(g_11_3,L0); 
    L_33_3(column,1)=fminunc(g_33_3,L0); 

    exp_11_3(column,:)=newExpFuncG(rad_y,L_11_3(column,1)); 
    exp_33_3(column,:)=newExpFuncG(rad_y,L_33_3(column,1)); 
end   

    % autocorrelation plots
    figure(1)
    subplot(2,2,1)
    plot(rad_x,obj.a_u_11_1(Ny/2,:),'b-*',rad_x,exp_11_1(Ny/2,:),'r-o');
    grid on; 
    legend('Actual','Predicted','Location','northeastoutside');
    title('11,1 Autocorrelation - horizontal center','Interpreter','Latex', 'FontSize',14)
    ylabel('a (r)'); xlabel('r (cm)')
    xlim([0 max(rad_x)]); ylim([-1 1]);

    subplot(2,2,2)
    plot(rad_x,obj.a_v_33_1(Ny/2,:),'b-*',rad_x,exp_33_1(Ny/2,:),'r-o');
    grid on;
    legend('Actual','Predicted','Location','northeastoutside');
    title('33,1 Autocorrelation - horizontal center','Interpreter','Latex', 'FontSize',14)
    ylabel('a (r)'); xlabel('r (cm)')
    xlim([0 max(rad_x)]); ylim([-1 1]);

    subplot(2,2,3)
    plot(rad_y,obj.a_u_11_3(Nx/2,:),'b-*',rad_y,exp_11_3(Nx/2,:),'r-o');
    grid on;
    legend('Actual','Predicted','Location','northeastoutside');
    title('11,3 Autocorrelation - vertical center','Interpreter','Latex', 'FontSize',14)
    ylabel('a (r)'); xlabel('r (cm)')
    xlim([0 max(rad_y)]); ylim([-1 1]);

    subplot(2,2,4)
    plot(rad_y,obj.a_v_33_3(Nx/2,:),'b-*',rad_y,exp_33_3(Nx/2,:),'r-o');
    grid on;
    legend('Actual','Predicted','Location','northeastoutside');
    title('33,3 Autocorrelation - vertical center','Interpreter','Latex', 'FontSize',14)
    ylabel('a (r)'); xlabel('r (cm)')
    xlim([0 max(rad_y)]); ylim([-1 1]);
    
    % integral length scale plots
    figure(2)
    subplot(2,2,1)
    plot(L_11_1,heights);
    grid on; 
    title('$L_{11,1}$','Interpreter','Latex', 'FontSize',14)
    ylabel('height (cm)'); xlabel('L (cm)')
    ylim([min(heights) max(heights)]);

    subplot(2,2,2)
    plot(L_33_1,heights);
    grid on;
    title('$L_{33,1}$','Interpreter','Latex', 'FontSize',14)
    ylabel('height (cm)'); xlabel('L (cm)')
    ylim([min(heights) max(heights)]);

    subplot(2,2,3)
    plot(widths, L_11_3);
    grid on;
    title('$L_{11,3}$','Interpreter','Latex', 'FontSize',14)
    ylabel('L (cm)'); xlabel('width (cm)')
    xlim([min(widths) max(widths)]);

    subplot(2,2,4)
    plot(widths, L_33_3);
    grid on;
    title('$L_{33,3}$','Interpreter','Latex', 'FontSize',14)
    ylabel('L (cm)'); xlabel('width (cm)')
    xlim([min(widths) max(widths)]);
            
    %%
    Resid=YhatG-a_v_33_1(lowb:uppb); % Calc. residuals
    SSR=sum(Resid.^2,2); % Sum of Squared Residuals
    SSE=sum((a_v_33_1(lowb:uppb)-mean(a_v_33_1(lowb:uppb))).^2); % Sum of Squared Errors
    RsquaredG(kp,row)=1-(SSR/SSE);
