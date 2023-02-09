classdef PIVanalysis < handle
    properties % test
        % list all variables that will be passed between methods in the
        % model here. think of them as global variables if a
        % value is included in the definition then it is treated as the
        % default value.

        % static model parameters with default values
  
        %g = 9.80665                    % (m/s2)                           
        mu = 2.5*1.81e-5                % (kg/m s) dynamic viscosity, u
        nu = 0.010034;                  % (cm2/s) % https://www.omnicalculator.com/physics/water-viscosity
        %nu =0.0000010034;              % m/s2 
        subwindow = 16;                 %subwindow size in pixels
        yaxis
        xaxis
        calibration
        heights
        widths

        u_original
        w_original
        datamax
        datamin
        original_percentleft
        agwpercentleft
        medianfilter_percentleft
        interpolated_percentleft
        m1             % mean flow strength in x direction
        m1_avg         % spatial average
        m1_median
        m3             % mean flow strength in z direction
        m3_avg         % spatial average
        m3_median
        mstar          % relative mean flow strength
        mstar_avg      % spatial average of relative mean flow strength
        mstar_median
        u_agw
        w_agw
        u_nanfilter
        w_nanfilter
        u_interpolated
        w_interpolated
        u_mean
        u_mean_savg
        u_mean_median
        w_mean
        w_mean_savg
        w_mean_median
        ugrad           % u velocity gradient

        u_f
        w_f
        u_rms
        w_rms
        u_rms_savg
        w_rms_savg
        u_rms_median
        w_rms_median

        tke             % turbulent kinetic energy
        tke_avg
        tke_median

        isotropy
        isotropy_avg    % spatial average
        isotropy_median % spatial median

        %integral scale values
        integral_avg    % integral length scale spatial average
        tau_integral_ls % integral time scale
        a_u_11_1
        a_w_33_1
        a_u_11_3
        a_w_33_3
        g_11_1
        g_33_3
        g_33_1
        g_11_3
        L_11_1
        L_33_1
        L_11_3
        L_33_3
        L_11_1_mean
        L_33_1_mean
        L_11_3_mean
        L_33_3_mean
        
        target
        lambda_tm_continuity       % Taylor microscale (centimeters) 
        lambda_tm_dvdy_dudx
        lambda_tm_dvdy_dwdz
        Re_lambda_continuity       % Taylor scale Reynolds number
        Re_lambda_dvdy_dudx
        Re_lambda_dvdy_dwdz

        % velocity gradients
        dudx
        dwdz
        dudz
        dwdx
        % gradient terms, squared and time averaged
        dudx_term
        dwdz_term
        dudz_term
        dwdx_term
        uxwz_term
        uzwx_term

        % dissipation values full equation, 3 different assumptions
        epsilon             % dissipation
        epsilon_avg         % dissipation spatial average
        tau_kt              % Kolmogorov time scale (sec) 
        eta_kl              % Kolmogorov length scale (time)
        epsilon_corrected
        epsilon_avg_corrected
        epsilon_median_corrected
        tau_kt_corrected
        eta_kl_corrected

        % dissipation values full equation, isotropy assumption 1
        epsilon_dvdy_dudx
        epsilon_dvdy_dudx_avg
        tau_kt_dvdy_dudx
        eta_kl_dvdy_dudx

        epsilon_dvdy_dudx_corrected
        epsilon_avg_dvdy_dudx_corrected 
        epsilon_median_dvdy_dudx_corrected
        tau_kt_dvdy_dudx_corrected 
        eta_kl_dvdy_dudx_corrected 
        
        % dissipation values full equation, isotropy assumption 2
        epsilon_dvdy_dwdz
        epsilon_dvdy_dwdz_avg
        tau_kt_dvdy_dwdz
        eta_kl_dvdy_dwdz

        epsilon_dvdy_dwdz_corrected
        epsilon_avg_dvdy_dwdz_corrected 
        epsilon_median_dvdy_dwdz_corrected
        tau_kt_dvdy_dwdz_corrected 
        eta_kl_dvdy_dwdz_corrected 
                          
    end      
    methods 
        % methods are defined to operate on class properties. the entire
        % framework for the numerical simulation will be constructed here.
        function obj = PIVanalysis()
            % include any initializations for properties here. It will be
            % ran whenever a new class instantiation is performed.

            %load('G:\J 20 All jets 2021\4_5 volts 1 sec 60 3 ms (new)\PIVlab_results','u_original','w_original')
            
            %%new tests Nov 2022 
            %load('G:\J20 Tests\4_5V_1sec_60percent_2_5ms_withlid\PIVlab_data','u_original','w_original')
            %load('/Users/almccutc/Desktop/PIVlab_results5V_1sec15p','u_original','v_original','calxy')
            load('/Users/almccutc/Desktop/PIVlab_results7V_1sec_15p','u_original','v_original','calxy')
            %load('PIVlab_results','u_original','v_original','calxy')

            obj.u_original = u_original; %(m/s)
            obj.w_original = v_original; %(m/s)
            obj.calibration = calxy*100; % Ex: 4e-05 cm/px
            
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        function allfunctions(obj)
            reshape(obj); % %this works well for PIVlab sessions
            checkHistogram(obj);
            applyAGWfilter(obj);
            medfilter(obj);
            velocityCalculations(obj);
            velocityPlots(obj)
            integrallength(obj)
            dissipation(obj);
            %delaunyinterpolation(obj); 
            %spatialspectra(obj);
            %temporalspectra(obj);
            
        end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        function reshape(obj)
         obj.u_original = cell2mat(permute(obj.u_original',[1,3,2])).*100; 
         obj.w_original = cell2mat(permute(obj.w_original',[1,3,2])).*100;
% 
         %uncomment for using agw and nan filters
         obj.u_original = obj.u_original(:,2:249,:); 
         obj.w_original = obj.w_original(:,2:249,:);

         %uncomment for using agw and nan filters
         obj.u_original(isnan(obj.u_original)) = NaN;
         obj.w_original(isnan(obj.w_original)) = NaN;
% 
%          obj.u_nanfilter = obj.u_original(:,1:248,1:5); 
%          obj.w_nanfilter = obj.w_original(:,1:248,1:5);
%          
% %          Matches NaN values for both u and w (i.e. if u has a NaN value at
% %          (1,1,1) and v does not, these lines will assign a NaN value at
% %          (1,1,1) for v to match u. 
%          obj.u_nanfilter(isnan(obj.u_nanfilter)) = NaN;
%          obj.w_nanfilter(isnan(obj.w_nanfilter)) = NaN;

        end
        function checkHistogram(obj)
            [Nt,Ny,Nx] = size(obj.u_original);
            total = Ny*Nx*Nt; 
            logicarray = ~isnan(obj.u_original);
            sumofnonNaN = sum(sum(sum(logicarray)));
            
            obj.original_percentleft = (sumofnonNaN/total)*100;

            ucheck = reshape(obj.u_original,1,Nx*Ny*Nt); wcheck = reshape(obj.w_original,1,Nx*Ny*Nt);

            figure(1)
            title('Histograms of the $u$ and w velocities - Pre-Filter')
            subplot(2,1,1)
            histogram(ucheck,100)
            title('Histogram of the $u$ velocities - Pre-Filter','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('cm/s')
            subplot(2,1,2)
            histogram(wcheck,100)
            title('Histogram of the $w$ velocities - Pre-Filter','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('cm/s')
            
            obj.datamax = str2num(cell2mat(inputdlg('Enter a number:',...
             'AGW Filter Max. and Min. Value', [1 50])));
            obj.datamin=-obj.datamax;
            close all
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function applyAGWfilter(obj)

            obj.u_original=permute(obj.u_original,[3 1 2]); obj.w_original=permute(obj.w_original,[3 1 2]);

            [Nt,Nx,Ny] = size(obj.u_original);

            ucheck = reshape(obj.u_original,1,Nx*Ny*Nt); wcheck = reshape(obj.w_original,1,Nx*Ny*Nt);

            uresh = reshape(obj.u_original,Nt*Nx,Ny); vresh = reshape(obj.w_original,Nt*Nx,Ny);
            %
                time=1:Nt*Nx;

                ufill = zeros(Ny,Nt*Nx); %creates double array with the same number of positions as the original data
                vfill = zeros(Ny,Nt*Nx);

                ufillvect=5; %not sure what this does or why it is 5

                for i = 1:Ny

                    utemp = uresh(:,i);
                    wtemp = vresh(:,i);
                    [udat utim] = agw_filter(utemp,time,obj.datamax,obj.datamin);
                    [vdat wtim] = agw_filter(wtemp,time,obj.datamax,obj.datamin);


                    %NaNs become 1000 and removed points become 1000.
                    %This uses utim/wtim to put the udat/vdat in the correct index of
                    %the original data
                    clear ufillvect
                    ufillvect = zeros(1,Nx*Nt);
                    ufillvect(utim)=udat;
                    ufill(i,:)=ufillvect; %the data that was not filtered, udat/vdat, is put into ufill/vfill
                    %keep their same index, where all other positions, which has data
                    %that filtered out, has a 1000 value

                    clear vfillvect
                    vfillvect = zeros(1,Nx*Nt);
                    vfillvect(wtim)=vdat;
                    vfill(i,:)=vfillvect;

                end

            %
            ufill = ufill'; obj.u_agw = reshape(ufill,Nt,Nx,Ny);
            vfill = vfill'; obj.w_agw = reshape(vfill,Nt,Nx,Ny);
            
            obj.u_agw(obj.u_agw==0)=NaN;
            obj.w_agw(obj.w_agw==0)=NaN;
            

            obj.u_agw(isnan(obj.w_agw)) = NaN;
            obj.w_agw(isnan(obj.u_agw)) = NaN;
            
            [Nt,Ny,Nx] = size(obj.w_agw);
            uaftercheck = reshape(obj.w_agw,1,Nx*Ny*Nt); waftercheck = reshape(obj.w_agw,1,Nx*Ny*Nt);

            % Percentage of velocities left after the AGW filter has been applied,
            % should be 90 percent or higher if the data is good
            [Ny,Nx,Nt] = size(obj.u_original);
            total = Ny*Nx*Nt; 
            logicarray = ~isnan(obj.u_agw);
            sumofnonNaN = sum(sum(sum(logicarray)));
            
            obj.agwpercentleft = (sumofnonNaN/total)*100; %it is the same for all components
            % Percent Remaining: Ideally 95% or more.')

            %----------------------
            % Histograms to compare pre and post filtering
            %----------------------
            figure(1)
            %title('Histograms of the $u$ and $w$ velocities---Pre and Post-Filter')
            subplot(2,2,1)
            
            histogram(ucheck,100)
            title('Histogram of the $u$ velocities','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('$u$ (cm/s)','Interpreter','Latex')
            annotation('textbox',[0 .8 .1 .2], ...
    'String','Pre-Filter    ','EdgeColor','none','Interpreter','Latex','fontsize', 18)
            annotation('textbox',[0 .3 .1 .2], ...
    'String','Post-Filter    ','EdgeColor','none','Interpreter','Latex','fontsize', 18)
            subplot(2,2,3)
            histogram(uaftercheck,100)
            title('Histogram of the $u$ velocities','Interpreter','Latex')
            ylabel('Frequency') 
            xlabel('$u$ (cm/s)','Interpreter','Latex')
            subplot(2,2,2)
            histogram(wcheck,100)
            title('Histogram of the $w$ velocities','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('$w$ (cm/s)','Interpreter','Latex')
            subplot(2,2,4)
            h = histogram(waftercheck,100);
            %set(h,'XData', obj.datamin:0.1:obj.datamax, 'YData',0:5000:max(h.Values))
%             xlim([obj.datamin-0.2 obj.datamax])
%             xticks(obj.datamin:0.1:obj.datamax)
            title('Histogram of the $w$ velocities','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('$w$ (cm/s)','Interpreter','Latex')
            set(gcf,'Position',[700 300 800 700])
%             ax = gca;
%             ax.YAxis.Exponent = 3;
            hold off
            
            pause;
            
        end
        function onlymedfilter(obj) % for K PIVlab data
            obj.target = 5; %initial guess 
            
            obj.u_original = obj.u_original(35:186, 1:161, :); %to do the needed reshape
            obj.w_original = obj.w_original(35:186, 1:161, :);
            
           %[Nt,Ny,Nx] = size(obj.u_original);
           [Ny,Nx,Nt] = size(obj.u_original);

            uclean_save = zeros(Ny,Nx,Nt);
            wclean_save = zeros(Ny,Nx,Nt);
           
            m=1;
            while(1)
            %selecttt = 0;
            if m==5
                ttstart = selecttt; 
            else 
                m=1;
                ttstart=1;
            end
            for tt=ttstart:Nt
                close all
                tt
                
                ut=squeeze(flip(obj.u_original(:,:,tt)));
                wt=squeeze(flip(obj.w_original(:,:,tt)));

                umed = mediannan(ut,3);
                vmed = mediannan(wt,3);

                flagu = abs(umed - ut) > obj.target;
                flagv = abs(vmed - wt) > obj.target;
                flag = flagu + flagv;
                flag(flag==2)=1;

                utclean = ut.*(1-flag);
                utclean(utclean==0) = NaN;
                wtclean = wt.*(1-flag);
                wtclean(wtclean==0) = NaN;
                
                uclean_save(:,:,tt) = utclean;
                wclean_save(:,:,tt) = wtclean;
                

                figure(2); %shows original data
                uoriginal = squeeze(flip(obj.u_original(:,:,tt)));
                woriginal = squeeze(flip(obj.w_original(:,:,tt)));
                
                scale_factor = 0.1;
                h1=quiver(uoriginal*scale_factor.*(flip(obj.A)),woriginal*scale_factor.*(flip(obj.A)),'r','AutoScale','off');
                xlim([0 Nx])
                ylim([0 Ny])
                set(gcf,'Position',[700 300 800 700])
                hold on;

                %the values after medfilter
                h3 = quiver(utclean*scale_factor.*(flip(obj.A)),wtclean*scale_factor.*(flip(obj.A)),'k','AutoScale','off'); 
                xlim([0 Nx])
                ylim([0 Ny])
                title(['Time step: ', num2str(tt), ' out of ', num2str(Nt), '. Target: ', num2str(obj.target)])
                legend('Original','AfterMedian','Location','northeastoutside');
                
                hold off
                scale=3;
                hU1 = get(h1,'UData');
                hV1 = get(h1,'WData');
                set(h1,'UData',scale*hU1,'WData',scale*hV1)
                hU3 = get(h3,'UData');
                hV3 = get(h3,'WData');
                set(h3,'UData',scale*hU3,'WData',scale*hV3)                
                
                if m==3
                    ttstart=1;
                    continue    
                else    
                    m = menu(['Yes if to continue through time, No for new' ...
                        ' target value (originally 5), All to apply' ...
                        ' filter at all time steps (reset to T.S. =1 first), ' ...
                        'Exit to stop program., Skip to TimeStep'],'Yes','No', 'All', 'Exit', 'Select T.S.');
                end 
                
                if m==2  % yes stored as 1, no stored as 2, all has a value of 3, exit is 4
                    break;
                end
                
                if m==4
                    break;
                end    
                
                if m==5
%                      selecttt = str2num(cell2mat(inputdlg('Enter new time step:',...
%             'Skip to this time step', [1 50])));
%                     m=3;
                    break;
                end 

                
            end
            if m==4
                break;
            end 
            if m==3
                break;
            end 
            if m==5
                selecttt = str2num(cell2mat(inputdlg('Enter new time step:',...
             'Skip to this time step', [1 50])));
            else
                obj.target = str2num(cell2mat(inputdlg('Enter new target:',...
            'Target', [1 50])));
            end 
 
            %check if the target should be updated then, depending on how
            %the plots looking 
            
            end
            obj.u_nanfilter = uclean_save;
            obj.w_nanfilter = wclean_save;
            
            obj.u_nanfilter = obj.u_nanfilter.*(flip(obj.A));
            obj.w_nanfilter = obj.w_nanfilter.*(flip(obj.A));
            
            %checks for percent remaining from the original data
            [Ny,Nx,Nt] = size(obj.u_original);
            total = Ny*Nx*Nt; 
            logicarray = ~isnan(obj.u_nanfilter);
            sumofnonNaN = sum(sum(sum(logicarray)));
            obj.medianfilter_percentleft = (sumofnonNaN/total)*100; 
            
            close all
            
            beep
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function medfilter(obj)
            obj.target = 8; %initial filter value guess 
            
            [Nt,Nx,Ny] = size(obj.u_agw);
            uclean_save = zeros(Nt,Nx,Ny); %preallocate 
            wclean_save = zeros(Nt,Nx,Ny); %preallocate 

            m=1;
            while(1)
            if m==5
                ttstart = selecttt; 
            else
                m=1;
                ttstart=1;
            end
            for tt=ttstart:Nt
                close all
                tt
                ut=squeeze(obj.u_agw(tt,:,:));
                wt=squeeze(obj.w_agw(tt,:,:));

                umed = mediannan(ut,3);
                vmed = mediannan(wt,3);

                flagu = abs(umed - ut) > obj.target;
                flagv = abs(vmed - wt) > obj.target;
                flag = flagu + flagv;
                flag(flag==2)=1;

                utclean = ut.*(1-flag);
                utclean(utclean==0) = NaN;
                wtclean = wt.*(1-flag);
                wtclean(wtclean==0) = NaN;

                uclean_save(tt,:,:) = utclean;
                wclean_save(tt,:,:) = wtclean;

                figure(2); %shows original data
                uoriginal = squeeze(obj.u_original(tt,:,:));
                woriginal = squeeze(obj.w_original(tt,:,:));
                scale_factor = 0.1;
                h1=quiver(uoriginal*scale_factor,woriginal*scale_factor,'b','AutoScale','off');
                xlim([0 Ny])
                ylim([0 Nx])
                set(gcf,'Position',[700 100 900 700])
                hold on;
                
                %shows what agw filter removes
                uagw = squeeze(obj.u_agw(tt,:,:));
                vagw = squeeze(obj.w_agw(tt,:,:));
                h2= quiver(uagw*scale_factor,vagw*scale_factor,'r','AutoScale','off');
                xlim([0 Ny])
                ylim([0 Nx])
                hold on;
                
                %the values after agw and medfilter
                h3 = quiver(utclean*scale_factor,wtclean*scale_factor,'k','AutoScale','off'); 
                xlim([0 Ny])
                ylim([0 Nx])
                title(['Time step: ', num2str(tt), ' out of ', num2str(Nt), '. Target: ', num2str(obj.target)])
                legend('AGW','Median','Remaining','Location','northeastoutside');
                
                hold off
                scale=1;
                hU1 = get(h1,'UData');
                hV1 = get(h1,'WData');
                set(h1,'UData',scale*hU1,'WData',scale*hV1)
                hU2 = get(h2,'UData');
                hV2 = get(h2,'WData');
                set(h2,'UData',scale*hU2,'WData',scale*hV2)
                hU3 = get(h3,'UData');
                hV3 = get(h3,'WData');
                set(h3,'UData',scale*hU3,'WData',scale*hV3)                
                
                if m==3
                    continue    
                else    
                    m = menu(['Yes if to continue through time, No for new target ' ...
                        'value (originally 0.15), All to apply filter at all ' ...
                        'time steps, Exit to stop program., Skip to TimeStep'], ...
                        'Yes','No', 'All', 'Exit', 'Select T.S.');
                end 
                
                if m==2  % yes stored as 1, no stored as 2, all has a value of 3, exit is 4
                    break;
                end
                
                if m==4
                    break;
                end    
                
                if m==5
                    break;
                end 
            end

            if m==4
                break;
            end 
            if m==3
                break;
            end 
            if m==5
                selecttt = str2double(cell2mat(inputdlg('Enter new time step:',...
             'Skip to this time step', [1 50])));
            else
                obj.target = str2double(cell2mat(inputdlg('Enter new target:',...
            'Target', [1 50])));
            end 
            
            %check if the target should be updated depending on how
            %the plot looks
            
            end
            obj.u_nanfilter = permute(uclean_save,[2 3 1]);
            obj.w_nanfilter = permute(wclean_save,[2 3 1]);

            %checks for percent remaining from the original data
            [Ny,Nx,Nt] = size(obj.u_original);
            total = Ny*Nx*Nt; 
            logicarray = ~isnan(obj.u_nanfilter);
            sumofnonNaN = sum(sum(sum(logicarray)));
            obj.medianfilter_percentleft = (sumofnonNaN/total)*100; 
            
            close all
            
            beep
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function bootstrap(obj)
            bs_ufull=sort(bootstrp(1000,'mean',obj.u_nanfilter));
            bs_ulower =  bs_ufull(25); bs_uupper =  bs_ufull(975);
            
            bs_wfull=sort(bootstrp(1000,'mean',obj.w_nanfilter));
            bs_wlower =  bs_wfull(25); bs_wupper =  bs_wfull(975);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function delaunyinterpolation(obj) 
            
            unan = obj.u_nanfilter;
            vnan = obj.w_nanfilter;
            
            [Nt, Ny, Nx] = size(unan);
            
            %replace any corner NaNs with the value of the median of existing corners
            %so that interpolation of inner parts can take place

            %upper left
            U_UL=nanmedian(unan(:,1,1));
            W_UL=nanmedian(vnan(:,1,1));
            %upper right
            U_UR=nanmedian(unan(:,1,Nx));
            W_UR=nanmedian(vnan(:,1,Nx));
            %lower left
            U_LL=nanmedian(unan(:,Ny,1));
            W_LL=nanmedian(vnan(:,Ny,1));
            %lower right
            U_LR=nanmedian(unan(:,Ny,Nx));
            W_LR=nanmedian(vnan(:,Ny,Nx));
            
            
            for tt=1:Nt
                tt
                unan_temp=squeeze(unan(tt,:,:));
                vnan_temp=squeeze(vnan(tt,:,:));
                [uuu,vvv] = fillNaNsTemp(unan_temp,vnan_temp, U_UL, U_UR, U_LL, U_LR, W_UL, W_UR, W_LL, W_LR);
                u_interp(tt,:,:)=uuu;
                w_interp(tt,:,:)=vvv;

                close all
            end      
            
            obj.u_interpolated = u_interp(:,:,:);
            obj.w_interpolated = w_interp(:,:,:);    
            
            [Ny,Nx,Nt] = size(obj.u_original);
            total = Ny*Nx*Nt; 
            logicarray = ~isnan(obj.u_interpolated);
            sumofnonNaN = sum(sum(sum(logicarray)));
            obj.interpolated_percentleft = (sumofnonNaN/total)*100;
            
            beep
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function velocityCalculations(obj) % and tke, isotropy, and mean flow strength
            u_o = obj.u_nanfilter; w_o = obj.w_nanfilter;
            
            [Ny,Nx,Nt] = size(u_o);
            
            obj.u_mean = mean(u_o,3,'omitnan'); obj.w_mean = mean(w_o,3,'omitnan');
            obj.u_mean_savg = mean(obj.u_mean,'all','omitnan'); obj.w_mean_savg = mean(obj.w_mean,'all','omitnan');
            obj.u_mean_median = median(obj.u_mean,'all','omitnan'); obj.w_mean_median = median(obj.w_mean,'all','omitnan');

            obj.u_f = u_o-obj.u_mean; obj.w_f = w_o-obj.w_mean;

            obj.u_rms = sqrt(mean(obj.u_f.^2,3,'omitnan')); obj.w_rms = sqrt(mean(obj.w_f.^2,3,'omitnan'));
            obj.u_rms_savg = mean(obj.u_rms,'all','omitnan'); obj.w_rms_savg = mean(obj.w_rms,'all','omitnan');
            obj.u_rms_median = median(obj.u_rms,'all','omitnan'); obj.w_rms_median = median(obj.w_rms,'all','omitnan');

            tke = 0.5*(2*(obj. u_rms.^2) + (obj. w_rms.^2));
            obj.tke = mean(tke,3,'omitnan');
            obj.tke_avg = mean(obj.tke,'all','omitnan'); 
            obj.tke_median = median(obj.tke,'all','omitnan'); 

            obj.isotropy = obj.u_rms./obj.w_rms; 
            obj.isotropy_avg = mean(obj.isotropy,'all','omitnan'); 
            obj.isotropy_median = median(obj.isotropy,'all','omitnan'); 

            obj.yaxis = (-Ny/2:1:Ny/2)*obj.calibration*obj.subwindow; %converts subwindow count to cm
            obj.xaxis = (-Nx/2:1:Nx/2)*obj.calibration*obj.subwindow; %converts subwindow count to cm

            %mean flow strength 
            obj.m1 = obj.u_mean./obj.u_rms; 
            obj.m1_avg = mean(obj.m1,'all','omitnan');
            obj.m1_median = median(obj.m1,'all','omitnan');

            obj.m3 = obj.w_mean./obj.w_rms;
            obj.m3_avg = mean(obj.m3,'all','omitnan');
            obj.m3_median = median(obj.m3,'all','omitnan');
            
            obj.mstar = (0.5*(2*(obj.u_mean.^2) + (obj.w_mean.^2)))./(obj.tke); %time average U, W, and tke
            obj.mstar_avg = mean(obj.mstar,'all','omitnan');
            obj.mstar_median = median(obj.mstar,'all','omitnan');

        end
        function velocityPlots(obj)     
            %Time average for each subwindow                            
            figure (1) 
            subplot(2,2,1)
            imagesc((obj.u_mean), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal') 
            title(['$\overline{U}:\:',(num2str(obj.u_mean_savg,3)),'\:cm/s\;\;M_d:\:',(num2str(obj.u_mean_median,3)),'\:cm/s$'],'Interpreter','latex', 'FontSize',14)
            ylabel('cm')
            xlabel('cm')

            subplot(2,2,2)
            imagesc((obj.w_mean), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            title(['$\overline{W}:\:',(num2str(obj.w_mean_savg,3)),'\:cm/s\;\;M_d:\:',(num2str(obj.w_mean_median,3)),'\:cm/s$'],'Interpreter','latex', 'FontSize',14)
            ylabel('cm')
            xlabel('cm')

            %RMS velocity for each subwindow
            subplot(2,2,3)
            imagesc((obj.u_rms), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            title(['$\overline{u_{rms}}:\:',(num2str(obj.u_rms_savg,3)),'\:cm/s\;\;M_d:\:',(num2str(obj.u_rms_median,3)),'\:cm/s$'],'Interpreter','latex', 'FontSize',14)
            ylabel('cm')
            xlabel('cm')

            subplot(2,2,4)
            imagesc((obj.w_rms), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            title(['$\overline{w_{rms}}:\:',(num2str(obj.w_rms_savg,3)),'\:cm/s\;\;M_d:\:',(num2str(obj.w_rms_median,3)),'\:cm/s$'],'Interpreter','latex', 'FontSize',14)
            ylabel('cm')
            xlabel('cm')
            set(gcf,'Position',[700 300 1000 700])

            figure (2)
            imagesc((obj.tke), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            title(['$\overline{k}:\:',(num2str(obj.tke_avg,3)),'\:cm^2/s^2\;\;M_d:\:',(num2str(obj.tke_median,3)),'\:cm^2/s^2$'],'Interpreter','latex', 'FontSize',14)
            ylabel('cm')
            xlabel('cm')
            
            figure (3)
            imagesc((obj.isotropy), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            title(['$Isotropy\:-\:\frac{u_{rms}}{w_{rms}}:\:',(num2str(obj.isotropy_avg,3)),'\;\;M_d:\:',(num2str(obj.isotropy_median,3)),'$'],'Interpreter','latex', 'FontSize',14)
            ylabel('cm')
            xlabel('cm')
            
            figure (4) 
            subplot(1,4,1)
            imagesc((obj.m1*100), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            title(['$\overline{M_1}:\:',(num2str(obj.m1_avg*100,3)),'\: \% \;\;M_d:\:',(num2str(obj.m1_median*100,3)),'\:\%$'],'Interpreter','latex', 'FontSize',14)
            ylabel('cm')
            xlabel('cm')

            subplot(1,4,2)
            imagesc((obj.m3*100), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')         
            title(['$\overline{M_3}:\:',(num2str(obj.m3_avg*100,3)),'\:\%\;\;M_d:\:',(num2str(obj.m3_median*100,3)),'\:\% $'],'Interpreter','latex', 'FontSize',14)
            ylabel('cm')
            xlabel('cm')
            
            subplot(1,4,3)
            imagesc((obj.mstar*100), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')       
            title(['$\overline{M^*}:\:',(num2str(obj.mstar_avg*100,3)),'\:\%\;\;M_d:\:',(num2str(obj.mstar_median*100,3)),'\:\% $'],'Interpreter','latex', 'FontSize',14)
            ylabel('cm')
            xlabel('cm')
            set(gcf,'Position',[700 300 1240 300])

            subplot(1,4,4)
            imagesc((obj.mstar*100), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            caxis([0 10])
            set(gca,'YDir','normal')       
            title(['$M^*\:10 \%$ value limit'],'Interpreter','latex', 'FontSize',14)
            ylabel('cm')
            xlabel('cm')
            set(gcf,'Position',[700 300 1240 300])

        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function taylorscales(obj) % which average to use for the integral length scale

            % Taylor microscale (centimeters) 
            obj.lambda_tm_continuity = sqrt(10)*obj.eta_kl_corrected^(2/3)* obj.integral_avg^(1/3);
            obj.lambda_tm_dvdy_dudx = sqrt(10)*obj.eta_kl_dvdy_dudx_corrected^(2/3)* obj.integral_avg^(1/3);
            obj.lambda_tm_dvdy_dwdz = sqrt(10)*obj.eta_kl_dvdy_dwdz_corrected ^(2/3)* obj.integral_avg^(1/3);
            
            % Taylor scale Reynolds number
            obj.Re_lambda_continuity = (2/3)*obj.tke_avg*sqrt(15/(obj.nu*obj.epsilon_avg_corrected));
            obj.Re_lambda_dvdy_dudx  = (2/3)*obj.tke_avg*sqrt(15/(obj.nu*obj.epsilon_avg_dvdy_dudx_corrected));
            obj.Re_lambda_dvdy_dwdz = (2/3)*obj.tke_avg*sqrt(15/(obj.nu*obj.epsilon_avg_dvdy_dwdz_corrected));

        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function spatialspectra(obj)
            
            u_o = obj.u_interpolated; w_o = obj.w_interpolated;
            
            u_o=permute(u_o,[3 2 1]); w_o=permute(w_o,[3 2 1]);             
            
            %53 by 79 subwindows
            Suutime=[];
            a = 35; %x direction subwindow strip
            b= 25; %y direction subwindow strip
            %calibration = s; %pixel to cm conversion
            deltax=16; %cm% smallest subwindow overlap. 32 by 32 subwindow with a 50% overlap is 16 same as delta y
            deltay = deltax; % in this case with a symmetric subwindow
            Lx = 1280-deltax; %cm% distance from the middle of the first to last subwindow
            Ly = 864-deltax; %cm% distance from the middle of the first to last subwindow
            ks = 2*pi/deltax; %cm% 
            deltakx = 2*pi/Lx; %between frequency points
            deltaky = 2*pi/Ly;
            kx= deltakx/2:deltakx:(ks-deltakx/2); %x axis values
            ky= deltaky/2:deltaky:(ks-deltaky/2); %y axis values
            
            for i = 1:nlay-10000 %number of image pairs

                uspatial = u_o(a,:,i);
                vspatial = w_o(b,:,i);  
                Suu = fft(uspatial).*conj(fft(uspatial));
                Svv = fft(vspatial).*conj(fft(vspatial));
            end 
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function temporalspectra(obj) %this mostly works but there are issues with the number of time steps
            %right it is set to ensemble average of 10 runs
            u_o = obj.u_interpolated; w_o = obj.w_interpolated;
            
            u_o=permute(u_o,[3 2 1]); w_o=permute(w_o,[3 2 1]); 
            
            %u_o(isnan(u_o)) = 0; w_o(isnan(w_o)) = 0;

            %436(left to right), 288(top going dowwards)
            %so approx subwindow (18,27)
            xlocation = 18;
            ylocation = 27;

            u = u_o;
            v = w_o;
            [Nx,Ny,Nt] = size(u);
            f_s = 105; %sample frequency
            N = Nt/10; 
            %N=1150;
            T = N/f_s;
            df = f_s/N;
            f = df/2:df:f_s-df/2;
            fny = f(1:(length(f)/2));

            tempspectraU = [];
            tempspectraW = [];

            %this is done as an ensemble average
            for i =1:10 
               U = u(xlocation,ylocation,N*i-N+1:N*i);
               W = v(xlocation,ylocation,N*i-N+1:N*i);
               tempspectraU(:,:,i) = ((1/(f_s*T*f_s))*abs((fft(U)).^2)); 
               tempspectraW(:,:,i) = ((1/(f_s*T*f_s))*abs((fft(W)).^2)); 
            end

            Suu=mean(tempspectraU,3);
            Suu=permute(Suu, [2 1]); 
            Suu = 2.*Suu(1:length(fny));

            Sww=mean(tempspectraW,3);
            Sww=permute(Sww, [2 1]); 
            Sww = 2.*Sww(1:length(fny));

            figure(1)
            loglog(fny,Suu, 'g')
            xlabel('f (Hz)','FontSize',12)
            ylabel('S_{uu}, S_{ww} (m^{2}/s^{3})','FontSize',12)
            xlim([-inf inf]);
            ylim([10^-5 10^-2]);
            grid on;
            hold on
            loglog(fny,Sww, 'r')
            %title({'\fontsize{14} \it Autospectral Density Function';' '},'FontWeight','Normal','Interpreter','Latex')
             b=-8;
             c=exp(b);
            hold on
            y = c*fny.^(-5/3);
            loglog(fny,y, 'k')

            legend('S_{uu}','S_{ww}', '-5/3' )
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        function integrallength(obj)    

             u_o = permute(obj.u_f,[3 1 2]); w_o = permute(obj.w_f,[3 1 2]); %cm/s
            
            [~,Ny, Nx]=size(u_o); 
            
            x_c=(Nx+1)/2; % determines the centerline position
            y_c=(Ny+1)/2; % determines the centerline position
            rad_x =[0,(obj.subwindow:obj.subwindow*2:Nx*obj.subwindow).*obj.calibration]; %calibrate to cm
            rad_y =[0,(obj.subwindow:obj.subwindow*2:Ny*obj.subwindow).*obj.calibration]; %calibrate to cm
            obj.heights=([-Ny/2:1:-1 1:1:Ny/2]).*obj.calibration.*obj.subwindow; %calibrate to cm
            obj.widths=([-Nx/2:1:-1 1:1:Nx/2]).*obj.calibration.*obj.subwindow; %calibrate to cm
            
            obj.a_u_11_1 = [ones([Ny 1]) NaN([Ny Nx/2])]; obj.a_w_33_1 = [ones([Ny 1]) NaN([Ny Nx/2])];
            obj.a_u_11_3 = [ones([Nx 1]) NaN([Nx Ny/2])]; obj.a_w_33_3 = [ones([Nx 1]) NaN([Nx Ny/2])];
            
             
            %spatial, horizontal autocorrelation calculation, for the case of an even # of vertical and horizontal subwindows          
            for row=1:Ny % calculated at every height for each subwindow 
                for radius=0.5:1:x_c-1
            
                    %11,1 - longitudinal, Horizontal velocity, horizontal separation
                    obj.a_u_11_1(row,radius+1.5)=mean(u_o(:,row,x_c-radius).*u_o(:,row,x_c+radius),'omitnan')...
                        ./sqrt(mean(u_o(:,row,x_c-radius).^2,'omitnan').*mean(u_o(:,row,x_c+radius).^2,'omitnan'));
                    
                    %33,1 - transverse, Vertical velocity, horizontal separation
                    obj.a_w_33_1(row,radius+1.5)=mean(w_o(:,row,x_c-radius).*w_o(:,row,x_c+radius),'omitnan')...
                        ./sqrt(mean(w_o(:,row,x_c-radius).^2,'omitnan').*mean(w_o(:,row,x_c+radius).^2,'omitnan'));
                end
            end
            
            %spatial, vertical autocorrelation calculation, for the case of an even # of vertical and horizontal subwindows          
            for column=1:Nx % calculated at every width for each subwindow 
                for radius=0.5:1:y_c-1
            
                    %11,3 - transverse, Horizontal velocity, vertical separation
                    obj.a_u_11_3(column, radius+1.5)=mean(u_o(:,y_c-radius,column).*u_o(:,y_c+radius,column),'omitnan')...
                        ./sqrt(mean(u_o(:,y_c-radius,column).^2,'omitnan').*mean(u_o(:,y_c+radius,column).^2,'omitnan'));
                    
                    %33,3 - longitudinal, Vertical velocity, vertical separation
                    obj.a_w_33_3(column,radius+1.5)=mean(w_o(:,y_c-radius,column).*w_o(:,y_c+radius,column),'omitnan')...
                        ./sqrt(mean(w_o(:,y_c-radius,column).^2,'omitnan').*mean(w_o(:,y_c+radius,column).^2,'omitnan'));
                end
            end
            
            %calculate L at all heights
            for row=1:Ny
                
                % Length scale calculation
                newExpFuncG=@(rr,ll) exp(-rr./ll); 
                newExpFuncG_transverse=@(rr,ll) exp(-rr./ll).*(1-(rr./(2.*ll))); 
            
                L0=6; % Starting guess for L
            
                obj.g_11_1=@(ll)sum((obj.a_u_11_1(row,:)-newExpFuncG(rad_x,ll)).^2,2);
                obj.g_33_1=@(ll)sum((obj.a_w_33_1(row,:)-newExpFuncG_transverse(rad_x,ll)).^2,2);
            
                % Minimize to find Lstar
                obj.L_11_1(row,1)=fminunc(obj.g_11_1,L0); 
                obj.L_11_1_mean = mean(obj.L_11_1);
                L_33_1(row,1)=fminunc(obj.g_33_1,L0); obj.L_33_1 = L_33_1*.5;
                obj.L_33_1_mean = mean(obj.L_33_1);
            
                exp_11_1(row,:)=newExpFuncG(rad_x,obj.L_11_1(row,1)); 
                exp_33_1(row,:)=newExpFuncG_transverse(rad_x,obj.L_33_1(row,1)); 
            end    
            
            %calculate L at all widths
            for column=1:Nx

                % Length scale calculation
                newExpFuncG=@(rr,ll) exp(-rr./ll); 
                newExpFuncG_transverse=@(rr,ll) exp(-rr./ll).*(1-(rr./(2.*ll))); 

                L0=6; % Starting guess for L
            
                obj.g_11_3=@(ll)sum((obj.a_u_11_3(column,:)-newExpFuncG_transverse(rad_y,ll)).^2,2);
                obj.g_33_3=@(ll)sum((obj.a_w_33_3(column,:)-newExpFuncG(rad_y,ll)).^2,2);
            
                % Minimize to find Lstar
                L_11_3(column,1)=fminunc(obj.g_11_3,L0); obj.L_11_3 = L_11_3*.5;
                obj.L_11_3_mean = mean(obj.L_11_3);
                obj.L_33_3(column,1)=fminunc(obj.g_33_3,L0); 
                obj.L_33_3_mean = mean(obj.L_33_3);
               
            
                exp_11_3(column,:)=newExpFuncG_transverse(rad_y,obj.L_11_3(column,1)); 
                exp_33_3(column,:)=newExpFuncG(rad_y,obj.L_33_3(column,1)); 
            end   
            
            % autocorrelation plots
            figure(6)
            t = tiledlayout(2,2,'TileSpacing','Compact');
            title(t,'Autocorrelation values vs. r','Interpreter','Latex', 'FontSize',14)
            xlabel(t,'r (cm)'); ylabel(t,'a (r)')
            grid on;
            % Tile 1
            nexttile, plot(rad_x,obj.a_u_11_1(Ny/2,:),'k.',rad_x,exp_11_1(Ny/2,:),'k');
            title('$a_{11,1}(r)$ - horizontal center','Interpreter','Latex', 'FontSize',14)
            grid on;
            % Tile 2
            nexttile, plot(rad_x,obj.a_w_33_1(Ny/2,:),'k.',rad_x,exp_33_1(Ny/2,:),'k');
            title('$a_{33,1}(r)$ - horizontal center','Interpreter','Latex', 'FontSize',14)
            grid on;
            % Tile 3
            nexttile, plot(rad_y,obj.a_u_11_3(Nx/2,:),'k.',rad_y,exp_11_3(Nx/2,:),'k');
            title('$a_{11,3}(r)$ - vertical center','Interpreter','Latex', 'FontSize',14)
            grid on;
            % Tile 4
            nexttile, plot(rad_y,obj.a_w_33_3(Nx/2,:),'k.',rad_y,exp_33_3(Nx/2,:),'k');
            title('$a_{33,3}(r)$ - vertical center','Interpreter','Latex', 'FontSize',14)
            grid on;
            lg  = legend('Actual','Predicted'); lg.Layout.Tile = 'East'; % <-- Legend placement with tiled layout
                
            % integral length scale plots
            figure(7)
            t = tiledlayout(2,2,'TileSpacing','Compact');
            title(t,'Integral Length Scale, $L$ (cm)','Interpreter','Latex', 'FontSize',14)
            % Tile 1
            nexttile
            plot(obj.L_11_1,obj.heights);
            title(['$\overline{L_{11,1}}\::\:',(num2str(obj.L_11_1_mean,3)),'\:cm$'],'Interpreter','latex', 'FontSize',14)
            ylabel('y (cm)'); xlabel('L (cm)'); grid on;
            ylim([min(obj.heights) max(obj.heights)]);
            % Tile 2
            nexttile
            plot(obj.L_33_1,obj.heights);
            title(['$\overline{L_{33,1}}\::\:',(num2str(obj.L_33_1_mean,3)),'\:cm$'],'Interpreter','latex', 'FontSize',14)
            ylabel('y (cm)'); xlabel('L (cm)'); grid on;
            ylim([min(obj.heights) max(obj.heights)]);
            % Tile 3
            nexttile
            plot(obj.widths, obj.L_11_3);
            title(['$\overline{L_{11,3}}\::\:',(num2str(obj.L_11_3_mean,3)),'\:cm$'],'Interpreter','latex', 'FontSize',14)
            %
            ylabel('L (cm)'); xlabel('x (cm)'); grid on;
            xlim([min(obj.widths) max(obj.widths)]);
            % Tile 4
            nexttile
            plot(obj.widths, obj.L_33_3);
            title(['$\overline{L_{33,3}}\::\:',(num2str(obj.L_33_3_mean,3)),'\:cm$'],'Interpreter','latex', 'FontSize',14)
            ylabel('L (cm)'); xlabel('x (cm)'); grid on;
            xlim([min(obj.widths) max(obj.widths)]);             

            %obj.tau_intergral_ts = obj.integral_avg/sqrt(obj.tke_avg);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        function dissipation(obj)  

            %<10 cm2/s3
            %time, 0.01 s?

            deltax = obj.subwindow*obj.calibration; %distance between velocity vectors (cm)
            
            [Ny,Nx,Nt] = size(obj.u_f); %Ny-# of rows, Nx-# of columns
            obj.yaxis = (-Ny/2:20:Ny/2)*obj.calibration*obj.subwindow; %converts subwindow count to cm
            obj.xaxis = (-Nx/2:20:Nx/2)*obj.calibration*obj.subwindow; %converts subwindow count to cm
            obj.heights = ([-Ny/2:1:-1 1:1:Ny/2]).*obj.calibration.*obj.subwindow; %calibrate to cm
            
            %zero matrices
            z_u = zeros(Ny,1,Nt); 
            z_v = zeros(1,Nx,Nt);
            
            %size adjusted matrices
            u2adjx = [obj.u_f z_u];
            w2adjz = [obj.w_f; z_v];
            u2adjz = [obj.u_f; z_v];
            w2adjx = [obj.w_f z_u];
            u1adjx = [z_u obj.u_f];
            w1adjz = [z_v; obj.w_f];
            u1adjz = [z_v; obj.u_f];
            w1adjx = [z_u obj.w_f];
            
            %velocity gradients
            obj.dudx = u2adjx - u1adjx;
            obj.dwdz = w2adjz - w1adjz;
            obj.dudz = u2adjz - u1adjz;
            obj.dwdx = w2adjx - w1adjx;
            
            %gradient terms, squared and time averad
            obj.dudx_term = obj.nu.*mean(((obj.dudx(2:Ny,2:Nx,:)./deltax).^2),3,'omitnan'); %time average
            obj.dwdz_term = obj.nu.*mean(((obj.dwdz(2:Ny,2:Nx,:)./deltax).^2),3,'omitnan'); %time average
            obj.dudz_term = obj.nu.*mean(((obj.dudz(2:Ny,2:Nx,:)./deltax).^2),3,'omitnan'); %time average
            obj.dwdx_term = obj.nu.*mean(((obj.dwdx(2:Ny,2:Nx,:)./deltax).^2),3,'omitnan'); %time average
            obj.uxwz_term = obj.nu.*mean((((obj.dudx(2:Ny,2:Nx,:)/deltax).*(obj.dwdz(2:Ny,2:Nx,:)./deltax))),3,'omitnan'); %time average
            obj.uzwx_term = obj.nu.*mean((((obj.dudz(2:Ny,2:Nx,:)./deltax).*(obj.dwdx(2:Ny,2:Nx,:)./deltax))),3,'omitnan'); %time average

            figure (8)
            plot((mean(obj.dudx_term,2,'omitnan')),obj.heights(1:Ny-1),'-r',...
                (mean(obj.dwdz_term,2,'omitnan')),obj.heights(1:Ny-1),'--r',...
                (mean(obj.dudz_term,2,'omitnan')),obj.heights(1:Ny-1),'-k',...
                (mean(obj.dwdx_term,2,'omitnan')),obj.heights(1:Ny-1),'--k',...
                (mean(obj.uxwz_term,2,'omitnan')),obj.heights(1:Ny-1),'-b',(mean...
                (obj.uzwx_term,2,'omitnan')),obj.heights(1:Ny-1),'-g');
            title('Dissipation Components','Interpreter','latex', 'FontSize',16)
            legend('$\overline{(\frac{\partial u}{\partial x})^2}$','$\overline{(\frac{\partial w}{\partial z})^2}$','$\overline{(\frac{\partial u}{\partial z})^2}$','$\overline{(\frac{\partial w}{\partial x})^2}$',...
                '$\overline{(\frac{\partial u}{\partial x} \frac{\partial w}{\partial z})}$','$\overline{(\frac{\partial u}{\partial x} \frac{\partial w}{\partial z})}$','Interpreter','Latex', 'FontSize',14)
            ylabel('y (cm)'); xlabel('cm^2/s^3'); grid on;
            
            %%%%%%%%%% continuity assumption
            %%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % dissipation rate time average, m2/s3
            obj.epsilon =2.*(4.*obj.dudx_term+obj.dudz_term+...
                obj.dwdx_term+2.*obj.dwdz_term+2.*obj.uxwz_term+2.*obj.uzwx_term);
            % dissipation rate spatial average
            obj.epsilon_avg = mean(obj.epsilon,'all','omitnan');   
            obj.tau_kt = (obj.nu/obj.epsilon_avg)^0.5; % time (s)
            obj.eta_kl = (obj.nu^3/obj.epsilon_avg)^0.25; %length (m)
            % integrated dissipation spectrum
            R = (deltax/obj.eta_kl); 
            x_pos = 2*pi/R;

            % figure that shows the Intergration of Universal Spectrum
            % enter the correction value
            f = figure('Units','Normalized',...
                 'Position',[.1 .1 .05 .05],...
                 'NumberTitle','off',...
                 'Name','Dissipation Correction - Continuity assumption');
            imshow('universalspectrum.png')
            title('Intergration of Universal Spectrum, Enter % value.');
            e = uicontrol('Style','Edit',...
                 'Units','Normalized',...
                 'Position',[.2 .4 .1 .1],...
                 'Tag','myedit');
            p = uicontrol('Style','PushButton',...
                 'Units','Normalized',...
                 'Position',[.3 .4 .1 .1],...
                 'String','Enter',...
                 'CallBack','uiresume(gcbf)');
            annotation('textbox',[.2 .5 .4 .1],'String',['2*\pi/R: '...
                ,num2str(x_pos,3)],'FitBoxToText','on','BackgroundColor','w');

            uiwait(f)
            cvalue = str2double(get(e,'String'));
            close all           
            
            %corrected values
            obj.epsilon_corrected = obj.epsilon.*(2-cvalue/100);
            obj.epsilon_avg_corrected = mean(obj.epsilon_corrected,'all','omitnan'); 
            obj.epsilon_median_corrected = median(obj.epsilon_corrected,'all','omitnan'); 

            % Kolmogorov length and time scale
            obj.tau_kt_corrected = (obj.nu/obj.epsilon_avg_corrected)^0.5; % time (s)
            obj.eta_kl_corrected = (obj.nu^3/obj.epsilon_avg_corrected)^0.25; %length (cm)
            
%             % using x position to go up and across
%             figure(2)
%             plot(f,Gain, '-.b') % dissipation spectra
%             hold on
%             plot(f,G2, '--r') % integrated dissipation spectrum
%             xlabel('k_\eta')
%             ylabel('D(k)/u^3_\eta, \epsilon(0,k)/\epsilon_m')
%             title('Normalized dissipation spectrum and cumulative dissipation')
%             legend('Normalized dissipation spectrum','Cumulative dissipation')

            %%%%%%%%%% Isotropic assumption 1 - dv/dy = du/dx %%%%%%%%%% 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            obj.epsilon_dvdy_dudx = 2.*(4.*obj.dudx_term+obj.dudz_term+...
            obj.dwdx_term+obj.dwdz_term+2.*obj.uzwx_term);
            obj.epsilon_dvdy_dudx_avg = mean(obj.epsilon_dvdy_dudx,'all','omitnan'); % dissipation rate spatial average

            % Kolmogorov length and time scale
            obj.tau_kt_dvdy_dudx = (obj.nu/obj.epsilon_dvdy_dudx_avg)^0.5; % time (s)
            obj.eta_kl_dvdy_dudx = (obj.nu^3/obj.epsilon_dvdy_dudx_avg)^0.25; %length (m)
            
            %integrated dissipation spectrum
            R = (deltax/obj.eta_kl_dvdy_dudx); 
            x_pos = 2*pi/R;
            
%             %correction value
%             cvalue = str2num(cell2mat(inputdlg('Value on chart for the correction in %?',...
%              'Dissipation Correction', [1 50])));

             f = figure('Units','Normalized',...
                 'Position',[.1 .1 .05 .05],...
                 'NumberTitle','off',...
                 'Name','Dissipation Correction');
            imshow('universalspectrum.png')
            title('Interation of Universal Spectrum, Enter % value.');
            e = uicontrol('Style','Edit',...
                 'Units','Normalized',...
                 'Position',[.2 .4 .1 .1],...
                 'Tag','myedit');
            p = uicontrol('Style','PushButton',...
                 'Units','Normalized',...
                 'Position',[.3 .4 .1 .1],...
                 'String','Enter',...
                 'CallBack','uiresume(gcbf)');
            annotation('textbox',[.2 .5 .4 .1],'String',['2*\pi/R: ',num2str(x_pos,3)],'FitBoxToText','on','BackgroundColor','w');

            uiwait(f)
            cvalue = str2double(get(e,'String'));
            close all           
            
            %corrected values
            obj.epsilon_dvdy_dudx_corrected = obj.epsilon_dvdy_dudx*(2-cvalue/100);
            obj.epsilon_avg_dvdy_dudx_corrected = mean(obj.epsilon_dvdy_dudx_corrected,'all','omitnan');
            obj.epsilon_median_dvdy_dudx_corrected = median(obj.epsilon_dvdy_dudx_corrected,'all','omitnan'); 

            % Kolmogorov length and time scale
            obj.tau_kt_dvdy_dudx_corrected = (obj.nu/obj.epsilon_avg_dvdy_dudx_corrected)^0.5; % time
            obj.eta_kl_dvdy_dudx_corrected = (obj.nu^3/obj.epsilon_avg_dvdy_dudx_corrected)^0.25; %length (cm)
            
            %%%%%%%%%% Isotropic assumption 2 - dv/dy = dw/dz %%%%%%%%%% 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            obj.epsilon_dvdy_dwdz = 2.*(3.*obj.dudx_term+obj.dudz_term+...
            obj.dwdx_term+2.*obj.dwdz_term+2.*obj.uzwx_term);
            obj.epsilon_dvdy_dwdz_avg = mean(obj.epsilon_dvdy_dwdz,'all','omitnan'); % dissipation rate spatial average
            obj.tau_kt_dvdy_dwdz = (obj.nu/obj.epsilon_avg)^0.5; % time (s)
            obj.eta_kl_dvdy_dwdz = (obj.nu^3/obj.epsilon_avg)^0.25; %length (m)
            
            %integrated dissipation spectrum
            R = (deltax/obj.eta_kl_dvdy_dwdz); 
            x_pos = 2*pi/R;
            
%             %correction value
%             cvalue = str2num(cell2mat(inputdlg('Value on chart for the correction in %?',...
%              'Dissipation Correction', [1 50])));

             f = figure('Units','Normalized',...
                 'Position',[.1 .1 .05 .05],...
                 'NumberTitle','off',...
                 'Name','Dissipation Correction');
            imshow('universalspectrum.png')
            title('Interation of Universal Spectrum, Enter % value.');
            e = uicontrol('Style','Edit',...
                 'Units','Normalized',...
                 'Position',[.2 .4 .1 .1],...
                 'Tag','myedit');
            p = uicontrol('Style','PushButton',...
                 'Units','Normalized',...
                 'Position',[.3 .4 .1 .1],...
                 'String','Enter',...
                 'CallBack','uiresume(gcbf)');
            annotation('textbox',[.2 .5 .4 .1],'String',['2*\pi/R: ',num2str(x_pos,3)],'FitBoxToText','on','BackgroundColor','w');

            uiwait(f)
            cvalue = str2double(get(e,'String'));
            uiresume(f)
            
            close all           
            
            %corrected values
            obj.epsilon_dvdy_dwdz_corrected = obj.epsilon_dvdy_dwdz*(2-cvalue/100);
            obj.epsilon_avg_dvdy_dwdz_corrected = mean(obj.epsilon_dvdy_dwdz_corrected,'all','omitnan'); 
            obj.epsilon_median_dvdy_dwdz_corrected = median(obj.epsilon_dvdy_dwdz_corrected,'all','omitnan'); 

            % Kolmogorov length and time scale
            obj.tau_kt_dvdy_dwdz_corrected = (obj.nu/obj.epsilon_avg_dvdy_dwdz_corrected)^0.5; % time
            obj.eta_kl_dvdy_dwdz_corrected = (obj.nu^3/obj.epsilon_avg_dvdy_dwdz_corrected)^0.25; %length (cm)      
            
            %%%%% dissipation plot of all 3 assumptions
            figure (9) 
            subplot(1,3,1)
            imagesc((obj.epsilon_corrected), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            title(['$Continuity\:-\: \overline{\epsilon}\::\:',(num2str(obj.epsilon_avg_corrected,3)),'\:cm^2/s^3\;\;M_d:\:',(num2str(obj.epsilon_median_corrected,3)),'\:cm^2/s^3$'],'Interpreter','latex', 'FontSize',13)
            ylabel('cm')
            xlabel('cm')
            subplot(1,3,2)
            imagesc((obj.epsilon_dvdy_dudx_corrected), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')      
            title(['$\frac{\partial v}{\partial y}=\frac{\partial u}{\partial x}\:-\: \overline{\epsilon}\::\:',(num2str(obj.epsilon_avg_dvdy_dudx_corrected,3)),'\:cm^2/s^3\;\;M_d:\:',(num2str(obj.epsilon_median_dvdy_dudx_corrected,3)),'\:cm^2/s^3$'],'Interpreter','latex', 'FontSize',13)
            ylabel('cm')
            xlabel('cm')
            subplot(1,3,3)
            imagesc((obj.epsilon_dvdy_dwdz_corrected), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')    
            title(['$\frac{\partial v}{\partial y}=\frac{\partial w}{\partial x}\:-\: \overline{\epsilon}\::\:',(num2str(obj.epsilon_avg_dvdy_dwdz_corrected,3)),'\:cm^2/s^3\;\;M_d:\:',(num2str(obj.epsilon_median_dvdy_dwdz_corrected,3)),'\:cm^2/s^3$'],'Interpreter','latex', 'FontSize',13)
            ylabel('cm')
            xlabel('cm')
            set(gcf,'Position',[700 300 1240 300])


        end
        end
end