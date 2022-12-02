classdef PIVanalysis < handle
    properties
        % list all variables that will be passed between methods in the
        % model here. think of them as global variables within FFD. if a
        % value is included in the definition then it is treated as the
        % default value.
        
        % static model parameters with default values
  
        %g = 9.80665                     % (m/s2)                           
        mu = 2.5*1.81e-5                % (kg/m s) dynamic viscosity, u
        nu = 0.8927*10-6;               % (m2/s) % https://www.omnicalculator.com/physics/water-viscosity

        calibration
        u_original
        v_original
        datamax
        datamin
        original_percentleft
        agwpercentleft
        medianfilter_percentleft
        interpolated_percentleft
        m1             % mean flow strength in x direction
        m1_avg         % spatial average
        m3             % mean flow strength in z direction
        m3_avg         % spatial average
        mstar          % relative mean flow strength
        mstar_avg      % spatial average of relative mean flow strength
        u_agw
        v_agw
        u_nanfilter
        v_nanfilter
        u_interpolated
        v_interpolated
        u_mean
        u_mean_savg
        v_mean
        v_mean_savg
        ugrad               % u velocity gradient
        u_f
        v_f
        u_rms
        v_rms
        u_rms_savg
        v_rms_savg
        tke                 % turbulent kinetic energy
        tke_avg
        isotropy
        isotropy_avg        % spatial average
        integral_avg        % integral length scale spatial average
        epsilon             % dissipation
        epsilon_avg         % dissipation spatial average
        target
        tau_kt              % Kolmogorov time scale (sec) 
        eta_kl              % Kolmogorov length scale (time)
        lambda_tm           % Taylor microscale (centimeters) 
        Re_lambda           % Taylor scale Reynolds number

        yaxis
        xaxis
    
                          
        % nondimensional terms
        Re                  % Reynolds number w.r.t. Uinf
        Fr                  % Freud number w.r.t Uinf
        % figures and tables
        vtbl            % table containing all static variables
        ubarfig         % figure showing velocity in top boundary
        wbarfig         % figure showing velocity in center boundary
        computeBUfig    % figure showing approximated bulk velocities
    end      
    methods 
        % methods are defined to operate on class properties. the entire
        % framework for the numerical simulation will be constructed here.
        function obj = PIVanalysis()
            % include any initializations for properties here. It will be
            % ran whenever a new class instantiation is performed.
            
            
            %load('G:\J 20 All jets 2021\3_5 volts 1 sec 40 3_5 ms (redo)\PIVlab_results','u_original','v_original','calxy')
            %load('G:\J 20 All jets 2021\3_5 volts 1 sec 50 3_5 ms (new)\PIVlab_results','u_original','v_original')
            %load('G:\J 20 All jets 2021\3_5 volts 1 sec 60 3_5 ms (new)\PIVlab_results','u_original','v_original')
            %load('G:\J 20 All jets 2021\4 volts 0_6 sec 40 3_5 ms (new)\PIVlab_results','u_original','v_original')
            %load('G:\J 20 All jets 2021\4 volts 0_6 sec 50 3 ms (new)\PIVlab_results','u_original','v_original')
            %load('G:\J 20 All jets 2021\4 volts 0_6 sec 60 3_5 ms (new)\PIVlab_results','u_original','v_original')
            %load('G:\J 20 All jets 2021\4 volts 1 sec 40 3_5 ms (redo)\PIVlab_results','u_original','v_original')
            %load('G:\J 20 All jets 2021\4 volts 1 sec 50 3 ms (new)\PIVlab_results','u_original','v_original')
            %load('G:\J 20 All jets 2021\4 volts 1 sec 60 3 ms (new)\PIVlab_results','u_original','v_original')
            %load('G:\J 20 All jets 2021\4_5 volts 0_6 sec 40 3 ms (new)\PIVlab_results','u_original','v_original')
            %load('G:\J 20 All jets 2021\4_5 volts 0_6 sec 50 3 ms (new)\PIVlab_results','u_original','v_original')
            %load('G:\J 20 All jets 2021\4_5 volts 0_6 sec 60 3 ms (new)\PIVlab_results','u_original','v_original')
            %load('G:\J 20 All jets 2021\4_5 volts 1 sec 40 3 ms (redo)\PIVlab_results','u_original','v_original')
            %load('G:\J 20 All jets 2021\4_5 volts 1 sec 50 3 ms (new)\PIVlab_results','u_original','v_original')
            %load('G:\J 20 All jets 2021\4_5 volts 1 sec 60 3 ms (new)\PIVlab_results','u_original','v_original')
            
            %%new tests Nov 2022
            %load('G:\J20 Tests\4_5V_1sec_60percent_2_5ms_withlid\PIVlab_data','u_original','v_original')
            %load('G:\J20 Tests\4_5V_1s_60percent_2_5ms_no_lid\PIVlab_results','u_original','v_original')
            %load('H:\Aubrey Data\J20_Test\3_5V_1s_40percent_3_5ms_no_lid\PIVlab_results','u_original','v_original')
            %load('G:\J20 Tests\3_5V_40percent_1sec_4ms_cornerjets\PIVlab_results','u_original','v_original')
            %load('H:\Aubrey Data\J20_Test\3_5V_1s_40percent_4_5ms_no_lid_cornerjetsfixed\PIVlab_results','u_original','v_original')
            %load('H:\Aubrey Data\J20_Test\3_5V_1s_40percent_4_5ms_no_lid_cornerjetsfixed_nojetmesh\PIVlab_results','u_original','v_original','calxy')
            %load('G:\J20 Tests\4_5V_1s_40percent_4_5ms_nomeshcornerjets\PIVlab_results','u_original','v_original','calxy')
            %load('H:\Aubrey Data\J20_Test\4_5V_1s_40percent_5ms_nomeshcornerjets_somejetschanged\PIVlab_results','u_original','v_original','calxy')
            load('H:\Aubrey Data\J20_Test\4_5V_1s_40percent_4_5ms_nomeshcornerjets_somejetschanged2_10mintues\PIVlab_results','u_original','v_original','calxy')
           
            obj.u_original = u_original;
            obj.v_original = v_original;
            obj.calibration = calxy; % Ex: 4e-05 m/px 
            
            
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        function allfunctions(obj)
            reshapes(obj); % %this works well for PIVlab sessions
            checkHistogram(obj);
            applyAGWfilter(obj);
            medfilter(obj);
            velocityCalculations(obj);
            %delaunyinterpolation(obj); 
            %spatialspectra(obj);
            %temporalspectra(obj);
            
        end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        function reshapes(obj)
         obj.u_original = cell2mat(permute(obj.u_original',[1,3,2])); %No filter is done within PIVLAB
         obj.v_original = cell2mat(permute(obj.v_original',[1,3,2])); %Has NaNs
         
%          obj.u_original = obj.u_original(:,32:217,:); 
%          obj.v_original = obj.v_original(:,32:217,:);

%          obj.u_original = obj.u_original(17:169,49:200,:); 
%          obj.v_original = obj.v_original(17:169,49:200,:);
% 
         obj.u_original = obj.u_original(:,:,1:2); 
         obj.v_original = obj.v_original(:,:,1:2);
         
         %Matches NaN values for both u and w (i.e. if u has a NaN value at
         %(1,1,1) and v does not, these lines will assign a NaN value at
         %(1,1,1) for v to match u. 
         obj.u_original(isnan(obj.v_original)) = NaN;
         obj.v_original(isnan(obj.u_original)) = NaN;
         
        end
        function checkHistogram(obj)
            [Ny,Nx,Nt] = size(obj.u_original);
            total = Ny*Nx*Nt; 
            logicarray = ~isnan(obj.u_original);
            sumofnonNaN = sum(sum(sum(logicarray)));
            
            obj.original_percentleft = (sumofnonNaN/total)*100;

            ucheck = reshape(obj.u_original,1,Nx*Ny*Nt); vcheck = reshape(obj.v_original,1,Nx*Ny*Nt);

            figure(1)
            title('Histograms of the $u$ and w velocities - Pre-Filter')
            subplot(2,1,1)
            histogram(ucheck,100)
            title('Histogram of the $u$ velocities - Pre-Filter','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('m/s')
            subplot(2,1,2)
            histogram(vcheck,100)
            title('Histogram of the $w$ velocities - Pre-Filter','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('m/s')
            
            obj.datamax = str2num(cell2mat(inputdlg('Enter a number:',...
             'AGW Filter Max. and Min. Value', [1 50])));
            obj.datamin=-obj.datamax;
            close all
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function applyAGWfilter(obj)
            
            [Ny,Nx,Nt] = size(obj.u_original); 

            ucheck = reshape(obj.u_original,1,Nx*Ny*Nt); vcheck = reshape(obj.v_original,1,Nx*Ny*Nt);

            obj.u_original=permute(obj.u_original,[3 2 1]); obj.v_original=permute(obj.v_original,[3 2 1]); 

            [Nt,Ny,Nx] = size(obj.u_original);

            uresh = reshape(obj.u_original,Nt*Nx,Ny); vresh = reshape(obj.v_original,Nt*Nx,Ny);
            %
                time=1:Nt*Nx;

                ufill = zeros(Ny,Nt*Nx); %creates double array with the same number of positions as the original data
                vfill = zeros(Ny,Nt*Nx);

                ufillvect=5; %not sure what this does or why it is 5

                for i = 1:Ny

                    utemp = uresh(:,i);
                    vtemp = vresh(:,i);
                    [udat utim] = agw_filter(utemp,time,obj.datamax,obj.datamin);
                    [vdat vtim] = agw_filter(vtemp,time,obj.datamax,obj.datamin);


                    %NaNs become 1000 and removed points become 1000.
                    %This uses utim/vtim to put the udat/vdat in the correct index of
                    %the original data
                    clear ufillvect
                    ufillvect = zeros(1,Nx*Nt);
                    ufillvect(utim)=udat;
                    ufill(i,:)=ufillvect; %the data that was not filtered, udat/vdat, is put into ufill/vfill
                    %keep their same index, where all other positions, which has data
                    %that filtered out, has a 1000 value

                    clear vfillvect
                    vfillvect = zeros(1,Nx*Nt);
                    vfillvect(vtim)=vdat;
                    vfill(i,:)=vfillvect;

                end
            %
            ufill = ufill'; obj.u_agw = reshape(ufill,Nt,Ny,Nx);
            vfill = vfill'; obj.v_agw = reshape(vfill,Nt,Ny,Nx);
            
            obj.u_agw(obj.u_agw==0)=NaN;
            obj.v_agw(obj.v_agw==0)=NaN;
            

            obj.u_agw(isnan(obj.v_agw)) = NaN;
            obj.v_agw(isnan(obj.u_agw)) = NaN;
            
            [Nt,Ny,Nx] = size(obj.v_agw);
            uaftercheck = reshape(obj.v_agw,1,Nx*Ny*Nt); vaftercheck = reshape(obj.v_agw,1,Nx*Ny*Nt);

            % Percentage of velocities left after the AGW filter has been applied, I think Blair said this should be 90 percent or higher
            % if the data is good
            [Ny,Nx,Nt] = size(obj.u_original);
            total = Ny*Nx*Nt; 
            logicarray = ~isnan(obj.u_agw);
            sumofnonNaN = sum(sum(sum(logicarray)));
            
            obj.agwpercentleft = (sumofnonNaN/total)*100; %it is the same for all components
            
            %f = msgbox(num2str(obj.agwpercentleft*100), 'Percent Remaining. Ideally 95% or more.')

            %----------------------
            % Histograms to compare pre and post filtering
            %----------------------
            figure(1)
            %title('Histograms of the $u$ and $w$ velocities---Pre and Post-Filter')
            subplot(2,2,1)
            
            histogram(ucheck,100)
            title('Histogram of the $u$ velocities','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('$u$ (m/s)','Interpreter','Latex')
            annotation('textbox',[0 .8 .1 .2], ...
    'String','Pre-Filter    ','EdgeColor','none','Interpreter','Latex','fontsize', 18)
            annotation('textbox',[0 .3 .1 .2], ...
    'String','Post-Filter    ','EdgeColor','none','Interpreter','Latex','fontsize', 18)
            subplot(2,2,3)
            histogram(uaftercheck,100)
            title('Histogram of the $u$ velocities','Interpreter','Latex')
            ylabel('Frequency') 
            xlabel('$u$ (m/s)','Interpreter','Latex')
            subplot(2,2,2)
            histogram(vcheck,100)
            title('Histogram of the $w$ velocities','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('$w$ (m/s)','Interpreter','Latex')
            subplot(2,2,4)
            h = histogram(vaftercheck,100);
            %set(h,'XData', obj.datamin:0.1:obj.datamax, 'YData',0:5000:max(h.Values))
%             xlim([obj.datamin-0.2 obj.datamax])
%             xticks(obj.datamin:0.1:obj.datamax)
            title('Histogram of the $w$ velocities','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('$w$ (m/s)','Interpreter','Latex')
            set(gcf,'Position',[700 300 800 700])
%             ax = gca;
%             ax.YAxis.Exponent = 3;
            hold off
            
            pause;
            
        end
        function onlymedfilter(obj) % for K PIVlab data
            obj.target = 4; %initial guess 
            
            obj.u_original = obj.u_original(35:186, 1:161, :); %to do the needed reshape
            obj.v_original = obj.v_original(35:186, 1:161, :);
            
           %[Nt,Ny,Nx] = size(obj.u_original);
           [Ny,Nx,Nt] = size(obj.u_original);
           
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
                
%                 ut=squeeze(obj.u_original(tt,:,:));
%                 vt=squeeze(obj.v_original(tt,:,:));
                
                ut=squeeze(flip(obj.u_original(:,:,tt)));
                vt=squeeze(flip(obj.v_original(:,:,tt)));

                umed = mediannan(ut,3);
                vmed = mediannan(vt,3);

                flagu = abs(umed - ut) > obj.target;
                flagv = abs(vmed - vt) > obj.target;
                flag = flagu + flagv;
                flag(flag==2)=1;

                utclean = ut.*(1-flag);
                utclean(utclean==0) = NaN;
                vtclean = vt.*(1-flag);
                vtclean(vtclean==0) = NaN;
                
%                 uclean_save(tt,:,:) = utclean;
%                 vclean_save(tt,:,:) = vtclean;
                uclean_save(:,:,tt) = utclean;
                vclean_save(:,:,tt) = vtclean;
                

                figure(2); %shows original data
%                 uoriginal = squeeze(obj.u_original(tt,:,:));
%                 voriginal = squeeze(obj.v_original(tt,:,:));
                uoriginal = squeeze(flip(obj.u_original(:,:,tt)));
                voriginal = squeeze(flip(obj.v_original(:,:,tt)));
                
                scale_factor = 0.1;
                h1=quiver(uoriginal*scale_factor.*(flip(obj.A)),voriginal*scale_factor.*(flip(obj.A)),'r','AutoScale','off');
                xlim([0 Nx])
                ylim([0 Ny])
                set(gcf,'Position',[700 300 800 700])
                hold on;
                
%                 %shows what agw filter removes
%                 uagw = squeeze(obj.u_agw(tt,:,:));
%                 vagw = squeeze(obj.v_agw(tt,:,:));
%                 h2= quiver(uagw*scale_factor,vagw*scale_factor,'r','AutoScale','off');
%                 xlim([0 Nx])
%                 ylim([0 Ny])
%                 hold on;
                
                %the values after medfilter
                h3 = quiver(utclean*scale_factor.*(flip(obj.A)),vtclean*scale_factor.*(flip(obj.A)),'k','AutoScale','off'); 
                xlim([0 Nx])
                ylim([0 Ny])
                title(['Time step: ', num2str(tt), ' out of ', num2str(Nt), '. Target: ', num2str(obj.target)])
                lgd = legend('Original','AfterMedian','Location','northeastoutside');
                
                hold off
                scale=3;
                hU1 = get(h1,'UData');
                hV1 = get(h1,'VData');
                set(h1,'UData',scale*hU1,'VData',scale*hV1)
%                 hU2 = get(h2,'UData');
%                 hV2 = get(h2,'VData');
%                 set(h2,'UData',scale*hU2,'VData',scale*hV2)
                hU3 = get(h3,'UData');
                hV3 = get(h3,'VData');
                set(h3,'UData',scale*hU3,'VData',scale*hV3)                
                
                if m==3
                    ttstart=1;
                    continue    
                else    
                    m = menu('Yes if to continue through time, No for new target value (originally 0.15), All to apply filter at all time steps (reset to T.S. =1 first), Exit to stop program., Skip to TimeStep','Yes','No', 'All', 'Exit', 'Select T.S.');
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
            obj.v_nanfilter = vclean_save;
            
            obj.u_nanfilter = obj.u_nanfilter.*(flip(obj.A));
            obj.v_nanfilter = obj.v_nanfilter.*(flip(obj.A));
            
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
            obj.target = 0.15; %initial guess 
            
            [Nt,Ny,Nx] = size(obj.u_agw);
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
                ut=squeeze(obj.u_agw(tt,:,:));
                vt=squeeze(obj.v_agw(tt,:,:));

                umed = mediannan(ut,3);
                vmed = mediannan(vt,3);

                flagu = abs(umed - ut) > obj.target;
                flagv = abs(vmed - vt) > obj.target;
                flag = flagu + flagv;
                flag(flag==2)=1;

                utclean = ut.*(1-flag);
                utclean(utclean==0) = NaN;
                vtclean = vt.*(1-flag);
                vtclean(vtclean==0) = NaN;
                
%                 uclean_save = nan(Nt,Ny,Nx); %preallocate 
%                 vclean_save = nan(Nt,Ny,Nx);
                uclean_save(tt,:,:) = utclean;
                vclean_save(tt,:,:) = vtclean;
                

                figure(2); %shows original data
                uoriginal = squeeze(obj.u_original(tt,:,:));
                voriginal = squeeze(obj.v_original(tt,:,:));
                scale_factor = 1;
                h1=quiver(uoriginal*scale_factor,voriginal*scale_factor,'b','AutoScale','off');
                xlim([0 Nx])
                ylim([0 Ny])
                set(gcf,'Position',[700 300 800 700])
                hold on;
                
                %shows what agw filter removes
                uagw = squeeze(obj.u_agw(tt,:,:));
                vagw = squeeze(obj.v_agw(tt,:,:));
                h2= quiver(uagw*scale_factor,vagw*scale_factor,'r','AutoScale','off');
                xlim([0 Nx])
                ylim([0 Ny])
                hold on;
                
                %the values after agw and medfilter
                h3 = quiver(utclean*scale_factor,vtclean*scale_factor,'k','AutoScale','off'); 
                xlim([0 Nx])
                ylim([0 Ny])
                title(['Time step: ', num2str(tt), ' out of ', num2str(Nt), '. Target: ', num2str(obj.target)])
                lgd = legend('AGW','Median','Remaining','Location','northeastoutside');
                
                hold off
                scale=3;
                hU1 = get(h1,'UData');
                hV1 = get(h1,'VData');
                set(h1,'UData',scale*hU1,'VData',scale*hV1)
                hU2 = get(h2,'UData');
                hV2 = get(h2,'VData');
                set(h2,'UData',scale*hU2,'VData',scale*hV2)
                hU3 = get(h3,'UData');
                hV3 = get(h3,'VData');
                set(h3,'UData',scale*hU3,'VData',scale*hV3)                
                
                if m==3
%                     if tt == selecttt
%                         m = menu('Yes if to continue through time, No for new target value (originally 0.15), All to apply filter at all time steps, Exit to stop program., Skip to TimeStep','Yes','No', 'All', 'Exit', 'Select T.S.');
%                     else
%                         continue
%                     end
                    continue    
                else    
                    m = menu('Yes if to continue through time, No for new target value (originally 0.15), All to apply filter at all time steps, Exit to stop program., Skip to TimeStep','Yes','No', 'All', 'Exit', 'Select T.S.');
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
                selecttt = str2double(cell2mat(inputdlg('Enter new time step:',...
             'Skip to this time step', [1 50])));
            else
                obj.target = str2double(cell2mat(inputdlg('Enter new target:',...
            'Target', [1 50])));
            end 

            
            %check if the target should be updated then, depending on how
            %the plots looking 
            
            
            end
            obj.u_nanfilter = uclean_save;
            obj.v_nanfilter = vclean_save;
            
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
            
            bs_wfull=sort(bootstrp(1000,'mean',obj.v_nanfilter));
            bs_wlower =  bs_wfull(25); bs_wupper =  bs_wfull(975);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function delaunyinterpolation(obj) 
            
            unan = obj.u_nanfilter;
            vnan = obj.v_nanfilter;
            
            [Nt, Ny, Nx] = size(unan);
            
            %replace any corner NaNs with the value of the median of existing corners
            %so that interpolation if inner parts can take place

            %upper left
            U_UL=nanmedian(unan(:,1,1));
            V_UL=nanmedian(vnan(:,1,1));
            %upper right
            U_UR=nanmedian(unan(:,1,Nx));
            V_UR=nanmedian(vnan(:,1,Nx));
            %lower left
            U_LL=nanmedian(unan(:,Ny,1));
            V_LL=nanmedian(vnan(:,Ny,1));
            %lower right
            U_LR=nanmedian(unan(:,Ny,Nx));
            V_LR=nanmedian(vnan(:,Ny,Nx));
            
            
            for tt=1:Nt
                tt
                unan_temp=squeeze(unan(tt,:,:));
                vnan_temp=squeeze(vnan(tt,:,:));
                [uuu,vvv] = fillNaNsTemp(unan_temp,vnan_temp, U_UL, U_UR, U_LL, U_LR, V_UL, V_UR, V_LL, V_LR);
                u_interp(tt,:,:)=uuu;
                v_interp(tt,:,:)=vvv;
                
%                 scale_factor = 1;
%                 figure(3)
% 
%                 h1 = quiver(uuu*scale_factor,vvv*scale_factor,'r','AutoScale','off'); 
%                 xlim([0 Nx])
%                 ylim([0 Ny]) 
%                 hold on;
%                 
%                 unanf = squeeze(obj.u_nanfilter(tt,:,:));
%                 vnanf = squeeze(obj.v_nanfilter(tt,:,:));
%                 h2 = quiver(unanf*scale_factor,vnanf*scale_factor,'k','AutoScale','off'); 
%                 xlim([0 Nx])
%                 ylim([0 Ny])
%                 title(['Time step: ', num2str(tt), ' out of ', num2str(Nt)])
%                 lgd = legend('Interpolated','Filtered','Location','northeastoutside');
%                 
%                 hold off
%                 scale=3;
%                 hU1 = get(h1,'UData');
%                 hV1 = get(h1,'VData');
%                 set(h1,'UData',scale*hU1,'VData',scale*hV1)
%                 hU2 = get(h2,'UData');
%                 hV2 = get(h2,'VData');
%                 set(h2,'UData',scale*hU2,'VData',scale*hV2)
%                 
%                 pause
                close all
            end      
            
            obj.u_interpolated = u_interp(:,:,:);
            obj.v_interpolated = v_interp(:,:,:);    
            
            [Ny,Nx,Nt] = size(obj.u_original);
            total = Ny*Nx*Nt; 
            logicarray = ~isnan(obj.u_interpolated);
            sumofnonNaN = sum(sum(sum(logicarray)));
            obj.interpolated_percentleft = (sumofnonNaN/total)*100;
            
            beep
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function velocityCalculations(obj) % and tke, isotropy, and mean flow strength
            
            % Bring in data, the is for a 3D double array, mine is 53 (height) by 79 (length) by 10,500 (in time)
            u_o = obj.u_nanfilter; v_o = obj.v_nanfilter;
            
            u_o=permute(u_o,[3 2 1]); v_o=permute(v_o,[3 2 1]); 

            
            [Ny,Nx,Nt] = size(u_o);
            
            obj.u_mean = nanmean(u_o,3); obj.v_mean = nanmean(v_o,3);
            obj.u_mean_savg = nanmean(nanmean(obj.u_mean)); obj.v_mean_savg = nanmean(nanmean(obj.v_mean));
            obj. u_f = u_o-obj.u_mean; obj.v_f = v_o-obj.v_mean;
            obj.u_rms = sqrt(nanmean((obj.u_f.^2),3)); obj.v_rms = sqrt(nanmean((obj.v_f.^2),3));
            obj.u_rms_savg = nanmean(nanmean(obj.u_rms)); obj.v_rms_savg = nanmean(nanmean(obj.v_rms));
            tke = 0.5*(2*(obj. u_f.^2) + (obj. v_f.^2));
            obj.tke = nanmean(tke,3);
            obj.tke_avg = nanmean(obj.tke,'all'); 
            obj.isotropy = obj.u_rms./obj.v_rms; 
            obj.isotropy_avg = nanmean(obj.isotropy,'all'); 

            obj.yaxis = [-Ny/2:20:Ny/2]*obj.calibration*100*16; %converts subwindow count to cm
            obj.xaxis = [-Nx/2:20:Nx/2]*obj.calibration*100*16; %converts subwindow count to cm
            
            %Time average for each subwindow                            
            figure (1) 
            subplot(2,2,1)
            imagesc(flipud(obj.u_mean*100), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            %caxis([-0.02,0.02])
            title('Colorbar U-Velocity (cm/s)','Interpreter','Latex')
            ylabel('cm')
            xlabel('cm')

            subplot(2,2,2)
            imagesc(flipud(obj.v_mean*100), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            %caxis([-0.02,0.02])
            title('Colorbar W-Velocity (cm/s)','Interpreter','Latex')
            ylabel('cm')
            xlabel('cm')

            %RMS velocity for each subwindow
            subplot(2,2,3)
            imagesc(flipud(obj.u_rms*100), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            %caxis([-0.02,0.02])
            title(['u_{rms} (cm/s), Spatial Avg: ',(num2str(obj.u_rms_savg,3))])
            ylabel('cm')
            xlabel('cm')

            subplot(2,2,4)
            imagesc(flipud(obj.v_rms*100), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            %caxis([-0.02,0.02])
            title(['w_{rms} (cm/s), Spatial Avg: ',(num2str(obj.v_rms_savg,3))])
            ylabel('cm')
            xlabel('cm')
            set(gcf,'Position',[700 300 1000 700])

            figure (2)
            imagesc(flipud(obj.tke*10000), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            %caxis([-0.02,0.02])
            title(['k (cm^2/s^2), Spatial Avg: ',(num2str(obj.tke_avg,3))])
            ylabel('cm')
            xlabel('cm')
            
            figure (3)
            imagesc(flipud(obj.isotropy), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            %caxis([-0.02,0.02])
            title(['Isotropy (u_{rms}/w_{rms}), Spatial Avg: ',(num2str(obj.isotropy_avg,3))])
            ylabel('cm')
            xlabel('cm')
            
            %mean flow strength 
            obj.m1 = obj.u_mean./obj.u_rms; 
            obj.m3 = obj.v_mean./obj.v_rms;
            obj.m1_avg = nanmean(obj.m1,'all');
            obj.m3_avg = nanmean(obj.m3, 'all');
            
            obj.mstar = (0.5*(2*(obj.u_mean.^2) + (obj.v_mean.^2)))./(obj.tke); %time average U, W, and tke
            obj.mstar_avg = nanmean(obj.mstar,'all');
            
            figure (4) 
            subplot(1,3,1)
            imagesc(flipud(obj.m1), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            %caxis([-0.02,0.02])
            title(['M_1, Spatial Avg: ',num2str(obj.m1_avg*100,3),' %'])
            ylabel('cm')
            xlabel('cm')

            subplot(1,3,2)
            imagesc(flipud(obj.m3), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            %caxis([-0.02,0.02])
            title(['M_3, Spatial Avg: ',num2str(obj.m3_avg*100,3),' %'])
            ylabel('cm')
            xlabel('cm')
            
            subplot(1,3,3)
            imagesc(flipud(obj.mstar), 'XData', obj.xaxis, 'YData',obj.yaxis)
            colorbar
            set(gca,'YDir','normal')
            %caxis([-0.02,0.02])
            title(['M^*, Spatial Avg: ',num2str(obj.mstar_avg*100,3),' %'])
            ylabel('cm')
            xlabel('cm')

            set(gcf,'Position',[700 300 1240 300])
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function kolmogorovscales(obj) %dissipation needs to fixed
            obj.tau_kt = (obj.nu/obj.epsilon_avg)^0.5; % time
            
            obj.eta_kl = (obj.nu^3/obj.epsilon_avg)^0.25; %length (cm)
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function taylorscales(obj) %fix dissipation and intregral length scale
            %possibly fix tke too

            % Taylor microscale (centimeters) 
            obj.lambda_tm  = sqrt(10)*obj.eta_kl^(2/3)* obj.integral_avg^(1/3);
            
            % Taylor scale Reynolds number
            obj.Re_lambda = (2/3)*obj.tke_avg*sqrt(15/(obj.nu*obj.epsilon_avg));
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function spatialspectra(obj)
            
            u_o = obj.u_interpolated; v_o = obj.v_interpolated;
            
            u_o=permute(u_o,[3 2 1]); v_o=permute(v_o,[3 2 1]);             
            
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
                vspatial = v_o(b,:,i);  
                Suu = fft(uspatial).*conj(fft(uspatial));
                Svv = fft(vspatial).*conj(fft(vspatial));
            end 
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function temporalspectra(obj) %this mostly works but there are issues with the number of time steps
            %right it is set to ensemble average of 10 runs
            u_o = obj.u_interpolated; v_o = obj.v_interpolated;
            
            u_o=permute(u_o,[3 2 1]); v_o=permute(v_o,[3 2 1]); 
            
            %u_o(isnan(u_o)) = 0; v_o(isnan(v_o)) = 0;

            %436(left to right), 288(top going dowwards)
            %so approx subwindow (18,27)
            xlocation = 18;
            ylocation = 27;

            u = u_o;
            v = v_o;
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
            u_o = obj.u_nanfilter; v_o = obj.v_nanfilter;
            
            u_o=permute(u_o,[1 3 2]); v_o=permute(v_o,[1 3 2]); 
            %u_o=permute(u_o,[2 3 1]); v_o=permute(v_o,[2 1 3]);
            
            subwindow=32;    %pixels between subwindow centers 
            res=0.0105;
            lowb=1;
            uppb=64; %not sure?
            
            [time height width]=size(u_o);
            height_subw=height:-1:1;
            height_subw=height_subw.*res.*subwindow; %calibrate to cm
            
            center=(width+1)/2;
            radd=1:(center-2);
            radd=radd.*subwindow.*res;
            radd2(2:uppb)=radd.*2;
            radd2(1)=0;

            %autocorrelation calculation          
            for row=1:height-7 %height-9 1:height-2, can't do height b/c otherwise it crashes...
                for radius=1:center-2
                    %calculate a(r) from VarianoCowen2008 formula 4.1
                    au(row,radius)=nanmean(udata(:,row,center-radius).*udata(:,row,center+radius))...
                        ./sqrt(nanmean(udata(:,row,center-radius).^2).*nanmean(udata(:,row,center+radius).^2));

                    av(row,radius)=nanmean(vdata(:,row,center-radius).*vdata(:,row,center+radius))...
                        ./sqrt(nanmean(vdata(:,row,center-radius).^2).*nanmean(vdata(:,row,center+radius).^2));

                end
            end
                auplot(2:uppb)=au(row,:);
                auplot(1)=1;
                figure(120)
                hold off;
                plot(radd2,auplot,'b*');
        %         title(row);
                hold on;
                grid on;
        %         xlim([0 10])
        %         ylim([0 1]);
            
                newExpFuncG=@(rr,ll) exp(-rr./ll).*(1-(rr./(2.*ll)));          % Defines anon. func.
                L0=6;                                         % Starting guess for L
                g=@(ll)sum((avplot(lowb:uppb)-newExpFuncG(radd2(lowb:uppb),ll)).^2,2);            % Anon. function to be minimized
                LstarG(kp,row)=fminunc(g,L0);                            % Minimize to find Lstar!
                YhatG=newExpFuncG(radd2(lowb:uppb),LstarG(kp,row));                       % Predict values of Y
                Resid=YhatG-avplot(lowb:uppb);                                   % Calc. residuals
                SSR=sum(Resid.^2,2);                            % Sum of Squared Residuals
                SSE=sum((avplot(lowb:uppb)-mean(avplot(lowb:uppb))).^2);                        % Sum of Squared Errors
                RsquaredG(kp,row)=1-(SSR/SSE);

            
            
        end
        function dissipation(obj)  
            %dissipation is calculated using the direct method
            
            %integrated dissipation spectrum
            % eta is in mm, PIV resolution in mm
            R = deltax/eta;
            x_pos = 2*pi/R;
            % using x position to go up and across
            figure(2)
            plot(f,Gain, '-.b') % dissipation spectra
            hold on
            plot(f,G2, '--r') % integrated dissipation spectrum
            xlabel('k_\eta')
            ylabel('D(k)/u^3_\eta, \epsilon(0,k)/\epsilon_m')
            title('Normalized dissipation spectrum and cumulative dissipation')
            legend('Normalized dissipation spectrum','Cumulative dissipation')
            
        end
    end
end