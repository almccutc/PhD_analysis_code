classdef PIVanalysis < handle
    properties
        % list all variables that will be passed between methods in the
        % model here. think of them as global variables within FFD. if a
        % value is included in the definition then it is treated as the
        % default value.
        
        % static model parameters with default values
  
        g = 9.80665                     % (m/s2)                           
        mu = 2.5*1.81e-5                % (kg/m s) dynamic viscosity, u
        nu                              % (m2/s) 
 
        u_o
        v_o
        u_original
        v_original
        datamax
        datamin
        percentleft
        u_agw
        v_agw
        u_nanfilter
        v_nanfilter
        
        % matrices for velocity values
        mean                % intermediate r-momentum state matrix (KP-04/25)
        fluctuating         % intermediate z-momentum state matrix (KP-04/25)
        rms                 % intermediate r-momentum advection operator 
        tke                    % state matrix (KP-04-25)                            
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
            load('F:\032420_AM_0.6_40_105_630_10min_4.5V_Correct.mat','u_original','v_original') %Change% 
            %load('032420_AM_0.6_40_105_630_10min_4.5V_Correct.mat','u_original','v_original') %Change% 
            %load('032420_AM_01_80_120_720_10min_4.5V_Correct.mat','u_original','v_original')
            %load('4Vworkspace.mat','Utotal','Wtotal');
            %load('6Vworkspace.mat','Utotal','Wtotal');
            obj.u_original = u_original;
            obj.v_original = v_original;
            
        end
        function reInitObj(obj)
            % recomputes static variables. should be ran if any of the
            % system properties are changed externally to ensure that 
            % everything is consistent.
            
            %obj.u_o = obj.u_original;
            %obj.v_o = obj.v_original;
            
%             obj.ztop = obj.H;
%             obj.rbar = 1e-6:obj.drbar:obj.b;          
%             obj.zbar = 0:obj.dzbar:1;
%             
%             obj.rMaxIndex = size(obj.rbar,2);
%             obj.zMaxIndex = size(obj.zbar,2);
        end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        function reshapes(obj)
         obj.u_o = cell2mat(obj.u_original); %No filter is done within PIVLAB
         obj.v_o = cell2mat(obj.v_original); %These values have NaNs in them

        % % Rotate Matrix
         nlay = length(obj.u_original);
         [r,c] = size(obj.u_o);
        % 
         obj.u_o   = permute(reshape(obj.u_o',[c,r/nlay,nlay]),[2,1,3]);
         obj.v_o   = permute(reshape(obj.v_o',[c,r/nlay,nlay]),[2,1,3]);
         
        end
        function checkHistogram(obj)
            reshapes(obj);
            [Ny,Nx,Nt] = size(obj.u_o);

            ucheck = reshape(obj.u_o,1,Nx*Ny*Nt); vcheck = reshape(obj.v_o,1,Nx*Ny*Nt);

            figure(1)
            title('Histograms of the $u$ and w velocities---Pre-Filter')
            subplot(2,1,1)
            histogram(ucheck,100)
            title('Histogram of the $u$ velocities---Pre-Filter','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('m/s')
            subplot(2,1,2)
            histogram(vcheck,100)
            title('Histogram of the $w$ velocities---Pre-Filter','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('m/s')
            
            obj.datamax = str2num(cell2mat(inputdlg('Enter a number:',...
             'AGW Filter Max. and Min. Value', [1 50])));
            obj.datamin=-obj.datamax;
            close all
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function applyAGWfilter(obj)
            checkHistogram(obj)
            
            [Ny,Nx,Nt] = size(obj.u_o); 

            ucheck = reshape(obj.u_o,1,Nx*Ny*Nt); vcheck = reshape(obj.v_o,1,Nx*Ny*Nt);

            obj.u_o=permute(obj.u_o,[3 2 1]); obj.v_o=permute(obj.v_o,[3 2 1]); 

            [Nt,Ny,Nx] = size(obj.u_o);

             uresh = reshape(obj.u_o,Nt*Nx,Ny); vresh = reshape(obj.v_o,Nt*Nx,Ny);


                time=1:Nt*Nx;
                %0.6_40_105_630_10min_4.5V: 0.4 and -0.4, 0.8970
                %0.6_20_105_630_10min_4.5V: 0.4 and -0.4, 0.9341
                %0.6_60_105_630_10min_4.5V: 0.4 and -0.4, 0.8693
                %0.6_80_110_660_10min_4.5V: 0.4 and -0.4, 0.8773
                %01_20_115_690_10min_4.5V: 0.4 and -0.4, 0.8857
                %01_40_120_720_10min_4.5V: 0.4 and -0.4, 0.8980
                %01_60_120_720_10min_4.5V: 0.4 and -0.4, 0.8979
                %01_80_120_720_10min_4.5V: 0.4 and -0.4, 0.9418 and 0.9424

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
                    ufillvect = zeros(1,Nx*Nt)+1000;
                    ufillvect(utim)=udat;
                    ufill(i,:)=ufillvect; %the data that was not filtered, udat/vdat, is put into ufill/vfill
                    %keep their same index, where all other positions, which has data
                    %that filtered out, has a 1000 value

                    clear vfillvect
                    vfillvect = zeros(1,Nx*Nt)+1000;
                    vfillvect(vtim)=vdat;
                    vfill(i,:)=vfillvect;

                end

                ufill = ufill'; unew = reshape(ufill,Nt,Ny,Nx);

                vfill = vfill'; vnew = reshape(vfill,Nt,Ny,Nx);

                %if value is 1000, gets assigned a value of 1 at each index, otherwise is 0
                flaguAGW0=unew==1000; flagvAGW0=vnew==1000;
                flaggedAGW0 = flaguAGW0 + flagvAGW0;
                flaggedAGW0(flaggedAGW0==2)=1;

                unewnan = unew; vnewnan = vnew;

                %sets value of index to NaN if it has a flagged value of 1
                for tt=1:Nt
                    for row=1:Ny
                        for col=1:Nx
                            if flaggedAGW0(tt,row,col) > 0
                                unewnan(tt,row,col)=NaN;
                                vnewnan(tt,row,col)=NaN;
                            end
                        end
                    end
                end

            obj.u_agw=unewnan; obj.v_agw=vnewnan;

            % Percentage of velocities left after the AGW filter has been applied, I think Blair said this should be 90 percent or higher
            % if the data is good
            obj.percentleft = length(utim)/length(time); %it is the same for all components
            
            f = msgbox(num2str(obj.percentleft*100), 'Percent Remaining. Ideally 95% or more.')

            %----------------------
            % Histograms to compare pre and post filtering
            %----------------------
            figure(1)
            %title('Histograms of the $u$ and $w$ velocities---Pre and Post-Filter')
            subplot(2,2,1)
            hist(ucheck,100)
            title('Histogram of the $u$ velocities---Pre-Filter','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('$u$ (m/s)','Interpreter','Latex')
            subplot(2,2,2)
            hist(udat,100)
            title('Histogram of the $u$ velocities---Post-Filter','Interpreter','Latex')
            ylabel('Frequency') 
            xlabel('$u$ (m/s)','Interpreter','Latex')
            ax = gca;
            ax.YAxis.Exponent = 3;
            subplot(2,2,3)
            hist(vcheck,100)
            title('Histogram of the $w$ velocities---Pre-Filter','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('$w$ (m/s)','Interpreter','Latex')
            subplot(2,2,4)
            hist(vdat,100)
            title('Histogram of the $w$ velocities---Post-Filter','Interpreter','Latex')
            ylabel('Frequency')
            xlabel('$w$ (m/s)','Interpreter','Latex')
            
            medfilter(obj)
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function medfilter(obj)
            target = 10; %initial guess 
            
            [Nt,Ny,Nx] = size(obj.u_agw);
            
            while(1)
     
            m=1;
            
            for tt=1:Nt
                

                ut=squeeze(obj.u_agw(tt,:,:));
                vt=squeeze(obj.v_agw(tt,:,:));

                umed = mediannan(ut,3);
                vmed = mediannan(vt,3);

                flagu = abs(umed - ut) > target;
                flagv = abs(vmed - vt) > target;
                flag = flagu + flagv;
                flag(flag==2)=1;

                utclean = ut.*(1-flag);
                utclean(utclean==0) = NaN;
                vtclean = vt.*(1-flag);
                vtclean(vtclean==0) = NaN;

                uclean_save(tt,:,:) = utclean;
                vclean_save(tt,:,:) = vtclean;
         
                flagged(tt,:,:) = flag;
                
                figure(1); %shows original data
                hold off;
                unewtemp = squeeze(obj.u_original(tt,:,:));
                vnewtemp = squeeze(obj.v_original(tt,:,:));
                quiver(unewtemp,vnewtemp,2,'b');
                xlim([0 col])
                ylim([0 row])
                
                figure(1); %shows what agw filter removes
                hold off;
                unewtemp = squeeze(obj.u_agw(tt,:,:));
                vnewtemp = squeeze(obj.v_agw(tt,:,:));
                unewtemp(unewtemp==1000)=0;
                vnewtemp(vnewtemp==1000)=0;
                quiver(unewtemp,vnewtemp,2,'k');
                xlim([0 col])
                ylim([0 row])
                
                figure(1);
                hold on;
                quiver(utclean,vtclean,2,'r'); %the values after agw and medfilter
                xlim([0 col])
                ylim([0 row])
                
                if m==3
                    continue
                    else    
                    m = menu('Press yes no all','Yes','No', 'All');
                end 
                
                if m==2  % yes stored as 1, no stored as 2, no has a value of 3
                break;
                end
                
            end
            
            %check if the target should be updated then, depending on how
            %the plots looking 
            target = str2num(cell2mat(inputdlg('Enter new target or 10 (to exit medfilter function):',...
            'Target', [1 50])));
        
            if target==10
               break
            end
            
            end
            obj.u_nanfilter = uclean_save;
            obj.v_nanfilter = vclean_save;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function permute(obj)
            u_o = u_agw; v_o = v_agw;

            u_o=permute(u_o,[1 3 2]); v_o=permute(v_o,[1 3 2]); 

            u_o=permute(u_o,[2 3 1]); v_o=permute(v_o,[2 3 1]);

        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function velocityCalculations(obj)
            u_mean = nanmean(u_o,3); v_mean = nanmean(v_o,3);

            u_f = u_o-u_mean; v_f = v_o-v_mean;

            u_rms = sqrt(nanmean((u_f.^2),3)); v_rms = sqrt(nanmean((v_f.^2),3));

            obj.tke = 0.5*(2*(u_rms.^2) + (v_rms.^2));
            
            %Time average for each subwindow                            
            figure (1) 
            subplot(3,2,1)
            imagesc(u_mean)
            colorbar
            caxis([-0.02,0.02])
            title('Colorbar U-Velocity (m/s)','Interpreter','Latex')
            subplot(2,2,2)
            imagesc(v_mean)
            colorbar
            caxis([-0.02,0.02])
            title('Colorbar W-Velocity (m/s)','Interpreter','Latex')

            %Time rms average for each subwindow
            subplot(2,2,3)
            imagesc(u_rms)
            colorbar
            %caxis([-0.02,0.02])
            title('Colorbar rms (m/s)')
            subplot(2,2,4)
            imagesc(v_rms)
            colorbar
            %caxis([-0.02,0.02])
            title('Colorbar rms (m/s)')

            figure (2)
            imagesc(tke)
            colorbar
            %caxis([-0.02,0.02])
            title('TKE (m^2/s^2)')
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add any additional functions here
        % as a general rule, try to keep the number of operations for 
        % each method at a minimum. if a function is more than ~10 lines
        % then it should probably be split into multiple methods.
        % in general, the methods can be treated as a typical matlab
        % function with the acception of the method for passing in
        % properties. the general form for functions that have an output is
      
        % where the inputs are just arbitrary variables that aren't defined
        % as properties. it is typically a good practice to limit the usage
        % of these and instead just operate with class properties. if a
        % non-property variable is being passed into multiple functions 
        % then it should instead be defined as a property.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function animateSpeed(obj, plot_filename, rr, snapK)
%             % generates gif that displays velocity magnitude at every time
%             % step
%             if nargin < 4
%                 snapK = 1e10;
%             end
%             if nargin < 3
%                 rr = 0.5;
%             end
%             if nargin < 2
%                 plot_filename = 'speedPlot.gif';
%             end
%             figure('Units', 'normalized', ...
%             'Position', [0 0 0.4*obj.b 0.4], 'Visible', 'off');            
%             [~, tfigs] = plotZR(obj,  obj.uMag(obj.uk(1), obj.vk(1)), true);
%             ylabel(colorbar, '$\vert u/U_\infty\vert$', 'interpreter', 'latex', ...
%                 'FontSize', 14); 
%             caxis([0, sqrt(2)]);
%             gif(plot_filename, 'frame', gcf);
%             % update sequentially
%             for ii = 1:length(obj.tau)                                      
%                 gif;
%                 t = obj.tau2t(obj.tau(ii));
%                 u_ = obj.uMag(obj.uk(ii), obj.vk(ii));
%                 tfigs.ZData = u_; 
%                 title(sprintf('$t$ = %1.0f s', t), 'interpreter', 'latex', ...
%                 'FontSize', 14);
%                 pause(rr);
%                 % save plot if on save time-step
%                 if mod(ii, snapK) == 0
%                     captureSpeed(obj, ii);
%                 end
%             end            
%         end
        
%         function animateStress(obj, plot_filename, rr, snapK)
%             generates gif that displays velocity magnitude at every time
%             step
%             if nargin < 4
%                 snapK = 1e10;
%             end
%             if nargin < 3
%                 rr = 0.5;
%             end
%             if nargin < 2
%                 plot_filename = 'stressPlot.gif';
%             end
%             ur = obj.uk(1); uz = obj.vk(1);
%             [~, ~, tauR, tauZ] = stress(obj, reshape(ur', [], 1), ...
%                                                  reshape(uz', [], 1));
%             s = sqrt(tauR.^2 + tauZ.^2);                                          
%             figure('Units', 'normalized', ...
%             'Position', [0 0 0.4*obj.b 0.4], 'Visible', 'off');            
%             [~, tfigs] = plotZR(obj, s, true);
%             ylabel(colorbar, '$\vert \tau \vert$ (Pa)', 'interpreter', 'latex', ...
%                 'FontSize', 14); 
%             caxis([min(s(:)), max(s(:))]);
%             gif(plot_filename, 'frame', gcf);
%             update sequentially
%             for ii = 1:length(obj.tau)                                      
%                 gif;
%                 ur = obj.uk(ii); uz = obj.vk(ii);
%                 [~, ~, tauR, tauZ] = stress(obj, reshape(ur', [], 1), ...
%                                                  reshape(uz', [], 1));
%                 s = sqrt(tauR.^2 + tauZ.^2); 
%                 t = obj.tau2t(obj.tau(ii));
%                 tfigs.ZData = s; 
%                 title(sprintf('$t$ = %1.0f s', t), 'interpreter', 'latex', ...
%                 'FontSize', 14);
%                 pause(rr);
%                 save plot if on save time-step
%                 if mod(ii, snapK) == 0
%                     captureStress(obj, ii);
%                 end
%             end            
%         end
        
                  
%         function captureSpeed(obj, k)
%             % plots the non-dimensional speed profile over the entire
%             % domain at the indicated time step, k
%             % compute velocity magnitudes and time that will be ploted                                    
%             u_ = obj.uMag(obj.uk(k), obj.vk(k));
%             t = obj.tau2t(obj.tau(k));
%             [R, Z] = meshgrid(obj.rbar, obj.zbar); 
%             % plot contour
%             fzri = figure('Units', 'normalized', ...
%                 'Position', [0 0 0.4*obj.b 0.4], 'Visible', 'off');
%             pzri = surf(R, Z, u_);         
%             xlim([0, max(R(:))])
%             ylim([min(Z(:)), max(Z(:))])
%             pbaspect([obj.b, 1, 1]);  % figure sized proportional to aspect ratio
%             view(0, 90);
%             caxis([0, sqrt(2)]);
%             colormap(flipud(spring));
%             cb = colorbar;
%             cb.Ruler.MinorTick = 'on';
%             set(pzri, 'linestyle', 'none');
%             ylabel(cb, '$\vert u/U_\infty\vert$', 'interpreter', 'latex', 'FontSize', 14);
%             xlabel('$r/H$', 'interpreter', 'latex', 'FontSize', 14);
%             ylabel('$z/H$', 'interpreter', 'latex', 'FontSize', 14);
%             title(sprintf('$t$ = %1.0f s', t), 'interpreter', 'latex', ...
%                 'FontSize', 14);
%             if nargin > 2 && ShowPlot
%                 fzri.Visible = 'on';
%             end      
%             % save figure
%             saveas(pzri, sprintf('speed_%1.0d.png', k));
%         end   
%         function captureStress(obj, k)
%             % plots the dimensional stress profile over the entire
%             % domain at the indicated time step, k
%             % compute stress magnitudes and time that will be ploted    
%             ur = obj.uk(k); uz = obj.vk(k);
%             [~, ~, tauR, tauZ] = stress(obj, reshape(ur', [], 1), ...
%                                                  reshape(uz', [], 1));
%             s = sqrt(tauR.^2 + tauZ.^2);                                
%             t = obj.tau2t(obj.tau(k));
%             [R, Z] = meshgrid(obj.rbar, obj.zbar); 
%             % plot contour
%             fzri = figure('Units', 'normalized', ...
%                 'Position', [0 0 0.4*obj.b 0.4], 'Visible', 'off');
%             pzri = surf(R, Z, s);         
%             xlim([0, max(R(:))])
%             ylim([min(Z(:)), max(Z(:))])
%             pbaspect([obj.b, 1, 1]);  % figure sized proportional to aspect ratio
%             view(0, 90);
%             caxis([0, 1]);
%             colormap(flipud(spring));
%             cb = colorbar;
%             cb.Ruler.MinorTick = 'on';
%             set(pzri, 'linestyle', 'none');
%             ylabel(cb, '$\vert \tau\vert$ (Pa)', 'interpreter', 'latex', 'FontSize', 14);
%             xlabel('$r/H$', 'interpreter', 'latex', 'FontSize', 14);
%             ylabel('$z/H$', 'interpreter', 'latex', 'FontSize', 14);
%             title(sprintf('$t$ = %1.0f s', t), 'interpreter', 'latex', ...
%                 'FontSize', 14);
%             if nargin > 2 && ShowPlot
%                 fzri.Visible = 'on';
%             end      
%             % save figure
%             saveas(pzri, sprintf('stress_%1.0d.png', k));
%         end
%         function Main(obj)
%             computeArStar(obj);
%             computeAzStar(obj);
%             computeApStar(obj);
%             for iter_index = 1:size(obj.tau,2)
%                 cur_tau = obj.tau(iter_index); cur_t = obj.tau2t(cur_tau);
%                 computeUStar(obj);
%                 computeu(obj);
%                 patchVelocity(obj,iter_index);
%             end
%             captureSpeed(obj,size(obj.tau,2));
%             captureSpeed(obj,size(obj.tau,2));
%             animateSpeed(obj, 'plotSpeed.gif', 1, 1e10);
%             animateStress(obj, 'plotSpeed.gif', 1, 1e10);
%         end
    end
             
end
