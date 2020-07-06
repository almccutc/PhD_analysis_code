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
 
        u_original
        v_original
        datamax
        datamin
        original_percentleft
        medianfilter_percentleft
        agwpercentleft
        u_agw
        v_agw
        u_nanfilter
        v_nanfilter
        u_mean
        v_mean
        u_f
        v_f
        u_rms
        v_rms
        target
        
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
            load('G:\J20matfiles\032420_AM_0.6_40_105_630_10min_4.5V_Correct.mat','u_original','v_original') %Change% 
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
        function allfunctions(obj)
            reshapes(obj);
            checkHistogram(obj);
            applyAGWfilter(obj);
            medfilter(obj);
            velocityCalculations(obj);
            %spatialspectra(obj);
            %temporalspectra(obj);
            
        end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        function reshapes(obj)
         u_o = cell2mat(obj.u_original); %No filter is done within PIVLAB
         v_o = cell2mat(obj.v_original); %These values have NaNs in them

        % Rotate Matrix
         nlay = length(obj.u_original);
         [r,c] = size(u_o);
        % 
         obj.u_original = permute(reshape(u_o',[c,r/nlay,nlay]),[2,1,3]);
         obj.v_original = permute(reshape(v_o',[c,r/nlay,nlay]),[2,1,3]);
         
         %matches NaN values for both u and w
         obj.u_original(isnan(obj.v_original)) = NaN;
         obj.v_original(isnan(obj.u_original)) = NaN;
         
        end
        function checkHistogram(obj)
            [Ny,Nx,Nt] = size(obj.u_original);
            total = Ny*Nx*Nt; 
            logicarray = ~isnan(obj.u_original);
            sumofnonNaN = sum(sum(sum(logicarray)));
            
            obj.original_percentleft = sumofnonNaN/total;

            ucheck = reshape(obj.u_original,1,Nx*Ny*Nt); vcheck = reshape(obj.v_original,1,Nx*Ny*Nt);

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

            % Percentage of velocities left after the AGW filter has been applied, I think Blair said this should be 90 percent or higher
            % if the data is good
            [Ny,Nx,Nt] = size(obj.u_original);
            total = Ny*Nx*Nt; 
            logicarray = ~isnan(obj.u_agw);
            sumofnonNaN = sum(sum(sum(logicarray)));
            
            obj.agwpercentleft = sumofnonNaN/total; %it is the same for all components
            
            %f = msgbox(num2str(obj.agwpercentleft*100), 'Percent Remaining. Ideally 95% or more.')

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
            hold off
            
            pause;
            
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function medfilter(obj)
            obj.target = 1; %initial guess 
            
            [Nt,Ny,Nx] = size(obj.u_agw);
            
            while(1)
     
            m=1;
            
            for tt=1:Nt
                close all

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
                
                uclean_save(tt,:,:) = utclean;
                vclean_save(tt,:,:) = vtclean;
                

                figure(2); %shows original data
                uoriginal = squeeze(obj.u_original(tt,:,:));
                voriginal = squeeze(obj.v_original(tt,:,:));
                quiver(uoriginal,voriginal,2,'b');
                xlim([0 Nx])
                ylim([0 Ny])
                set(gcf,'Position',[700 300 800 700])
                hold on;
                
                %shows what agw filter removes
                uagw = squeeze(obj.u_agw(tt,:,:));
                vagw = squeeze(obj.v_agw(tt,:,:));
                quiver(uagw,vagw,2,'r');
                xlim([0 Nx])
                ylim([0 Ny])
                hold on;
                
                %the values after agw and medfilter
                quiver(utclean,vtclean,2,'k'); 
                xlim([0 Nx])
                ylim([0 Ny])
                
                if m==3
                    continue
                    else    
                    m = menu('Yes if to continue through time, No for new target value, All to apply filter at all time steps, Exit to stop program.','Yes','No', 'All', 'Exit');
                end 
                
                if m==2  % yes stored as 1, no stored as 2, all has a value of 3, exit is 4
                    break;
                end
                
                if m==4
                    break;
                end    

                
            end
            
            if m==4
                break;
            end 
            if m==3
                break;
            end 
            
            %check if the target should be updated then, depending on how
            %the plots looking 
            obj.target = str2num(cell2mat(inputdlg('Enter new target:',...
            'Target', [1 50])));
            
            end
            obj.u_nanfilter = uclean_save;
            obj.v_nanfilter = vclean_save;
            
            [Ny,Nx,Nt] = size(obj.u_original);
            total = Ny*Nx*Nt; 
            logicarray = ~isnan(obj.u_nanfilter);
            sumofnonNaN = sum(sum(sum(logicarray)));
            obj.medianfilter_percentleft = sumofnonNaN/total; 
            
            
            close all
            
            beep
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function velocityCalculations(obj)
            
            % Bring in data, the is for a 3D double array, mine is 53 (height) by 79 (length) by 10,500 (in time)
            u_o = obj.u_nanfilter; v_o = obj.v_nanfilter;
            
            u_o=permute(u_o,[1 3 2]); v_o=permute(v_o,[1 3 2]); 
            u_o=permute(u_o,[2 3 1]); v_o=permute(v_o,[2 3 1]);
            
            obj.u_mean = nanmean(u_o,3); obj.v_mean = nanmean(v_o,3);
            obj. u_f = u_o-obj.u_mean; obj.v_f = v_o-obj.v_mean;
            obj.u_rms = sqrt(nanmean((obj.u_f.^2),3)); obj.v_rms = sqrt(nanmean((obj.v_f.^2),3));
            obj.tke = 0.5*(2*(obj. u_f.^2) + (obj. v_f.^2));
            obj.tke = nanmean(obj.tke,3);
            
            %Time average for each subwindow                            
            figure (1) 
            subplot(3,2,1)
            imagesc(obj.u_mean)
            colorbar
            caxis([-0.02,0.02])
            title('Colorbar U-Velocity (m/s)','Interpreter','Latex')
            subplot(2,2,2)
            imagesc(obj.v_mean)
            colorbar
            caxis([-0.02,0.02])
            title('Colorbar W-Velocity (m/s)','Interpreter','Latex')

            %Time rms average for each subwindow
            subplot(2,2,3)
            imagesc(obj.u_rms)
            colorbar
            %caxis([-0.02,0.02])
            title('Colorbar rms (m/s)')
            subplot(2,2,4)
            imagesc(obj.v_rms)
            colorbar
            %caxis([-0.02,0.02])
            title('Colorbar rms (m/s)')

            figure (2)
            imagesc(obj.tke)
            colorbar
            %caxis([-0.02,0.02])
            title('TKE (m^2/s^2)')
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function spatialspectra(obj)
            %53 by 79 subwindows
            Suutime=[];
            a = 35; %x direction subwindow strip
            b= 25; %y direction subwindow strip
            calibration = s; %pixel to cm conversion
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
        function temporalspectra(obj)
            u_o = obj.u_nanfilter; v_o = obj.v_nanfilter;
            
            u_o(isnan(u_o)) = 0; v_o(isnan(v_o)) = 0;

            %436(left to right), 288(top going dowwards)
            %so approx subwindow (18,27)
            xlocation = 18;
            ylocation = 27;

            u = u_o;
            v = v_o;
            [Ny,Nx,Nt] = size(u);
            f_s = 115; %sample frequency
            N = Nt/10; 
            %N=1150;
            T = N/f_s;
            df = f_s/N;
            f = df/2:df:f_s-df/2;
            fny = f(1:(length(f)/2));

            tempspectraU = [];
            tempspectraW = [];

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
