

            deltax = 16*0.0000247; %m        %distance between velocity vectors
            obj.nu = 0.000000977; %m2/s
            
            [Ny,Nx,Nt] = size(obj.u_nanfilter); %Ny-# of rows, Nx-# of columns
            
            %zero matrices
            z_u = zeros(Ny,1,Nt); 
            z_v = zeros(1,Nx,Nt);
            
            %size adjusted matrices
            u2adjx = [obj.u_nanfilter z_u];
            v2adjy = [obj.v_nanfilter; z_v];
            u2adjy = [obj.u_nanfilter; z_v];
            v2adjx = [obj.v_nanfilter z_u];
            u1adjx = [z_u obj.u_nanfilter];
            v1adjy = [z_v; obj.v_nanfilter];
            u1adjy = [z_v; obj.u_nanfilter];
            v1adjx = [z_u obj.v_nanfilter];
            
            %velocity gradients
            obj.ugradx = u2adjx - u1adjx;
            obj.vgrady = v2adjy - v1adjy;
            obj.ugrady = u2adjy - u1adjy;
            obj.vgradx = v2adjx - v1adjx;
            
            %gradient terms, squared and time averad
            obj.ugrad_termx = nanmean(((obj.ugradx(2:Ny,2:Nx,:)./deltax).^2),3); %time average
            obj.vgrad_termy =nanmean(((obj.vgrady(2:Ny,2:Nx,:)./deltax).^2),3); %time average
            obj.ugrad_termy = nanmean(((obj.ugrady(2:Ny,2:Nx,:)./deltax).^2),3); %time average
            obj.vgrad_termx = nanmean(((obj.vgradx(2:Ny,2:Nx,:)./deltax).^2),3); %time average
            obj.uxvy_term = nanmean((((obj.ugradx(2:Ny,2:Nx,:)./deltax).*(obj.vgrady(2:Ny,2:Nx,:)./deltax))),3); %time average
            obj.uyvx_term = nanmean((((obj.ugrady(2:Ny,2:Nx,:)./deltax).*(obj.vgradx(2:Ny,2:Nx,:)./deltax))),3); %time average
            obj.uy_plus_vx_term = nanmean((((obj.ugrady(2:Ny,2:Nx,:)./deltax)+(obj.vgradx(2:Ny,2:Nx,:)./deltax))),3); %time average
            % dissipation rate time average, m2/s3
            obj.epsilon = obj.nu.*(4.*obj.ugrad_termx+obj.ugrad_termy+obj.vgrad_termx+2.*obj.vgrad_termy+2.*obj.uxvy_term+2.*obj.uyvx_term);
            % dissipation rate spatial average
            obj.epsilon_avg = nanmean(obj.epsilon(:));          
            
            obj.tau_kt = (obj.nu/obj.epsilon_avg)^0.5; % time (s)
            obj.eta_kl = (obj.nu^3/obj.epsilon_avg)^0.25; %length (m)
            
            %integrated dissipation spectrum
            R = (deltax/obj.eta_kl); 
            x_pos = 2*pi/R
            
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
            cvalue = str2num(get(e,'String'));
            close all           
            
            %corrected values
            obj.epsilon_corrected = obj.epsilon*(2-cvalue/100);
            obj.epsilon_avg_corrected = nanmean(obj.epsilon_corrected(:)); 
            obj.tau_kt_corrected = (obj.nu/obj.epsilon_avg_corrected)^0.5; % time (s)
            obj.eta_kl_corrected = (obj.nu^3/obj.epsilon_avg_corrected)^0.25; %length (cm)
            
            figure (1)
            imagesc(obj.epsilon_corrected)
            colorbar
            %caxis([-0.02,0.02])
            title(['Dissipation, Avg: ',num2str(obj.epsilon_avg_corrected,5), ' m^2/s^3'])
%             % using x position to go up and across
%             figure(2)
%             plot(f,Gain, '-.b') % dissipation spectra
%             hold on
%             plot(f,G2, '--r') % integrated dissipation spectrum
%             xlabel('k_\eta')
%             ylabel('D(k)/u^3_\eta, \epsilon(0,k)/\epsilon_m')
%             title('Normalized dissipation spectrum and cumulative dissipation')
%             legend('Normalized dissipation spectrum','Cumulative dissipation')

              pause
              close all
            
%% HIT turbulence dissipation equation    

            obj.epsilon_iso = obj.nu.*(4.*obj.ugrad_termx+2.*obj.vgrad_termy+2.*obj.uy_plus_vx_term);
            % dissipation rate spatial average
            obj.epsilon_iso_avg = nanmean(obj.epsilon_iso(:));

            obj.tau_kt_iso = (obj.nu/obj.epsilon_avg)^0.5; % time (s)
            obj.eta_kl_iso = (obj.nu^3/obj.epsilon_avg)^0.25; %length (m)
            
            %integrated dissipation spectrum
            R = (deltax/obj.eta_kl); 
            x_pos = 2*pi/R
            
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
            cvalue = str2num(get(e,'String'));
            close all           
            
            %corrected values
            obj.epsilon_iso_corrected = obj.epsilon*(2-cvalue/100);
            obj.epsilon_avg_iso_corrected = nanmean(obj.epsilon_iso_corrected(:)); 
            obj.tau_kt_iso_corrected = (obj.nu/obj.epsilon_avg_iso_corrected)^0.5; % time (s)
            obj.eta_kl_iso_corrected = (obj.nu^3/obj.epsilon_avg_iso_corrected)^0.25; %length (cm)
            
            figure (1)
            imagesc(obj.epsilon_iso_corrected)
            colorbar
            %caxis([-0.02,0.02])
            title(['Dissipation, Iso equation, Avg: ',num2str(obj.epsilon_avg_iso_corrected,5), ' m^2/s^3'])

            