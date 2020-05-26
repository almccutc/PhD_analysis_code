classdef FFD < handle
    properties
        % list all variables that will be passed between methods in the
        % model here. think of them as global variables within FFD. if a
        % value is included in the definition then it is treated as the
        % default value.
        
        % static model parameters with default values
  
        g = 9.80665                     % (m/s2)
        rhoPack = 2000                  % (kg/m3) particle packed bulk 
                                        % density

        mu = 2.5*1.81e-5                % (kg/m s) dynamic viscosity, u,  of particles 
                                        % moving in air(Bicerano, Douglas 
                                        % and Brune, 1999)
        nu                              % (m2/s) kinematic viscosity of 
                                        % particles moving in air  
        Uinf = 0.01                     % (m/s) normalization velocity
        
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
        end
        function reInitObj(obj)
            % recomputes static variables. should be ran if any of the
            % system properties are changed externally to ensure that 
            % everything is consistent.
            obj.ztop = obj.H;
            obj.rbar = 1e-6:obj.drbar:obj.b;          
            obj.zbar = 0:obj.dzbar:1;
            
            obj.rMaxIndex = size(obj.rbar,2);
            obj.zMaxIndex = size(obj.zbar,2);
            
             % Updated 05-07 JH from KP: storing ghost points
            obj.Urbar = zeros(obj.zMaxIndex+1, obj.rMaxIndex);
            obj.Uzbar = zeros(obj.zMaxIndex, obj.rMaxIndex+1);
            obj.Ubar = [reshape(obj.Urbar', [], 1); ...
                        reshape(obj.Uzbar', [], 1)];
            obj.Pbar  = zeros(obj.zMaxIndex+1,obj.rMaxIndex+1);
            
            obj.tau = 0:obj.dtau:obj.tauEnd;         
            obj.dt = obj.dtau*obj.H/obj.Uinf;
            obj.nu = obj.mu/obj.rhoLoose;  
            obj.Re = obj.H*obj.Uinf/obj.nu;              
            obj.Fr = obj.Uinf/sqrt(obj.g*obj.H);
            
            % matrices that contain velocity values at every time step for
            % all non-ghost points
            obj.u = sparse(kron(eye(length(obj.tau)), ... 
                             ones(length(obj.zbar), ...
                             length(obj.rbar))));
            obj.v = sparse(kron(eye(length(obj.tau)), ... 
                             ones(length(obj.zbar), ...
                             length(obj.rbar))));
        end 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add any additional functions here
        % as a general rule, try to keep the number of operations for 
        % each method at a minimum. if a function is more than ~10 lines
        % then it should probably be split into multiple methods.
        % in general, the methods can be treated as a typical matlab
        % function with the acception of the method for passing in
        % properties. the general form for functions that have an output is
        function x = exampleFunction(obj, input1, input2)
            x = input1*obj.Re + input2*obj.Fr;
        end
        % and the general form for operations is
        function exampleOperation(obj, input1, input2)
            % updates meshed domain
            obj.rbar = 1e-6:obj.drbar:input1;
            obj.zbar = 0:obj.dzbar:input2;
            reInitObj(obj);
        end
        % where the inputs are just arbitrary variables that aren't defined
        % as properties. it is typically a good practice to limit the usage
        % of these and instead just operate with class properties. if a
        % non-property variable is being passed into multiple functions 
        % then it should instead be defined as a property.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function computeUStar(obj)
            % computes the intermediate non-divergence free r velocity
            % component by solving the decoupled r-momentum equation with
            % an implicit-explicit technique.
            n = length(obj.zbar); m = length(obj.rbar);
            n_m = (n+1)*m;
            m_n = (m+1)*n;

            % compute intermediate state matrices
            %   computeArStar(obj);
            %   computeAzStar(obj);

            % compute intermediate advection operators
            computeNr(obj);
            computeNz(obj);
            % compute intermediate r-velocity
            obj.Ustar = [obj.ArStar, zeros(n_m,m_n); zeros(m_n,n_m),...
            obj.AzStar] \[obj.Nr; obj.Nz];


        end
        function R = Fill(obj,vec,ncol,value,rowflag,locator,begin_i,stop_i)
            if rowflag ==true
                begin = (locator-1)*ncol+begin_i;
                stop = (locator-1)*ncol + stop_i;
                vec(begin:stop)=value;
                R = vec;
            else
                begin = (begin_i-1)*ncol + locator;
                stop  = (stop_i-1) * ncol   + locator;
                vec(begin:ncol:stop) = value;
                R = vec;
            end
        end
        function Rfill = FillBoundary(obj,vec,nrow,ncol,...
                                      lval,rval,upval,lowval)
            Rfill = vec;
            if ~isnan(lval) 
                Rfill = Fill(obj,Rfill,ncol,lval,false,1,1,nrow);
            end
            if ~isnan(rval)
                Rfill = Fill(obj,Rfill,ncol,rval,false,ncol,1,nrow);
            end
            if ~isnan(upval)
                Rfill = Fill(obj,Rfill,ncol,upval,true,1,1,ncol);
            end
            if ~isnan(lowval)
                Rfill = Fill(obj,Rfill,ncol,lowval,true,nrow,1,ncol);
            end
        end
        
        function ArStar = SetArBoundaries(obj,Omegas,nr,mr)
            omega1 = Omegas{1};
            omega2 = Omegas{2};
            omega3 = Omegas{3};
            omega4 = Omegas{4};
            omega5 = Omegas{5};
            
            tnmr = (nr+2)*mr;
            
            % left boundary: Ur_ij = 0 (symmetry) Omega1=1.0
            % right boundary: Ur_ij = 0 (wall) Omega1 = 1.0
            omega1 = FillBoundary(obj,omega1,nr,mr,...
                                       1.0,1.0,NaN,NaN);
             
            % full Omega1 vectors to be put in matrix as diagnal
            omega1 = [ones(mr,1);omega1;ones(mr,1)];
        
            
            % left boundary condition: Ur_ij = 0 (symmetry) Omega2 = 0
            % right boundary condition Ur_ij = 0 (wall) Omega2 = 0
            omega2 = FillBoundary(obj,omega2,nr,mr,0,0,NaN,NaN);
            
            % Omega2 for ghost points: all should be zero
            % the 1st upper subdiagonal; first element popped off
            omega2 = [zeros(mr+1,1);omega2;zeros(mr-1,1)];                    
            
            % upper boundary condition Ur_ij + Ur_(i+1)j = 0; Omega3 = 1.0
            Omega3_UpperGhostPoint = ones(mr,1);
            
            % two upper corner ghost points won't have Omega3;
            Omega3_UpperGhostPoint([1,mr])=0.0;
            
                       
            % left and right boundary points: Omega3 = 0
            omega3 = FillBoundary(obj,omega3,nr,mr,...
                                       0.0,0.0,NaN,NaN);
            
            % first (m,1) vectors are to be cut off by spdiags function
            omega3 = [zeros(mr,1);Omega3_UpperGhostPoint;omega3];
                       

            % left and right boundary points: Omega4 = 0  
            omega4 = FillBoundary(obj,omega4,nr,mr,...
                                       0.0,0.0,NaN,NaN);
              
            % last 0.0 is to be cut off by spdiags function
            omega4 = [zeros(mr-1,1); omega4 ;zeros(mr+1,1)];
            
            % outlet num
            OutletNodeNum = sum((obj.rbar)<=obj.a0);
            % lower Boundary conditions: Ur_ij + Ur_(i-1)j = 0 Omega5 = 1.0
            Omega5_LowerGhostPoint = ones(mr,1);
            % for outlet ghost points: Ur_ij - Ur_(i-1)j = 0 (no gradient)
            % Omega5 = -1.0
            Omega5_LowerGhostPoint([1:OutletNodeNum])=-1.0;
            % two lower corner ghost points don't have meaningful eqn
            Omega5_LowerGhostPoint([1,mr])=0.0;
            
            % left and right boundary points: Omega5 = 0
            omega5 = FillBoundary(obj,omega5,nr,mr,...
                                       0.0,0.0,NaN,NaN);
            % last mx1 vectors are to be cut off by spdiags function
            omega5 = [omega5;Omega5_LowerGhostPoint;zeros(mr,1)];
            ArStar = spdiags([omega1, omega2, omega3,...
                         omega4, omega5], ...
                         [0, 1, mr, -1, -mr], tnmr,  tnmr);  
                               
        end
        function computeArStar(obj)
            % evaluates state matrix for intermediate r-velocity comp.
            
            % n the number of z nodes; m the number of r nodes
            n = length(obj.zbar); m = length(obj.rbar); 
            
            % the number of real points m*(n-1)
            n_1m = (n-1)*m;
            
            % Omega1*u_ij + Omega2*u_i(j+1) + Omega3*u_(i+1)j
            %  +Omega4*u_i(j-1) + Omega5*u_(i-1)j = N
            Omega1 = -1/obj.dtau - 2/(obj.Re*obj.drbar^2) ...
                     - 1./(obj.Re*repmat(obj.rbar', n-1, 1).^2) ...
                     - 2/(obj.Re*obj.dzbar^2); 
            Omega2 = 1./(2*obj.Re*repmat(obj.rbar',n-1,1)*obj.drbar)...
                    + 1/(obj.Re*obj.drbar^2);
            Omega3 = (1/(obj.Re*obj.dzbar^2))*ones(n_1m, 1);
            Omega4 = -1./(2*obj.Re*repmat(obj.rbar',n-1,1)*obj.drbar) ...
                    + 1/(obj.Re*obj.drbar^2);
            Omega5 = (1/(obj.Re*obj.dzbar^2))*ones(n_1m, 1);
            Omegas = {Omega1,Omega2,Omega3,Omega4,Omega5};
                        
            obj.ArStar = SetArBoundaries(obj,Omegas,n-1,m);    
            
        end
        
        
        function AzStar = SetAzBoundaries(obj,Omegas,nz,mz)
            omega1 = Omegas{1};
            omega2 = Omegas{2};
            omega3 = Omegas{3};
            omega4 = Omegas{4};
            omega5 = Omegas{5};
            
            tmnz = nz*mz;
            
            %Finds the nodes of the outlet 
            OutletNodeNum = sum((obj.rbar)<=obj.a0);
            
            % all boundary Omega1 = 1.0
            omega1 = FillBoundary(obj,omega1,nz,mz,...
                                      1.0,1.0,1.0,1.0);
            
            % left boundary except for two corner ghost points
            % left boundary condition: zero gradient 
            % Omega2 = -1.0
            % Omega2 for two left corner ghost points = 0
            % Upper boundary omega2 = 0
            % Lower boundary condition Omega2 = 0
            % Right boundary condition Omega2 = 0
            
            omega2 = FillBoundary(obj,omega2,nz,mz,...
                                      -1.0,0.0,0.0,0.0);
            
            % Last element doesn't have Omega2
            omega2(end)=[];
            
            % first zero is to be cut off
            omega2 = [0.0;omega2];
            
            % upper boundary omega3 = 0
            % left boundary Omega3 = 0
            % right boundary Omega3 = 0.0                                         
            
            omega3 = FillBoundary(obj,omega3,nz,mz,...
                                     0.0,0.0,0.0,0.0);
            
            % last row of nodes (lower boundary) don't have Omega3
            omega3(mz*(nz-1)+1:end)=[];
            % zeros are to be cut off
            omega3 = [zeros(mz,1);omega3];
           
            % left boundary Omega4 = 0
            % lower boundary Omega4 = 0           
            % right boundary Omega4 = 1.0
            % upper boundary Omega4 = 0.0
            omega4 = FillBoundary(obj,omega4,nz,mz,...
                                      0.0,1.0,0.0,0.0);
            
            
            % first element doesn't have Omega4
            omega4(1)=[];
            %zeros to be cut off
            omega4 = [omega4;0.0];
            
            % left boundary Omega5 = 0.0
            % right boundary Omega5 = 0.0
            % lower boundary outlet condition: zero gradient
            % U_ij = U_(i-1)j  Omega5 = -1.0
            % lower boundary outlet condition: wall
            % U_ij = 0; Omega5 = 0.0
            omega5 = FillBoundary(obj,omega5,nz,mz,...
                                      0.0,0.0,0.0,NaN);
            omega5 = Fill(obj,omega5,mz,-1.0,true,nz,2,OutletNodeNum);                      
            omega5(mz*(nz-1)+OutletNodeNum+1:end)=0.0;
            omega5(mz*(nz-1)+1)=0.0;
            % first (m+1) elements don't have Omega5
            omega5(1:mz)=[];
            omega5 = [omega5;zeros(mz,1)];
            
            AzStar = spdiags([omega1,omega2,...
                         omega3, omega4, ...
                         omega5],[0, 1, mz, -1, -mz],tmnz, tmnz); 
        end
        
        function computeAzStar(obj)
            % evaluates state matrix for intermediate r-velocity comp.
            
            % n the number of z nodes; m the number of r nodes
            % the number of real points (m-1)*n
            n = length(obj.zbar); m = length(obj.rbar); m_1n = n*(m-1);

            % the number of ghost points 2n
            TotalNum = n*(m+1);
            
            % internal real points Omega1
            Omega1 = (-1/obj.dtau - 2/(obj.Re*obj.drbar^2) ...
                     - 2/(obj.Re*obj.dzbar^2))*ones(TotalNum,1);
            
            % internal real points Omega2
            Omega2 = 1./(2*obj.Re*(repmat(([-obj.drbar,...
                     obj.rbar])',n, 1) ...
                     + obj.drbar/2)*obj.drbar) + 1/(obj.Re*obj.drbar^2);
            
            % real internal points Omega3
            Omega3 = (1/(obj.Re*obj.dzbar^2))*ones(TotalNum,1);
            
            % real internal points Omega4
            Omega4 = -1./(2*obj.Re*(repmat(([-obj.drbar,...
                     obj.rbar])',n, 1) ...
                     + obj.drbar/2)*obj.drbar) + 1/(obj.Re*obj.drbar^2);
            
            % real internal points Omega5
            Omega5 = 1/(obj.Re*obj.dzbar^2)*ones(TotalNum,1);
            
            
            Omegas = {Omega1,Omega2,Omega3,Omega4,Omega5};
            
            obj.AzStar = SetAzBoundaries(obj,Omegas,n,m+1);                     
        end
        function computeNr(obj)
            % computes advection operator for intermediate
            % non-divergence-free r-momentum equation
            n = size(obj.Urbar, 1); m = size(obj.Urbar, 2); nm = n*m;
            
            % initialize Nr if not already
            if isempty(obj.Nr)
                obj.Nr = NaN*ones(nm, 1);
            end   
            
            % match uz cells with ur
            [ur, uz, ~] = obj.matchCells('ur');                        
                        
            % compute diagonal elements of state matrices
            Lambda1 = 1/(2*obj.drbar); Lambda2 = -1/(2*obj.drbar);
            Pi1 = 1/(2*obj.dzbar); Pi2 = -1/(2*obj.dzbar); 
            
            % compile sparse diagonal state matrices
            obj.ArN = spdiags([Lambda1*ones(nm, 1), ...
                               Lambda2*ones(nm, 1)], [1, -1], nm, nm);              
            obj.BrN = spdiags([Pi1*ones(nm, 1), Pi2*ones(nm, 1)], ...
                              [m, -m], nm, nm);             

            % compute advection operator
            ur = reshape(ur', [], 1);
            uz = reshape(uz', [], 1);
            obj.Nr = diag(ur)*(obj.ArN*ur - ones(nm, 1)./obj.dtau) ...
                + diag(uz)*obj.BrN*ur;              
                                    
            % replace boundary elements
            setNrBoundaries(obj);                         
        end
                      
        function setNrBoundaries(obj)
            % replaces boundary elements in Nr according to prescribed
            % boundary conditions
            m = size(obj.Urbar, 2);

            % boundary 1 at z = 0 (free surface)
            obj.Nr(1:m-1) = 1; %0;
            
            % boundary 2 at r = 0 (centerline)
            obj.Nr(1:m:end-m) = 2; %0;
            
            % boundary 3 at z = 1 (bin opening)
            [~, ia] = min(abs(obj.a0 - obj.rbar));
            obj.Nr(end-m:end-m+ia) = 3; %0;
            
            % boundary 4 at z = 1 (bottom wall)
            obj.Nr(end-m+ia+1:end) = 4; %0;
            
            % boundary 5 at r = b (outer radius wall)
            obj.Nr(m:m:end) = 5; %0;                                   
        end
       
        function computeNz(obj)
            % computes advection operator for intermediate
            % non-divergence-free z-momentum equation
            n = size(obj.Uzbar, 1); m = size(obj.Uzbar, 2); nm = n*m;
            
            % initialize Nr if not already
            if isempty(obj.Nz)
                obj.Nz = NaN*ones(nm, 1);
            end   
            
            % match ur cells with uz
            [ur, uz, ~] = obj.matchCells('uz');
                        
            % compute diagonal elements of state matrices
            Lambda1 = 1/(2*obj.drbar); Lambda2 = -1/(2*obj.drbar);
            Pi1 = 1/(2*obj.dzbar); Pi2 = -1/(2*obj.dzbar); 
            
            % compile sparse diagonal state matrices
            obj.AzN = spdiags([Lambda1*ones(nm, 1), ...
                               Lambda2*ones(nm, 1)], [1, -1], nm, nm);              
            obj.BzN = spdiags([Pi1*ones(nm, 1), Pi2*ones(nm, 1)], ...
                              [m, -m], nm, nm);                           

            % compute advection operator
            ur = reshape(ur', [], 1);
            uz = reshape(uz', [], 1);
            obj.Nz = diag(ur)*obj.AzN*uz ...
                   + diag(uz)*(obj.BzN*uz - ones(nm, 1)./obj.dtau) ...
                   - ones(nm, 1)./obj.Fr^2;               
                                    
            % replace boundary elements
            setNzBoundaries(obj);            
        end
                      
        function setNzBoundaries(obj)
            % replaces boundary elements in Nr according to prescribed
            % boundary conditions
            m = size(obj.Uzbar, 2);

            % boundary 1 at z = 1 (free surface)
            obj.Nz(1:m-1) = 1; %obj.Uinf;
            
            % boundary 2 at r = 0 (centerline)
            obj.Nz(1:m:end-m) = 2; %0;
            
            % boundary 3 at z = 0 (bin opening)
            [~, ia] = min(abs(obj.a0 - obj.rbar));
            obj.Nz(end-m:end-m+ia) = 3; %0;
            
            % boundary 4 at z = 0 (bottom wall)
            obj.Nz(end-m+ia+1:end) = 4; %0;
            
            % boundary 5 at r = b (outer radius wall)
            obj.Nz(m:m:end) = 5; %0;                                   
        end    
       % function computeNz(obj)
       %     % computes advection operator for intermediate
       %     % non-divergence-free r-momentum equation
       %     n = length(obj.zbar); m = length(obj.rbar); nm = n*m;
       %     
       %     % compute diagonal elements of state matrices
       %     Lambda1 = 1/(2*obj.drbar); Lambda2 = -1/(2*obj.drbar);
       %     Pi1 = 1/(2*obj.dzbar); Pi2 = -1/(2*obj.dzbar); 
       %     
       %     % compile sparse diagonal state matrices
       %     obj.AzN = spdiags([Lambda1*ones(nm, 1), ...
       %                        Lambda2*ones(nm, 1)], [1, -1], nm, nm);              
       %     obj.BzN = spdiags([Pi1*ones(nm, 1), Pi2*ones(nm, 1)], ...
       %                       [m, -m], nm, nm);                           
       %     R = spdiags(0.25*ones(nm, 4), [0, -m, -m+1, 1], nm, nm);

       %     % compute advection operator
       %     ur = obj.Ubar(1:nm);
       %     uz = obj.Ubar(nm+1:end);
       %     obj.Nz = diag(R*ur)*(obj.ArN*ur - ones(nm, 1)./obj.dtau) ...
       %            + diag(uz)*(obj.BrN*ur - ones(nm, 1)./obj.dtau) ...
       %            - ones(nm, 1)./obj.Fr^2;                                 
       % end
        
        
        function [ur, uz, p] = matchCells(obj, base_variable)
            % matches the cells for primitive variables to that of the base
            % variable (specified by 'ur', 'uz', or 'p')
            nr = size(obj.Urbar, 1); mr = size(obj.Urbar, 2); nmr = nr*mr;
            nz = size(obj.Uzbar, 1); mz = size(obj.Uzbar, 2);            
            
            % compute averaging matrices
            [Zr, Zp] = averageUz(obj);
            [Rz, Rp] = averageUr(obj);
            [Pr, Pz] = averageP(obj);
            
            % shift and align to grid for base primitive variable
            switch base_variable
                case 'ur'
                    % snap uz to ur grid
                    ur = obj.Urbar;
                    uz = reshape(Zr*reshape(obj.Uzbar', [], 1), mz, nz)';
                    uz = [NaN*ones(1, mr); uz(:, 2:end)];
                    p = reshape(Pr*reshape(obj.Pbar', [], 1), mz, nr)';
                    p = p(:, 2:end);
                case 'uz'
                    % snap ur to uz grid
                    ur = reshape(Rz*reshape(obj.Urbar', [], 1), mr, nr)';
                    ur = [NaN*ones(nz, 1), ur(2:end, :)];
                    uz = obj.Uzbar;
                    p = reshape(Pz*reshape(obj.Pbar', [], 1), mz, nr)';
                    p = p(2:end, :);
                case 'p'
                    % snap ur and uz to p grid
                    ur = reshape(Rp*reshape(obj.Ustar(1:nmr)', [], 1), ...
                                                                 mr, nr)';
                    ur = [NaN*ones(nr, 1), ur];
                    uz = reshape(Zp*reshape(obj.Ustar(nmr+1:end)', ...
                                                         [], 1), mz, nz)';
                    uz = [NaN*ones(1, mz); uz];
                    p = obj.Pbar;
            end                                                                               
        end
        
        function [Zr, Zp] = averageUz(obj)
            % computes averaging matrix for ur and p nodes
            nz = size(obj.Uzbar, 1); mz = size(obj.Uzbar, 2); nmz = nz*mz;
            
            % averaging matrix for 4 uz nodes surrounding ur        
            Zr = spdiags(0.25*ones(nmz, 4), [0, mz, mz-1, -1], nmz, nmz);
            
            % averaging matrix for 2 uz nodes above/below p
            Zp = spdiags(0.5*ones(nmz, 2), [0, mz], nmz, nmz);                       
            
            % correct for z = 0 boundary
            Zr(1:mz, :) = 0;           
            Zr(1:mz, 1:mz) = 0.5*eye(mz);
            Zr(2:nmz+1:(nmz+1)*(mz-1)) = 0.5;   
            Zp(1:mz, :) = 0;
            Zp(1:mz, 1:mz) = eye(mz); 
            
            % correct for z = 1 boundary            
            Zr(end-mz+1:end, :) = 0;                  
            Zr(end-mz+1:end, end-mz+1:end) = 0.5*eye(mz);
            Zr(end-(mz-1)*(nmz+1)+1:nmz+1:end-nmz) = 0.5;
            Zp(end-mz+1:end, :) = 0;
            Zp(end-mz+1:end, end-mz+1:end) = eye(mz);            
        end
        
        function [Rz, Rp] = averageUr(obj)
            % computes averaging matrix for uz and p nodes
            nr = size(obj.Urbar, 1); mr = size(obj.Urbar, 2); nmr = nr*mr;
            
            % averaging matrix for 4 ur nodes surrounding uz
            Rz = spdiags(0.25*ones(nmr, 4), [0, -mr, -mr+1, 1], nmr, nmr);
            
            % averaging matrix for 2 ur nodes left/right of p
            Rp = spdiags(0.5*ones(nmr, 2), [0, 1], nmr, nmr); 
            
            % correct for r = 0 boundary            
            Rz(1:mr:end-mr, :) = 0;            
            Rz(1:(nmr+1)*mr:end-nmr*mr) = 0.5;            
            Rz((nmr+1)*(mr-1)+1:(nmr+1)*mr:end) = 0.5;
            Rp(1:mr:end-mr, :) = 0;    
            Rp(1:(nmr+1)*mr:end-nmr*mr) = 1;
            
            % correct for r = b boundary
            Rz(mr:mr:end, :) = 0;
            Rz(nmr*mr+1:(nmr+1)*mr:end) = 0.5;
            Rz(nmr*(2*mr-1)+mr:(nmr+1)*mr:end) = 0.5;                        
            Rp(mr:mr:end, :) = 0;            
            Rp((nmr+1)*(mr-1)+1:(nmr+1)*mr:end) = 1;           
        end
        
        function [Pr, Pz] = averageP(obj)
            % computes averaging matrix for uz and p nodes
            n = size(obj.Pbar, 1); m = size(obj.Pbar, 2); nm = n*m;
            
            % averaging matrix for 4 p nodes surrounding ur
            Pr = spdiags(0.5*ones(nm, 2), [0, -1], nm, nm); 
            
            % averaging matrix for 4 p nodes surrounding uz
            Pz = spdiags(0.5*ones(nm, 2), [0, -m], nm, nm);                                
        end
       
        function AzStar = SetApBoundaries(obj,Omegas,np,mp)
            omega1 = Omegas{1};
            omega2 = Omegas{2};
            omega3 = Omegas{3};
            omega4 = Omegas{4};
            omega5 = Omegas{5};
            
            tmnp = np*mp;
            
            OutletNodeNum = sum((obj.rbar)<=obj.a0);
            
            % all boundary Omega1 = 1.0
            omega1 = FillBoundary(obj,omega1,np,mp,...
                                      1.0,1.0,1.0,1.0);
            
            % left boundary except for 2 corner ghost points: zero gradient 
            % Omega2 = -1.0
            % Omega2 for two left corner ghost points = 0
            % Upper boundary omega2 = 0
            % Lower boundary condition Omega2 = 0
            % Right boundary condition Omega2 = 0
            
            omega2 = FillBoundary(obj,omega2,np,mp,...
                                      -1.0,0.0,0.0,0.0);
            
            % Last element doesn't have Omega2
            omega2(end)=[];
            
            % first zero is to be cut off
            omega2 = [0.0;omega2];
            
            % upper boundary except for 2 corner ghost points:zero gradient 
            % Omega3 = -1.0
            % Omega3 for two left corner ghost points = 0
            % Omega3 for the first valid ghost points = 0 (specify p = 0)
            % left  boundary Omega3 = 0   
            % right boundary Omega3 = 0 
            % lower boundary Omega3 = 0
                                     
            
            omega3 = FillBoundary(obj,omega3,np,mp,...
                                     0.0,0.0,-1.0,0.0);
            omega3([1,mp]) = 0.0;
            
            % last row of nodes (lower boundary) don't have Omega3
            omega3(mp*(np-1)+1:end)=[];
            % zeros are to be cut off
            omega3 = [zeros(mp,1);omega3];
            
            
            % right boundary except for 2 corner ghost points:zero gradient 
            % Omega4 = -1.0
            % Omega4 for two left corner ghost points = 0
            % left  boundary Omega4 = 0   
            % upper boundary Omega4 = 0 
            % lower boundary Omega4 = 0
            omega4 = FillBoundary(obj,omega4,np,mp,...
                                      0.0,-1.0,0.0,0.0);
            
            
            % first element doesn't have Omega4
            omega4(1)=[];
            %zeros to be cut off
            omega4 = [omega4;0.0];

            % lower boundary except for 2 corner ghost points:zero gradient 
            % Omega5 = -1.0
            % Omega5 for two left corner ghost points = 0
            % Omega5 for the lower right corner (-1) = 0 
            % left  boundary Omega5 = 0   
            % right boundary Omega5 = 0 
            % upper boundary Omega5 = 0
            omega5 = FillBoundary(obj,omega5,np,mp,...
                                      0.0,0.0,0.0,-1.0);
            omega5([(np-1)*mp+1,tmnp-1,tmnp]) = 0.0;
            % first (m+1) elements don't have Omega5
            omega5(1:mp)=[];
            omega5 = [omega5;zeros(mp,1)];
            
            AzStar = spdiags([omega1,omega2,...
                         omega3, omega4, ...
                         omega5],[0, 1, mp, -1, -mp],tmnp, tmnp); 
        end
        
        
        function computeApStar(obj)
            % evaluates state matrix for n+1 pressure comp.
            
            % n the number of z nodes; m the number of r nodes
            % the number of real pressure points (m-1)*(n-1)
            n = length(obj.zbar); m = length(obj.rbar); 
            mn_1 = (m-1)*(n-1);

            % the number of total points for boundary settings
            TotalNum = (n+1)*(m+1);
            
            % internal real points Omega1
            Omega1 = (-2.0/(obj.drbar^2)-2.0/(obj.dzbar^2))...
                      *ones(TotalNum,1);
            
            % internal real points Omega2
            Omega2 = 1./(2*(repmat(([-obj.drbar,...
                     obj.rbar])',n+1, 1) ...
                     + obj.drbar/2)*obj.drbar) + 1/(obj.drbar^2);
            
            % real internal points Omega3
            Omega3 = (1/(obj.dzbar^2))*ones(TotalNum,1);
            
            % real internal points Omega4
            Omega4 = -1./(2*(repmat(([-obj.drbar,...
                     obj.rbar])',n, 1) ...
                     + obj.drbar/2)*obj.drbar) + 1/(obj.drbar^2);
            
            % real internal points Omega5
            Omega5 = 1/(obj.Re*obj.dzbar^2)*ones(TotalNum,1);
            
            
            Omegas = {Omega1,Omega2,Omega3,Omega4,Omega5};
            
            obj.ApStar = SetApBoundaries(obj,Omegas,n+1,m+1);                     
        end
        
        function computeDstar(obj)
            % Computes velocity operator for the discretized poison
            % equation
            n = size(obj.Pbar, 1); m = size(obj.Pbar, 2); nm = n*m;
            
            % compute Ustar if not yet populated
            if isempty(obj.Ustar)
                obj.Ustar = ones(length(obj.Urbar(:))+length(obj.Uzbar(:)), 1);
%                 computeUstar(obj);
            end
            
            % match urstar and uzstar cells with p
            [ur, uz, ~] = obj.matchCells('p');
            r = [NaN, obj.rbar];
            
            % compute diagonal elements of state matrices
            Lambda1 = 1/(2*obj.dtau*(repmat(r', n, 1) + obj.drbar/2)) ...
                    - 1/(obj.dtau*obj.drbar); 
            Lambda2 = 1/(2*obj.dtau*(repmat(r', n, 1) + obj.drbar/2)) ...
                    + 1/(obj.dtau*obj.drbar); 
            Pi1 = -1/(obj.dtau*obj.dzbar); Pi2 = 1/(obj.dtau*obj.dzbar); 
            
            % compile sparse diagonal state matrices
            Dr = spdiags([Lambda1', Lambda2'], [0, 1], nm, nm);              
            Dz = spdiags([Pi1*ones(nm, 1), Pi2*ones(nm, 1)], ...
                              [0, m], nm, nm);                           

            % compute advection operator
            ur = reshape(ur', [], 1);
            uz = reshape(uz', [], 1);
            obj.DStar = [Dr, Dz]*[ur; uz];                          
            
            % replace boundary elements
            setDstarBoundaries(obj);                                     
        end
        
        function setDstarBoundaries(obj)
            % sets boundary conditions for RHS of pressure equation (Dstar)
            m = size(obj.Pbar, 2);           
            
            % boundary 1 at z = 1 (free surface)
            obj.DStar(1:m-1) = 0;%1; %0; (updated JH 05-13)
            
            % boundary 2 at r = 0 (centerline)
            obj.DStar(1:m:end-m) = 0;%2; %0;
            
            % boundary 3 at z = 0 (bin opening)
            [~, ia] = min(abs(obj.a0 - obj.rbar));
            obj.DStar(end-m:end-m+ia) = 0;%3; %0;
            
            % boundary 4 at z = 0 (bottom wall)
            obj.DStar(end-m+ia+1:end) = 0;%4 %0;
            
            % boundary 5 at r = b (outer radius wall)
            obj.DStar(m:m:end) = 0; %5;      
            
            % prescribed boundary point for controling matrix rank
            obj.DStar(end-1) = 0; %6;
        end
        
        function computepressure(obj)
           n = size(obj.Pbar, 1); m = size(obj.Pbar, 2); nm = n*m;
            obj.Pbar = obj.ApStar\obj.DStar;
    
        end
        
        
        function computeu(obj)
            
            n = size(obj.Urbar, 1); m = size(obj.Uzbar, 2); nm = n*m;
            n_1m = (n-1)*m;
            m_1n = (m-1)*n;
            
            %computeApStar(obj);
            computeDstar(obj);
            setDstarBoundaries(obj);
            
            obj.Pbar = [obj.ApStar]\[obj.DStar]; 
            
            
            [~, ~, pr] = obj.matchCells('ur');
            [~, ~, pz] = obj.matchCells('uz');
            
            pr_ = reshape(pr', [], 1);
            pz_ = reshape(pz', [], 1);
            
            Dr = obj.diffR(n,m-1);
            Dz = obj.diffZ(n-1,m);
            
            obj.Urbar = (reshape(obj.Ustar(1:m_1n)-Dr*pr_,m-1,n))';
            obj.Uzbar = (reshape(obj.Ustar(m_1n+1:end)-Dz*pz_,m,n-1))';
 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % first derivative operators (KP)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function D = diffR(obj, n, m)
            % creates a differential operator for first derivative in r 
            % dimension with central differencing at internal points 
            % and forward/backward differencing at boundaries
            nm = n*m;
            
            % set central differencing for center points
            D = spdiags([ones(nm, 1), -ones(nm, 1)], [1, -1], nm, nm);
            D = D./(2*obj.drbar);
            
            % delete boundary entries
            D(m:m:end, :) = 0;      % r = 0
            D(m+1:m:end, :) = 0;    % r = b
            
            % set forward differencing for left boundary
            D(1:m*(nm+1):end) = -1/obj.drbar;
            D(nm+1:m*(nm+1):end) = 1/obj.drbar;
            
            % set backward differencing for right boundary
            D((m-1)*(nm+1)+1:m*(nm+1):end) = 1/obj.drbar;
            D((m-2)*(nm+1)+2:m*(nm+1):end) = -1/obj.drbar; 
        end
        
        function D = diffZ(obj, n, m)
            % creates a differential operator for first derivative in z 
            % dimension with central differencing at internal points 
            % and forward/backward differencing at boundaries
            nm = n*m;
            
            % set central differencing for center points
            D = spdiags([ones(nm, 1), -ones(nm, 1)], [m, -m], nm, nm);
            D = D./(2*obj.dzbar);
            
            % delete boundary entries
            D(1:m, :) = 0;          % z = 0
            D(end-m+1:end, :) = 0;    % z = 1
            
            % set forward differencing for top boundary
            D(1:nm+1:m*(nm+1)) = -1/obj.dzbar;
            D(m*nm+1:nm+1:2*m*nm+1) = 1/obj.dzbar;
            
            % set backward differencing for bottom boundary
            D(end-(m-1)*(nm+1):nm+1:end) = 1/obj.dzbar;
            D(end-2*m*nm-m:nm+1:end-m*nm) = -1/obj.dzbar;               
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % post processing and conversions (KP)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function t = tau2t(obj, tau)
            % converts non-dimensional time, tau, to dimensional time, t
            t = tau*obj.H/obj.Uinf;
        end
        
        function patchVelocity(obj, k_)
            % combines thetaS, thetaT and thetaC for time step k_
            n = length(obj.zbar); m = length(obj.rbar);
            % overlay current r-velocity
            obj.u(n*(k_-1)+1:n*k_, m*(k_-1)+1:m*k_) ...
                = obj.Urbar(2:end, :);
            % overlay current r-velocity
            obj.v(n*(k_-1)+1:n*k_, m*(k_-1)+1:m*k_) ...
                = obj.Uzbar(:, 2:end);
        end
        
        function x = uk(obj, k)
            % returns the r-velocity at every point in the mesh for time k
            n = length(obj.zbar); m = length(obj.rbar);
            x = obj.u(n*(k-1)+1:n*k, m*(k-1)+1:m*k);                       
        end
        
        function x = vk(obj, k)
            % returns the z-velocity at every point in the mesh for time k
            n = length(obj.zbar); m = length(obj.rbar);
            x = obj.v(n*(k-1)+1:n*k, m*(k-1)+1:m*k);                       
        end
        
        function u = uMag(obj, ur, uz)
            % computes the velocity magnitude for input vectors of the same
            % size (ghost points need to be removed prior)
            u = sqrt(ur.^2 + uz.^2);           
        end
        
        function [Srr, Szz, tauR, tauZ] = stress(obj, ur, uz)
            % computes strain rates and viscous stress at every point in
            % given velocity vectors
            n = length(obj.zbar); m = length(obj.rbar);
            
            % compute derivative operators
            Dr = obj.diffR(n, m);
            Dz = obj.diffZ(n, m);
            
            % compute strain rate with differential operators
            Srr = Dr*ur;
            Szz = Dz*uz;            
            
            % compute viscous stress
            tauR = 2*obj.mu*Srr*obj.Uinf/obj.H;
            tauZ = 2*obj.mu*Szz*obj.Uinf/obj.H; 
            
            % reshape to physical domain
            Srr = reshape(Srr, m, n)';
            Szz = reshape(Szz, m, n)';
            tauR = reshape(tauR, m, n)';
            tauZ = reshape(tauZ, m, n)';
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotting (KP)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function animateSpeed(obj, plot_filename, rr, snapK)
            % generates gif that displays velocity magnitude at every time
            % step
            if nargin < 4
                snapK = 1e10;
            end
            if nargin < 3
                rr = 0.5;
            end
            if nargin < 2
                plot_filename = 'speedPlot.gif';
            end
            figure('Units', 'normalized', ...
            'Position', [0 0 0.4*obj.b 0.4], 'Visible', 'off');            
            [~, tfigs] = plotZR(obj,  obj.uMag(obj.uk(1), obj.vk(1)), true);
            ylabel(colorbar, '$\vert u/U_\infty\vert$', 'interpreter', 'latex', ...
                'FontSize', 14); 
            caxis([0, sqrt(2)]);
            gif(plot_filename, 'frame', gcf);
            % update sequentially
            for ii = 1:length(obj.tau)                                      
                gif;
                t = obj.tau2t(obj.tau(ii));
                u_ = obj.uMag(obj.uk(ii), obj.vk(ii));
                tfigs.ZData = u_; 
                title(sprintf('$t$ = %1.0f s', t), 'interpreter', 'latex', ...
                'FontSize', 14);
                pause(rr);
                % save plot if on save time-step
                if mod(ii, snapK) == 0
                    captureSpeed(obj, ii);
                end
            end            
        end
        
        function animateStress(obj, plot_filename, rr, snapK)
            % generates gif that displays velocity magnitude at every time
            % step
            if nargin < 4
                snapK = 1e10;
            end
            if nargin < 3
                rr = 0.5;
            end
            if nargin < 2
                plot_filename = 'stressPlot.gif';
            end
            ur = obj.uk(1); uz = obj.vk(1);
            [~, ~, tauR, tauZ] = stress(obj, reshape(ur', [], 1), ...
                                                 reshape(uz', [], 1));
            s = sqrt(tauR.^2 + tauZ.^2);                                          
            figure('Units', 'normalized', ...
            'Position', [0 0 0.4*obj.b 0.4], 'Visible', 'off');            
            [~, tfigs] = plotZR(obj, s, true);
            ylabel(colorbar, '$\vert \tau \vert$ (Pa)', 'interpreter', 'latex', ...
                'FontSize', 14); 
            caxis([min(s(:)), max(s(:))]);
            gif(plot_filename, 'frame', gcf);
            % update sequentially
            for ii = 1:length(obj.tau)                                      
                gif;
                ur = obj.uk(ii); uz = obj.vk(ii);
                [~, ~, tauR, tauZ] = stress(obj, reshape(ur', [], 1), ...
                                                 reshape(uz', [], 1));
                s = sqrt(tauR.^2 + tauZ.^2); 
                t = obj.tau2t(obj.tau(ii));
                tfigs.ZData = s; 
                title(sprintf('$t$ = %1.0f s', t), 'interpreter', 'latex', ...
                'FontSize', 14);
                pause(rr);
                % save plot if on save time-step
                if mod(ii, snapK) == 0
                    captureStress(obj, ii);
                end
            end            
        end
        
        function [fzri, pzri] = plotZR(obj, x, ShowPlot)
            % plots ZR plane for one time instance
            [R, Z] = meshgrid(obj.rbar, obj.zbar); 
            % plot contour
            fzri = figure('Units', 'normalized', ...
                'Position', [0 0 0.4*obj.b 0.4], 'Visible', 'off');
            pzri = surf(R, Z, x);         
            xlim([0, max(R(:))])
            ylim([min(Z(:)), max(Z(:))])
            pbaspect([obj.b, 1, 1]);  % figure sized proportional to aspect ratio
            view(0, 90);
            caxis([0, 1]);
            colormap(flipud(spring));
            cb = colorbar;
            cb.Ruler.MinorTick = 'on';
            set(pzri, 'linestyle', 'none');
            xlabel('$r/H$', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$z/H$', 'interpreter', 'latex', 'FontSize', 14);
            if nargin > 2 && ShowPlot
                fzri.Visible = 'on';
            end
        end            
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
