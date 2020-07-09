function [u_int_temp,v_int_temp] = fillNaNsTemp(U,V,U_UL,U_UR,U_LL,U_LR, V_UL, V_UR, V_LL, V_LR);
%inputs U and V as instantaneous FOV's of displacements
%inputs U_UL, U_UR, etc. are the nanmedian values of corners so that NaN
%corners can be replaced with "good" corner values, since
%triscatteredinterp does not extrapolate, only interpolates

%January 8, 2014

% old version, which takes ages:
% %     %put ufnans into a matrix like todd uses, t,y,x,u,v,#
% %     
% %     count=1;
% %     for tt=1:time
% %         tt
% %         for row=1:height
% %             for col=1:width
% %                 if isnan(ufnan(tt,row,col))==0
% %                     datamat(count,1)=tt;
% %                     datamat(count,2)=row;
% %                     datamat(count,3)=col;
% %                     datamat(count,4)=ufnan(tt,row,col);
% %                     datamat(count,5)=vfnan(tt,row,col);
% %                     datamat(count,6)=count;
% %                     count=count+1;
% %                 end
% %             end
% %         end
% %     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%build fake U matrix for testing

% U(1,1)=11;
% U(2,1)=21;
% U(3,1)=NaN;
% U(1,2)=12;
% U(2,2)=22;
% U(3,2)=32;
% U(1,3)=13;
% U(2,3)=23;
% U(3,3)=33;
% U(1,4)=14;
% U(2,4)=NaN;
% U(3,4)=34;
% U(1,5)=15;
% U(2,5)=25;
% U(3,5)=35;
% 
% U
% 
[Y,X]=size(U);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%replace corner NaNs with median values of existing corner data
%corner data values are computed outside of the code and input in as vars

% % %upper left
% % U_UL=nanmedian(U(1,1));
% % V_UL=nanmedian(V(1,1));
% % %upper right
% % U_UR=nanmedian(U(1,X));
% % V_UR=nanmedian(V(1,X));
% % %lower left
% % U_LL=nanmedian(U(Y,1));
% % V_LL=nanmedian(V(Y,1));
% % %lower right
% % U_LR=nanmedian(U(Y,X));
% % V_LR=nanmedian(V(Y,X));

if isnan(U(1,1))==1;
    U(1,1)=U_UL;
    V(1,1)=V_UL;
end

if isnan(U(1,X))==1;
    U(1,X)=U_UR;
    V(1,X)=V_UR;
end

if isnan(U(Y,1))==1;
    U(Y,1)=U_LL;
    V(Y,1)=V_LL;
end

if isnan(U(Y,X))==1;
    U(Y,X)=U_LR;
    V(Y,X)=V_LR;
end

% U;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%convert 2D snapshot to 4-column listing

[Y,X]=size(U);
blank=zeros(X*Y,5);
countvector=1:1:X*Y;
blank(:,5)=countvector;

    for y=1:Y
        blank(1+X*(y-1) : X+X*(y-1) , 1) = y;
        blank((y-1)*X+1 : (y-1)*X+X , 2) = 1:X;

        blank((y-1)*X+1 : (y-1)*X+X , 3) = U(y,1:X);
        blank((y-1)*X+1 : (y-1)*X+X , 4) = V(y,1:X);
    end

% blanksave=blank; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%delete NaN rows of U

[Lb Wb]=size(blank);
nindex=find(isnan(blank(:,3)));

while length(nindex) > 0
    nindex_temp=nindex(1);
    btop=blank(1:nindex-1,1:Wb);
    bbot=blank(nindex+1:Lb,1:Wb);
    bnew=[btop;bbot];
    
    blank=bnew;
    [Lb Wb]=size(blank);
    nindex=find(isnan(blank(:,4)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%do the Delauney interpolation to replace existing non-corner NaNs

    fx=blank(:, 2);
    fy=blank(:, 1);
    fu=blank(:, 3);
    fv=blank(:, 4);

    Fu = TriScatteredInterp(fx,fy,fu);
    Fv = TriScatteredInterp(fx,fy,fv);

    [qx,qy] = meshgrid(1:X,1:Y);
    qut = Fu(qx,qy);
    qvt = Fv(qx,qy);
    u_int_temp(:,:)=qut(:,:);
    v_int_temp(:,:)=qvt(:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%