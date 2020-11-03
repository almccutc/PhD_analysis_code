function [Data,Time]=agw_filter(Data, Time, DataMax, DataMin)
%
% [DataF,TimeF]=agw_filter(Data, Time, DataMax, DataMin)
%
% Function to threshold and adaptive Gaussian filter data that is a
% function of one variable (assumed to be time here but could be one
% spatial coordinate).
%
% DataF - returned matrix of filtered data
% TimeF - returned vecotr of times of filtered data
% Data  - Data matrix to be filtered - assumed to be a function of 1
%         variable and that dependence coincides with the longer dimension of the
%         matrix
% DataMax - Maximum allowable value in the data set
% DataMin - Minimum allowable value in the data set
%
% Written by Edwin A. Cowen III
% Copyright reserved 2006
% Contact Cowen at eac20@cornell.edu before distributing this software.

%fact=normcdf(1,0,1)-normcdf(-1,0,1); %updated below to use student's t,
%now t_fact
[NI, NJ]=size(Data);
%Test of data has time access running in column direction 
%Note: assume that longer axis is time axis
%Take transpose if test not met
if NI>NJ
    Data=Data';
    [NI, NJ]=size(Data);
end
Nstart=NJ;
%Threshold filter the data to lie between DataMin and DataMax
valid=(Data <= DataMax & Data >= DataMin);
if NI>1
  ind=reshape(valid,NI,NJ);
  %Remove data at a given time that does not meet threshold in any component
  valid=(sum(ind)==NI);
else
  ind=valid;
end
%Return filtered data and time of valid data
Data=Data(:,find(valid));
Time=Time(find(valid));
[NI, NJ]=size(Data);
disp(sprintf('Started with %d, filtered to %d, threshold',Nstart,NJ))
% Adaptive Gaussian Filter

coef=[2 1.25 1.15 1 1.06 1.03 1.01];

for i=1:length(coef)
  DOF=NJ-1; %degrees of freedom
  p=1/(2*NJ); %two-sided probability
  filtFact=-tinv(p,DOF); %coefficient of dynamic filter
  t_fact=tinv(0.75,DOF); %coefficient for normalzing width of IQR based on student's-t
  DataMedian=median(Data')';
  %s=std(Data')'; %estimate standard deviation based on std
  s=iqr(Data')'/(2*t_fact); %estimate standard deviation based on IQR
  HiLim=DataMedian+coef(i)*filtFact*s;
  LoLim=DataMedian-coef(i)*filtFact*s;
  HiLim=repmat(HiLim,1,NJ);
  LoLim=repmat(LoLim,1,NJ);
  valid=(Data <= HiLim & Data >= LoLim);
  
  if NI>1
    ind=reshape(valid,NI,NJ);
    %Remove data at a given time that does not meet threshold in any component
    valid=(sum(ind)==NI);
  else
    ind=valid;
  end
  
  %Return filtered data and time of valid data
  Data=Data(:,find(valid));
  Time=Time(find(valid));
  [NI, NJ]=size(Data);
  disp(sprintf('Started with %d, filtered to %d',Nstart,NJ))
end
NJprev=NJ+1;
while NJprev>NJ
  DOF=NJ-1; %degrees of freedom
  p=1/(2*NJ); %two-sided probability
  filtFact=-tinv(p,DOF); %coefficient of dynamic filter
  t_fact=tinv(0.75,DOF); %coefficient for normalzing width of IQR based on student's-t
  DataMedian=mean(Data')';  %note switched to mean in loop
  %s=std(Data')'; %estimate standard deviation based on std
  s=iqr(Data')'/(2*t_fact); %estimate standard deviation based on IQR
  HiLim=DataMedian+coef(i)*filtFact*s;
  LoLim=DataMedian-coef(i)*filtFact*s;
  HiLim=repmat(HiLim,1,NJ);
  LoLim=repmat(LoLim,1,NJ);
  valid=(Data <= HiLim & Data >= LoLim);
  
  if NI>1
    ind=reshape(valid,NI,NJ);
    %Remove data at a given time that does not meet threshold in any component
    valid=(sum(ind)==NI);
  else
    ind=valid;
  end
  
  %Return filtered data and time of valid data
  Data=Data(:,find(valid));
  Time=Time(find(valid));
  NJprev=NJ;
  [NI, NJ]=size(Data);
  disp(sprintf('Started with %d, filtered to %d, in iterative loop',Nstart,NJ))
end
