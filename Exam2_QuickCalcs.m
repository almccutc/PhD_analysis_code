clear all

%% A conventional least squares fit of the model y=bt to these data 
%(regression of y on t) gives the following value for parameter b closest 
%to?

%y=bx then do y(1/b) = x
x= [-1,-2,-3,4,0,2]'; %this was y
y = [0,1,2,3,4,5.]'; %this was time
b = y\x % this works

%%
%A conventional least squares fit of the model y=bt^2  to these data gives a 
%value for parameter b closest to
x= [1;1;-1;-1]; %this was y
y = [0;2;4;6]; %this was time
b = y.^2\x % this works
%%
%A conventional least squares fit of the model y=c3t gives a value for c 
%closest to,
x= [1;1;-1;-1]; %this was y
y = [0;2;4;6]; %this was time
b = 3.^y\x % this works

%%
%  If standard deviations of errors in the 4 y data are [1,2,1 ,1], the 
%estimated value of parameter b via weighted least squares fit to the model 
%y=bt is closest to

v= [1,2,1 ,1];
Cd = diag(v);

x= [1;1;-1;-1]; %this was y
y = [0;2;4;6]; %this was time
 
m = ((y'*Cd^-1*y)^-1)*y'*(Cd^-1)*x %this might work? Go with closest value
%%
x= [2,4,6,8,9];
x=diag(x);
y = [1;2;3;4;5];

b = x\y
%%
G = [1;1;1];

ans1=((transpose(G)*G)^-1)*transpose(G);

cd = diag([1/2;1/2;1/2]);
var = [2;2;2];

ans2 = (1/(transpose(G)*cd*G))*G.*cd

%% final Q's Prob 5
%the least squares estimate of parameter ?a? for  yt ~ayt-1 
%can't figure this out
%y=bx 
x= [-1,-2,-3,4,0];
x=diag(x);
yt=[-2,-3,4,0,2];

b = x\yt'
%%
%What is the autocorrelation of the time series [-3,3,0,-3,3]
a = [-3,3,0,-3,3];
b = xcorr(a)

%% padded series dft

a = [zeros(1,100), rand(1,100)];
S = fft(a)

%% transient convolution 
g = [1,2,3];
h = [4,5];
w = conv(g,h)

%%
x = [2,-2,1,-1];
means = mean(x)
w=1;
V = var(x,w)

%%
 gt = [-2, 8,-1,-1];
 ht = [1,-1,8,-1];
 
ans= ifft(fft(gt).*fft(ht))

%%
U = [1,1,1];
V = [2,3,4]';
ans = U*V;
ans2 = V*U

%%
ans1 = fft([1,1,3,2,6,7])
ans2 = fft([1,1,3,2,6,7,0,0,0,0,0,0])
%diff=ans1-ans2

