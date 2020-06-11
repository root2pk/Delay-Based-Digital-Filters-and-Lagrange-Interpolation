%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%    PROGRAM : VIBRATO 
%    
%    PRROGRAM TO TAKE IN AN INPUT AUDIO AND PERFORM A 
%    VIBRATO(PITCH FLUCTUATION) ON IT. THIS IS DONE BY IMPLEMENTING A 
%    MODULATING DELAY LINE USING AN LFO.
%    THIS USES THE DOPPLER EFFECT TO PRODUCE THE PITCH CHANGE BY 
%    CONSTANTLY CHANGING THE DISTANCE BETWEEN THE OBSERVER AND SOURCE.
%    
%    AUTHOR : RUTHU PREM KUMAR
%    DATE : 06/12/2019
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%read in our WAV file, and store sample rate in Fs
[x,Fs]=audioread('Guitar_dry.wav');

%in case of stereo, to mono
x = 0.5*sum(x,2);


f = 9;                         %frequency = 9 Khz  (5 - 14 Hz)
D = 0.0003;                    %Delay time = 0.3 Milliseconds  


D=round(D*Fs);                 % Delay value in samples

f=f/Fs;                        % Frequency value in samples
L=length(x);                   % length of input audio
del_length=2+D+D*2;            % Length of delay line  
delayline=zeros(del_length,1); % Creating empty vector for delay line
y=zeros(size(x));              % Creating empty vector for output signal


%Computing Output signal
for n=1:L
    
   %Delay LFO
   s = sin(2*pi*f*n);
   %Z picks out the delay sample point
   Z = 1+D+D*s;
   % Since Z isnt always an integer, the value must be estimated using
   % Linear interpolation
   Z_int = floor(Z);
   %Finding the fractional difference between Z and Z_int
   error = Z - Z_int; 
   %Computing Delay line
   delayline=[x(n);delayline(1:del_length-1)]; 
   %Linear Interpolation
   y(n)=delayline(Z_int+1)*error+delayline(Z_int)*(1-error);
   
end 

%Plotting the graphs
figure;

subplot(2,1,1)
plot([1:length(x)]/Fs,x,'r');
title('Input Signal')
subplot(2,1,2);
plot([1:length(y)]/Fs,y,'b');
title('Output Signal');

%Playing the output audio file
soundsc(y,Fs);
