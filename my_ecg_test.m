clear all;
%-----------------------------------
X = load ('118e00m_n.mat');
X = X.val;
%-----------------------------------
N = length(X);     %number of samples from original signal
z = nextpow2(N);
p = 2^z;
%update number samples if necessary by padding zeros at the end for a power
%of 2:
if(p > N)
    p = p-N;
    X(end + p)=0;
    N = length(X);
end
%-----------------------------------
%2^7=128, 2^8=256 2^9=512
%-----------------------------------
i_max = z;   %last level possible.  i = 0 is the first level
i=1;
%-----------------------------------
cD_cell = cell(1,1);  %cell unit of 1x1 dimension
cA_cell = cell(1,1);  %cell unit of 1x1 dimension
%-----------------------------------
%Decomposition
cA = X;
test=N;

while (i < i_max +1)
    [cA,cD]=dwt(cA,'coif4','mode','per');   %send current approximation and save cA and cD
    cD_cell{i,1}= cD;
    cA_cell{i,1}= cA;
    test = N/(2^i);
    i=i+1;
end

i=1;
while (i < i_max +1)
    cD_i = cD_cell{i,1};  %current vector of coefficient at level i

    l_d = length(cD_i);
    m_i = median(abs(cD_i));
    sig_i = m_i/0.6745;   %constant for threshold
    T = sig_i*sqrt(2.*log(l_d));
    k=1;
    while(k < l_d+1)
        if( abs(cD_i(k)) < T )
            cD_i(k)=0;
        end
            k=k+1;
    end
    cD_cell{i,1} = cD_i;
    i=i+1;
end
%-----------------------------------
%Reconstruction
i=i_max;
cA = cA_cell{i,1};
while (i > 0)
    cA = idwt(cA,cD_cell{i,1},'coif4','mode','per');
    i=i-1;
end
figure(2)
X = X(1:1024);
plot(X);
grid on;
title('Original Signal');
xlabel('Number of Samples');
ylabel('Amplitude in mV');

figure(3)
cA = cA(1:1024);
plot(cA);
grid on;
title('Denoised Signal');
xlabel('Number of Samples');
ylabel('Amplitude in mV');

%TODO: Refine algorithm for locating peaks..
% find R-peak every first peak

[~,locs_R]=findpeaks(X,'minpeakheight',50,'minpeakdistance',250);
X_inverted = -X;
[~,locs_S]=findpeaks(X_inverted,'minpeakheight',50,'minpeakdistance',250);

%--------------------------------------
[~,locs_R_D]=findpeaks(cA,'minpeakheight',50,'minpeakdistance',250);
cA_inverted = -cA;
[~,locs_S_D]=findpeaks(cA_inverted,'minpeakheight',50,'minpeakdistance',250);

%for detecting Q, from R-peak position find the minimum all the way back and save
%the spot
i= length(locs_R);
j=1;
locs_Q = zeros(length(i));

while(j< i+1)
a=locs_R(j)-35;
if(a<0)
    a=1;
end
temp = X_inverted(a:locs_R(j));   % size properly 

[~,index]=findpeaks(temp,'npeaks',1);
    b = length(temp)- index;
    locs_Q(j)= locs_R(j)-b;
    j=j+1;
end

figure(4)
subplot(2,1,1);
plot(X), hold on;
plot(locs_R,X(locs_R),'rv','MarkerFaceColor','r'), 
plot(locs_S,X(locs_S),'rs','MarkerFaceColor','b'),
plot(locs_Q,X(locs_Q),'rs','MarkerFaceColor','g'),hold off;
title('Original Signal');
xlabel('Number of Samples');
ylabel('Amplitude in mV');
grid on;
legend(' ECG signal','Q-wave','R-wave','S-wave');

i= length(locs_R_D);
j=1;
locs_Q_D = zeros(length(i));
while(j< i+1)
    a=locs_R_D(j)-50;
if(a<0)
    a=1;
end
temp = cA_inverted(a:locs_R_D(j));    

[~,index]=findpeaks(temp,'npeaks',1);
    b = length(temp)- index;
    locs_Q_D(j)= locs_R_D(j)-b;
    j=j+1;
end

subplot(2,1,2);
plot(cA), hold on;
plot(locs_R_D,cA(locs_R_D),'rv','MarkerFaceColor','r'), 
plot(locs_S_D,cA(locs_S_D),'rs','MarkerFaceColor','b'),
plot(locs_Q_D,cA(locs_Q_D),'rs','MarkerFaceColor','g'), hold off;
title('Denoised Signal');
xlabel('Number of Samples');
ylabel('Amplitude in mV');
grid on;
legend('ECG signal','Q-wave','R-wave','S-wave');


figure(10)
subplot(2,1,1);
plot(X)
title('Original Signal');
xlabel('Number of Samples');
ylabel('Amplitude in mV');
grid on;

subplot(2,1,2);
plot(cA)

title('Denoised Signal');
xlabel('Number of Samples');
ylabel('Amplitude in mV');
grid on;
