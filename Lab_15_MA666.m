%%  MA666 - Lab 15:  Wavelets, and public data for neuroscientists. 
%   In this tutorial we will have a brief look at the complex Morelet
%   wavelet, as well as public data resources for neuronal data.

% some math for wavelets:
% for the fft
    
%    Xj=sum(Xn* e^-2piiFjtn) from 1 to N discrete timesteps
%   for continuous time
% X(freq)=integral(X(t)e^-piifjtn)) dt

% now wavelet transform:

%X(scale,pos)=integral(Xt)* Y(s,p) over time dt from - infinity to infinity
%   when Y(s,p) is the wavelet function
%   
% essentially compare wavelet function to data at different scales and
% positions


% wavelettransform= 2/sqrt(pi*fb) * e^2pi i Fc X) * e^(-x^2)/fb)
% looks like a little wavelet, like a baby SWR
% this has real and imaginary parts
% Fc is the higher freq dominant rhythm frequency
% the envelope is just a filter put on that freq

% scale and position:
% the idea is move around the wavelet by adjusting the scale/position
% in time and amplitude axes and look for a good fit
%
% adjusting scale:
% as S increases, the width of our little wavelet increases- basically
% multiple x axis by 1 or 2 or .5
% so basically you use the wavelet template to do your math
% the wavelet is basically a good mathematic choosing of good analytic
% frequencies

% in practice:
% with fixed scale
% rake your wavelet across your dataset, compute integral, Xi(S0, K=1:end)
% now repeat for different scale(S)
% output is a 2d array- x(m,n) m=S(scale), n=k(positions)

%%  Preliminaries.
%   Text preceded by a '%' indicates a 'comment'.  This text should appear
%   green on the screen.  I will use comments to explain what we're doing 
%   and to ask you questions.  Also, comments are useful in your own code
%   to note what you've done (so it makes sense when you return to the code
%   in the future).  It's a good habit to *always* comment your code.  I'll
%   try to set a good example, but won't always . . . 

%%  Part 1.  The complex Morlet wavelet.
%   To start, let's have a look at the complex Morelet wavelet.  To do so, 
%   we'll first fix the two parameters we need to specify the wavelet, the
%   bandwidth parameter and center frequency,

for i=1:4
fb=i*4; % bandwidth parameter
fc=.5; % center frequency (like scaling factor)

wname = ['cmor' num2str(fb) '-' num2str(fc)];
figure;
[psi,xval] = wavefun(wname, 10);
plot3(xval,real(psi),imag(psi));
end
%   Note that, in the line above, we define a new variable "wname".
%  The documentation for wavelets is under wavelet families.  this is
%  complex morelet, which is good for a continuous wavelet transform, but
%  for discrete you want to use something like daubechies, or 
%IN LAB Q:  What is "wname"?

%   Now, let's look at the complex Morelet wavelet.  To do so, we'll call
%   the MATLAB function "wavefun",

[psi,xval] = wavefun(wname, 10);

% wname= ['sym2']; % num2str(fb)];
% [phi,psi,xval] = wavefun(wname, 10);

%[F1,F2] = wfilters('haar');

%   The 2nd argument to wavefun specifies the resolution at which we'll
%   evaluate the wavelet.  In our case, we set it to be big enough that
%   the wavelet looks smooth;  note that "wavefun" evaluates the wavelet at
%   2^(2nd argument) points. Remember that the complex Morelet wavelet has
%   both a real and imaginary part, and we need to plot both,
figure;
% plot(xval, real(psi));
% hold on
% plot(xval, imag(psi), 'r')
% hold off
%plot3(xval,phi,psi);
plot3(xval,real(psi),imag(psi));

%IN LAB Q:  Does the wavelet appear as you expect?
% this is your basic wavelet, its a spherical wavelet so that you have
% better resolution
%%  Part 2:  Continuous wavelet transform.
%   With the wavelet now chosen (and examined), let's apply the continuous
%   wavelet transform to some example data.  To do so, we'll return to a
%   data set from a previous lab (Lab 5).  Please download from Blackboard,
%
%load   Ch5_d3.mat

%plot(t,d);
% this uses
%spectrum1=PlotPowerSpectrum(d,t,'hanning',1);

%periodogram(d.*hann(length(d)));

params.Fs=1/(t(2)-t(1));
params.fpass=[5 100];
[Sb,ts,f]=mtspecgramc(d,[1 .01],params);
figure;
imagesc(ts,fliplr(f),mean(Sb,3)'); set(gca,'ydir','normal')
% for i=1:5
%     subplot(1,5,i)
%     dnow=d((length(d)/i)-1:d(length(d)/i))
%     periodogram(dnow)
% end



% this dataset has a rhythm that gets faster by time
%
%IN LAB Q:  Load these data in MATLAB, and describe what you find.  Please
%compute the power spectrum (with a Hann window), and compare it to your
%visual inspection.  Are the two consistent?  Is the power spectrum a
%good representation of these data?

%   Now, let's compute the continuous wavelet transform for these data.  To
%   do so, we first need to chose the scales.  Let's start with these,

scales = (10:.05:80);

%   Then, we compute the continuous wavelet transform, and ask MATLAB to
%   plot the results.  Notice in the code below, we specify the wavelet
%   using the variable "wname".

figure;
X = cwt(d,scales,wname,'plot');
% this is just like doing a fft for alot of short windows, the only
% difference is that if you were to do the fft for time windows, you woud
% have the same time resolution across... but for this you have better time
% resolution for higher frequency rhythms

%X=dwt(d,wname);
%IN LAB Q:  Describe the plotted results.  What do you see?

%   The plot above is useful, but for a more interpretable result, we'd
%   like to display the wavelet transform as a function of *frequency*
%   (instead of scale).  Luckily, there's a MATLAB function to approximate
%   the frequency corresponding to each scale.  To see this apporximation,
%   we can use the MATLAB function "centfreq".

centfrq(wname,10,'plot');

%   The red curve indicates the sinusoid that best matches the frequency of
%   the real and imaginary parts of the wavelet (in blue).  Notice that, in
%   this case, the sinusoid completes 1 cycle over a duration of 2;  so
%   the frequency of the wavelet is approxiamtely 0.5 Hz.
%   NOTE:  Compare this to the parameter fc.

%   Finally, let's plot the modulus of the wavelet transform as a function
%   of time and *frequency*.  To do so in MATLAB, we'll call one more
%   function "scal2frq" which returns the pseudo-frequencies corresponding
%   to the scales.  Notice that, for this function, we also input the
%   sampling period (dt).  Let's try it,

dt=t(2)-t(1);

figure
f = scal2frq(scales, wname, dt);
imagesc(t,f,abs(X))
xlabel('Time [s]')
ylabel('Freq [Hz]')
colormap jet; set(gca,'ydir','normal');


%% IN LAB Q:  Does this look correct?

%   There's an insidious plotting issue occuring here.  Let's look at the
%   variable "f",

plot(f)

%   Notice that "f" is a decreasing vector, and not linear.  To plot the
%   frequencies on the image of "X", we need to be careful about how we
%   specify the vertical axes.  One way to do so is by setting the y-axis
%   tickmarks explicitly.  For example,
subplot(1,2,1);
imagesc(t,scales,real(X))
set(gca, 'YTick', scales(1:2:end))
set(gca, 'YTickLabel', round(f(1:2:end)))
xlabel('Time [s]')
ylabel('Freq [Hz]')
subplot(1,2,2);
imagesc(t,scales,abs(X))
linkaxes;
set(gca, 'YTick', scales(1:2:end))
set(gca, 'YTickLabel', round(f(1:1:end/2)))
xlabel('Time [s]')
ylabel('Freq [Hz]')
colorbar


%IN LAB Q:  Compare the magnitude of the wavelet transform (plotted as a
%function of time and frequency) with the activity in the original signal.
%Are the two consistent?

%%  Part 3:  Some (interesting) data resources.
%   As a computational neuroscientist, it's often useful to have access to
%   brain data.  Sometimes these data are provided by a colleague or
%   collaborator.  In addition, there exist a variety of online resources
%   for neuronal data.  Here's a link to an extensive list of these online
%   resources:
%
%   http://datasharing.incf.org/ep/Resources#Inventory_of_Resources_for_Datasharing_in_Neuroimaging_and_Electrophysiology_Research
%
%   As an example, let's focus on one resource,
%
%    CRCNS.org
%
%   As an example, I've downloaded one data set from the lab of Karel
%   Svoboda.  These data correspond to this paper,
%
%   http://www.ncbi.nlm.nih.gov/pubmed/24361077
%
%   To get a sample of these data, download this file (Week 15 folder) from
%   Blackboard,
%
%   data_structure_NL_example20140905_ANM219037_20131117.mat
%
%   And load it into MATLAB,

clear
load('data_structure_NL_example20140905_ANM219037_20131117.mat')

%IN LAB Q:  Let's look at the variable "obj".  What is it?  Can we make any
%sense of it?

%   Luckily, the investigators have provided much more detailed information
%   about this data structure.  To check it out, visit Blackboard and
%   download the PDF,
%
%   crcns_alm-1_data_description.pdf
%
%   In additionm the investigators have also provided us some example MATLAB
%   code.  Let's try it out!  Please download this file from Blackboard,
%
%   Demo_get_trial_aligned_raster_PSTH.m
%
%   Let's try it out,

Demo_get_trial_aligned_raster_PSTH(obj)

%IN LAB Q:  What do you find?  Can you "see" what's happening in the code?
