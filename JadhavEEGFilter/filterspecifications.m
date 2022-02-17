
% Shantanu - Starting with code from mcarr for filter specifications

srate = 1500;
d = designeegfilt(srate,8,18);
spindlefilter = d;
filterfile = sprintf('/home/shantanu/Src/Matlab/Filters/spindlefilter.mat');
save(filterfile, 'spindlefilter');




%Here is where the settings for all filters are specified

%Design filters for Cross-Frequency Coherence analysis for amplitude
%Saved as cfcampfilt*
srate = 1500;
for f = 20:2:200
    d =designeegfilt(srate,f-2,f+2);
    cfcampfilt = d;
    filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/cfcampfilt%s.mat',num2str(f));
	    save(filterfile, 'cfcampfilt');
end

%% Design theta filter
srate = 1500;
d = designeegfilt(srate,5,9);
thetafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/thetafilter.mat');
save(filterfile, 'thetafilter');

%% Design beta-gamma filter

srate = 1500;
d = designeegfilt(srate,10,50);
betagammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/betagammafilter.mat');
save(filterfile, 'betagammafilter');


%% Design beta filter

srate = 1500;
d = designeegfilt(srate,15,30);
betagammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/betafilter.mat');
save(filterfile, 'betagammafilter');


%% Design broadband gamma filter
srate = 1500;
d = designeegfilt(srate,25,150);
gammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/broadbandgammafilter.mat');
save(filterfile, 'gammafilter');

%Design broadband gamma filter for higher srate
srate = 3000;
d = designeegfilt(srate,25,150);
gammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/broadbandgammafilter3000.mat');
save(filterfile, 'gammafilter');

srate = 4500;
d = designeegfilt(srate,25,150);
gammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/broadbandgammafilter4500.mat');
save(filterfile, 'gammafilter');

%% Design ripple filter for slow ripples
% Based on Csicsvari J, Hirase H, Czurko A, Mamiya A, & Buzsaki G (1999)
% Fast network oscillations in the hippocampal CA1 region of the behaving
% rat. J Neuroscience 19:RC20(1-4).
srate = 1500;
d = designeegfilt(srate,100,130);
slowripplefilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/slowripplefilter.mat');
save(filterfile,'slowripplefilter');

%% Design animal specific low and high gamma filters
srate = 1500;

%Bond
lowgamma = [20 50];
d = designeegfilt(srate,lowgamma(1),lowgamma(2));
lowgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/bonlowgammafilter.mat');
save(filterfile,'lowgammafilter');

highgamma = [60 115];
d = designeegfilt(srate,highgamma(1),highgamma(2));
highgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/bonhighgammafilter.mat');
save(filterfile,'highgammafilter');

%Conley
lowgamma = [25 55];
d = designeegfilt(srate,lowgamma(1),lowgamma(2));
lowgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/conlowgammafilter.mat');
save(filterfile,'lowgammafilter');

highgamma = [65 140];
d = designeegfilt(srate,highgamma(1),highgamma(2));
highgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/conhighgammafilter.mat');
save(filterfile,'highgammafilter');

%Corriander
lowgamma = [25 55];
d = designeegfilt(srate,lowgamma(1),lowgamma(2));
lowgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/Corlowgammafilter.mat');
save(filterfile,'lowgammafilter');

highgamma = [60 110];
d = designeegfilt(srate,highgamma(1),highgamma(2));
highgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/Corhighgammafilter.mat');
save(filterfile,'highgammafilter');

%Dudley
lowgamma = [30 55];
d = designeegfilt(srate,lowgamma(1),lowgamma(2));
lowgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/dudlowgammafilter.mat');
save(filterfile,'lowgammafilter');

highgamma = [60 100];
d = designeegfilt(srate,highgamma(1),highgamma(2));
highgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/dudhighgammafilter.mat');
save(filterfile,'highgammafilter');

%Eight
lowgamma = [25 50];
d = designeegfilt(srate,lowgamma(1),lowgamma(2));
lowgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/Eiglowgammafilter.mat');
save(filterfile,'lowgammafilter');

highgamma = [55 110];
d = designeegfilt(srate,highgamma(1),highgamma(2));
highgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/Eighighgammafilter.mat');
save(filterfile,'highgammafilter');

%Five
lowgamma = [25 50];
d = designeegfilt(srate,lowgamma(1),lowgamma(2));
lowgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/Fivlowgammafilter.mat');
save(filterfile,'lowgammafilter');

highgamma = [55 110];
d = designeegfilt(srate,highgamma(1),highgamma(2));
highgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/Fivhighgammafilter.mat');
save(filterfile,'highgammafilter');

%Frank
lowgamma = [20 50];
d = designeegfilt(srate,lowgamma(1),lowgamma(2));
lowgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/fralowgammafilter.mat');
save(filterfile,'lowgammafilter');

highgamma = [60 110];
d = designeegfilt(srate,highgamma(1),highgamma(2));
highgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/frahighgammafilter.mat');
save(filterfile,'highgammafilter');

%Miles
lowgamma = [25 55];
d = designeegfilt(srate,lowgamma(1),lowgamma(2));
lowgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/millowgammafilter.mat');
save(filterfile,'lowgammafilter');

highgamma = [65 110];
d = designeegfilt(srate,highgamma(1),highgamma(2));
highgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/milhighgammafilter.mat');
save(filterfile,'highgammafilter');

%Seven
lowgamma = [25 55];
d = designeegfilt(srate,lowgamma(1),lowgamma(2));
lowgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/Sevlowgammafilter.mat');
save(filterfile,'lowgammafilter');

highgamma = [70 130];
d = designeegfilt(srate,lowgamma(1),lowgamma(2));
highgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/Sevhighgammafilter.mat');
save(filterfile,'highgammafilter');

%Six
lowgamma = [25 55];
d = designeegfilt(srate,lowgamma(1),lowgamma(2));
lowgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/Sixlowgammafilter.mat');
save(filterfile,'lowgammafilter');

highgamma = [65 110];
d = designeegfilt(srate,lowgamma(1),lowgamma(2));
highgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/Sixhighgammafilter.mat');
save(filterfile,'highgammafilter');

%Ten
lowgamma = [20 50];
d = designeegfilt(srate,lowgamma(1),lowgamma(2));
lowgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/tenlowgammafilter.mat');
save(filterfile,'lowgammafilter');

highgamma = [65 110];
d = designeegfilt(srate,highgamma(1),highgamma(2));
highgammafilter = d;
filterfile = sprintf('/home/mcarr/Src/Matlab/Filters/tenhighgammafilter.mat');
save(filterfile,'highgammafilter');
