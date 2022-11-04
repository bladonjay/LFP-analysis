



low = 6;
high = 12;

srate = 1500;
d = designeegfilt(srate,low,high);
kernel = d;
filter.kernel = kernel;
filter.samprate = srate;
filter.descript = [num2str(low),'-',num2str(high),' Hz resp filter'];
filterfile = ['C:\Users\Jadhavlab\Documents\gitRepos\LFP-analysis\JadhavEEGFilter\thetafilter.mat'];
save(filterfile, 'filter');