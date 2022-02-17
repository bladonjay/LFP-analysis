low = 7;
high = 8;

srate = 1500;
d = designeegfilt(srate,low,high);
kernel = d;
respfilter.kernel = kernel;
respfilter.samprate = srate;
respfilter.descript = [num2str(low),'-',num2str(high),' Hz resp filter'];
filterfile = ['C:\Users\Jadhavlab\Documents\gitRepos\Jadhav-Lab-Codes\Claires Work\respfilter.mat'];
save(filterfile, 'respfilter');