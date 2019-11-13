function LFPstruct=LFP_NEXfile(dirName,outDir)
% function LFP_NEXfile(dirName,outDir)
% extracts only LFP data from a nex file, asks you which nex file you want
% to use.

% Extract LFP traces

% First I have to turn all my lfps into nex files.
% once I get the nex files, I have to reconstitute the timestamps and match
% them with the traces.  Then im probably going to have to find the best
% LFP from the dataset.  I think i'll also have to restrict my dataset to
% the first two jumpers (max tetrode is 61) because alot of the MEC
% recordings also have RSP, and LEC recordings have PER


if ~exist('dirName','var')
    dirName=uigetdir('','Choose Directory of NEX files');
elseif isempty(dirName)
    dirName=uigetdir('','Choose Directory of NEX files');
end

if ~exist('outDir','var')
    outDir=uigetdir('','Choose a Save Directory');
elseif isempty(dirName)
    outDir=uigetdir('','Choose a Save Directory');
end
dash=WhichDash;

dirData=dir(dirName);
fileNames={dirData(~[dirData.isdir]).name}';
fileList=cellfun(@(a) fullfile(dirName,a), fileNames,'UniformOutput',false);
choose=0;
% so you can choose multiple files
if choose==1
    if length(fileNames)<20
        checked=checkBox(fileNames);
    else
        checked1=checkBox(fileNames(1:20));
        check   ed2=checkBox(fileNames(21:end));
        checked=[checked1 checked2+20];
    end
else
    checked=ones(length(fileNames),1);
end
% now go through all files and find your LFP
for j=1:length(fileList)
    if checked(j)==1
        nex=readNexFileM(fileList{j});
        
        allfields=fieldnames(nex);
        % look for my contvars
        
        % get all the contvars
        LFPfield=strcmp(allfields,'contvars');
        if any(LFPfield)
            allLFPs=nex.(allfields{LFPfield});
        end
        
        
        % ditch the strobe and non lfp channels
        keeplfp=[];
        for i=1:length(allLFPs)
            keeplfp(i)=any(strfind(allLFPs{i}.name,'F'));
        end
        allLFPs=allLFPs(logical(keeplfp));
        
        % now we can go after real lfps:
        % for first ten
        
        figure
        q=0;
        % label and plot all the lfps
        for i=1:min([10 length(allLFPs)]);
            q=q+1;
            subplot(5,2,i);
            plot(allLFPs{i}.data(10000:50000));
            title(['lfp ' num2str(i) ' of ' num2str(length(allLFPs))]);
            
        end
        linkaxes;
        ylim([-.7 .7]);
        
        % if you want to choose a single LFP for this session
        chooselfp=0;
        if chooselfp
            % choose, if you dont choose go for next ten
            bestlfp=input('Choose the best lfp number (1:10):  ');
            close
            if isempty(bestlfp)
                figure
                for i=1:10
                    q=q+1;
                    subplot(5,2,i);
                    plot(allLFPs{q}.data(10000:50000));
                    title(['lfp ' num2str(q) ' of ' num2str(length(allLFPs))]);
                end
                linkaxes;
                ylim([-.7 .7]);
                bestlfp=input('Choose the best lfp number (11:20):  ');
            end
            fprintf('\n');
            if ~isempty(bestlfp)
                LFPstruct=allLFPs{bestlfp};
                save([outDir dash fileNames{j}(1:end-4) '_LFP'],'LFPstruct');
                close
                fprintf(['saved in ' outDir]);
            else
                fprintf('you didnt choose, so nothing saved \n');
            end
        else
            LFPstruct=allLFPs;
        end
    
    end
end





