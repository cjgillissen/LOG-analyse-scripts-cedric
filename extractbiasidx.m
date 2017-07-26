function [Per fMat,BiasMat] = extractbiasidx(log,timeline,ctrials,windowsize)
% calculates nback performance and bias (perfL - perfR)/(perfL+perfR)
BiasMat = nan(length(ctrials),max(cellfun(@length,ctrials)));
PerfMat = nan(length(ctrials),max(cellfun(@length,ctrials)),2); %2 for sides
SideOpt = unique(log.Side);

for cidx = 1:length(ctrials)
    trials = ctrials{cidx};
    trialcount = 1;
    for tidx = trials
        for sididx = 1:length(SideOpt)
            %find nback trials of this side
            done = 0;
            count = 1;
            while ~done
                trialidx = tidx-ceil((windowsize+count)/2):1:tidx+floor((windowsize+count)/2);
                if any(trialidx<1) || any(trialidx>=max(cellfun(@max,ctrials)))
                    done = 1;
                elseif sum(strcmp(log.Side(trialidx),SideOpt(sididx))) == windowsize
                    done = 1;
                end
                count = count+1;                
            end
            trialidx(trialidx<1|trialidx>max(cellfun(@max,ctrials))) = [];
            try
            trialidx(~strcmp(log.Side(trialidx),SideOpt(sididx))) = [];
            catch
                keyboard
            end
            PerfMat(cidx,trialcount,sididx) = sum(strcmp(log.Reaction(trialidx),'Hit'))./sum(strcmp(log.Reaction(trialidx),'Hit')+strcmp(log.Reaction(trialidx),'Error')); 
            
        end
        BiasMat(cidx,trialcount) = (PerfMat(cidx,trialcount,1) - PerfMat(cidx,trialcount,2))./(PerfMat(cidx,trialcount,1) + PerfMat(cidx,trialcount,2));
        trialcount = trialcount+1;

        
    end
end


end
