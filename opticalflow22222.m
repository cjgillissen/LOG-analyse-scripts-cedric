
load('C:\Users\gillissen\Desktop\FREYC2\ThrowAwayIdx.mat')
load('C:\Users\gillissen\Desktop\FREYC2\Frey1_RawData_C2.mat')
load('C:\Users\gillissen\Desktop\FREYC2\BASELINEMAT.mat')


tmp = single(conddata(:,:,:,1:10));
clear conddata

for i = 1:size(tmp,4) %loop over trials 
    
    V1 = tmp(:,:,:,i);
    V1 = imgaussfilt(V1,2); %Gaussian filter
    tmp(:,:,:,i) = V1;            
end



%% Slow trent correction

trialidx = single(ctrials{2});
% 
% for i = 1:100:size(tmp,1)
%      QQ = single(tmp(i:i+99,:,:,:));
%      QQ = QQ ./ permute(repmat(BASELINEMAT(i:i+99,:,trialidx(1:10)),[1,1,1,size(tmp,3)]),[1,2,4,3]);
%      tmp(i:i+99,:,:,:) = QQ;
% end
% 
% clear QQ



for i = 1:100:size(tmp,1)
    QQ = single(tmp(i:i+99,:,:,:));
    base = single(squeeze(nanmean(QQ(:,:,timeline>=-300 & timeline<0,:),3))); %Baseline
    base(base==0) = nan; %Remove 0 and make nan; cannot divide by 0
    QQ = single(QQ(:,:,timeline>=0 &timeline<=2000,:)); %F
    QQ = (QQ - permute(repmat(base,[1,1,1,size(QQ,3)]),[1,2,4,3]))./permute(repmat(base,[1,1,1,size(QQ,3)]),[1,2,4,3]); %dF/F
    trace(i:i+99,:,:,:) = QQ;
end






mask = ones(800,800);
mask(isnan(trace(:,:,5,1)))=0;

trace = squeeze(nanmean(trace,4));
trace = trace(:,:,4:end);
trace = trace(220:680,220:680,:);
trace = trace(1:2:end,1:2:end,:);


save('FreyC2','trace')
save('MaskfreyC2','mask')


for i = 1:size(trace,3)
    
    imagesc(trace(:,:,i));
end
