
load('C:\Users\gillissen\Desktop\FREYC2\ThrowAwayIdx.mat')
load('C:\Users\gillissen\Desktop\FREYC2\Frey1_RawData_C2.mat')
load('C:\Users\gillissen\Desktop\FREYC2\BASELINEMAT.mat')


tracethiscond = single(conddata(:,:,:,1:10));
clear conddata

p=2; q=2;
tmp = zeros(size(tracethiscond,1)/p,size(tracethiscond,2)/q,size(tracethiscond,3),size(tracethiscond,4));
for i  = 1:size(tmp,4)
    for j = 1:size(tmp,3)
        
        M = tracethiscond(:,:,j,i);
        [m,n]=size(M); %M is the original matrix
        M=sum( reshape(M,p,[]) ,1 );
        M=reshape(M,m/p,[]).'; %Note transpose
        M=sum( reshape(M,q,[]) ,1);
        M=reshape(M,n/q,[]).'; %Note transpose
        tmp(:,:,j,i) = M;
        
    end
end



for i = 1:size(tmp,4) %loop over trials 
    
    V1 = tmp(:,:,:,i);
    V1 = imgaussfilt(V1,2.5); %Gaussian filter
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

clear QQ

I = imagesc(trace(:,:,5,5));
BW = roipoly(I);



mask = ones(800,800);
mask(isnan(trace(:,:,5,1)))=0;

trace = squeeze(nanmean(trace,4));
trace = trace(:,:,4:end);
trace = trace(220:680,220:680,:);
trace = trace(1:2:end,1:2:end,:);

% normalize between 0 and 1
% don't  normalize use binning and then  *100% df/f 
% mintmp = nanmin(trace);
% maxtmp = nanmax(trace);
% difftmp = maxtmp-mintmp;
% tmp = (trace-repmat(mintmp,size(trace,1),1))./repmat(difftmp,size(trace,1),1);

save('FreyC2','tmp')
save('MaskfreyC2','mask')


for i = 1:size(trace,3)
    
    imagesc(trace(:,:,i));
end

%% plot quiver
% with interactive cool ui slider that can slide trough the time...

% [x,y] = meshgrid(1:size(uvCLG,1),1:size(uvCLG,2)); 

i = 1;
f = figure(1);
p = quiver(imag(uvCLG(:,:,i)),real(uvCLG(:,:,i)));

 %// initialize the slider
    h = uicontrol(...
        'parent'  , f,...        
        'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
        'style'   , 'slider',...        
        'position', [0.05 0.05 0.9 0.05],...
        'min'     , 1,...               %// Make the A between 1...
        'max'     , size(uvCLG,3),...              %// and 10, with initial value
        'value'   , i,...               %// as set above.
        'callback', @sliderCallback);   %// This is called when using the arrows
                                  
   hLstn = handle.listener(h,'ActionEvent',@sliderCallback);                                
                               
   function sliderCallback(~,~)
        delete(p);
        p = quiver(imag(uvCLG(:,:,get(h,'value'))),real(uvCLG(:,:,(get(h,'value')))))
        axis tight
        axis([0 2*pi -10 10])
   end 


                
%% try whole brain now
% with mask!

































