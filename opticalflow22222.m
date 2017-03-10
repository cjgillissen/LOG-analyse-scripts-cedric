
load('C:\Users\gillissen\Desktop\Optical flow data\Frey1216\ThrowAwayIdx.mat')
load('C:\Users\gillissen\Desktop\Optical flow data\Frey1216\Frey1_RawData_C2.mat')
% load('C:\Users\gillissen\Desktop\Optical flow data\FREYC2\BASELINEMAT.mat')
load('C:\Users\gillissen\Desktop\Optical flow data\Frey1216\Frey_20161216_B1.mat')

resize=1;

x = LOG.currentdelay(ctrials{2});

tracethiscond = single(conddata(:,:,:,x==0));
clear conddata

if resize
    
    p=2; q=2;
    tmp = zeros(size(tracethiscond,1)/p,size(tracethiscond,2)/q,size(tracethiscond,3),size(tracethiscond,4));
    for i  = 1:size(tmp,4)
        for j = 1:size(tmp,3)

            M = tracethiscond(:,:,j,i);
            [m,n]=size(M); %M is the original matrix
            M=sum( reshape(M,p,[]) ,1 ); % when using [] (max one in reshape)  reshape calculates the size of that dimension to ensure numel(a) ==numel(b)
            M=reshape(M,m/p,[]).'; %Note transpose
            M=sum( reshape(M,q,[]) ,1);
            M=reshape(M,n/q,[]).'; %Note transpose
            tmp(:,:,j,i) = M;

        end
    end

    else 

    tmp = tracethiscond;
    

end

clear tracethiscond

I = imagesc(tmp(:,:,5,5)); % select brainmask before multiplying df/f by 100 otherwise colormap is out of range to see the brain...
BW = roipoly;
brainmask = BW;
mask = double(brainmask);

for i = 1:size(tmp,4) %loop over trials 
    
    V1 = tmp(:,:,:,i);
%     V1(~repmat(brainmask,[1,1,size(tmp,3)]))= nan;
    V1 = imgaussfilt(V1,2.5); %Gaussian filter
    tmp(:,:,:,i) = V1;            
end


%% Slow trent correction

trialidx = single(ctrials{2}(x==0));
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
    dFFav(i:i+99,:,:) = squeeze(nanmean(QQ,4));
end

clear QQ


% normalize between 0 and 1
% don't  normalize use binning and then  *100% df/f 
% mintmp = nanmin(trace);
% maxtmp = nanmax(trace);
% difftmp = maxtmp-mintmp;
% tmp = (trace-repmat(mintmp,size(trace,1),1))./repmat(difftmp,size(trace,1),1);



dFFav = dFFav(:,:,4:end); %frame 3 is bad


xvoor = double(dFFav*10);
% tmp = zeros(size(dFFav,1)/2,size(dFFav,2)/2,size(dFFav,3));

% for i = 1:size(dFFav,3)
%     tmp(:,:,i) = imresize(dFFav(:,:,i),0.5);
% end


for i = 1:size(dFFav,1)
    for j = 1:size(dFFav,2)
        
        dFFav(i,j,:) = medfilt1(dFFav(i,j,:),3);
    end
end

xna = double(dFFav*10);

dFFav = double(dFFav*10);




save('FreyC2','xvoor')
save('FreyC2filteredt','xna')

save('MaskfreyC22','mask') % OFAMM requires the mask to be named mask !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!@@@##


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


                
%% create movie




%Creating Video Writer Object
writerObj = VideoWriter('sequence.avi');
writerObj.FrameRate = 5;
% Using the 'Open` method to open the file
open(writerObj);

% Creating a figure.
% Each frame will be a figure data

axis tight
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

for k = 1:size(dFFav,3)
   imagesc(dFFav(:,:,k))
   % Frame includes image data
   frame = getframe;
   % Adding the frame to the video object using the 'writeVideo' method
   writeVideo(writerObj,frame);
end

% Closing the file and the object using the 'Close' method
close(writerObj);




























