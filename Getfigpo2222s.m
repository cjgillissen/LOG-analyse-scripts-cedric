%% Get figure position
% construct logical map with the location of the figure position

% load('\\VC2NIN\wbimaging\pRF Results\Marsellus\4Mapping\AllenBorders')
load('\\VC2NIN\wbimaging\MainAnalysis\Marsellus\brainareamodel.mat')
load('\\\\VC2NIN\wbimaging\pRF Results\Marsellus\4Mapping\pRFmaps')

fgy = 7;
fgsz = 50;
r = [0:1:359];
area = 11;
z = 0;
edge = 2.5; %+/- this much around edge
x = fgx+cosd(r).*(fgsz/2);
y = fgy+sind(r).*(fgsz/2);
figmap = zeros(800,800);
AZIs = imgaussfilt(AZIi,3);
ELEs = imgaussfilt(ELEi,3);

 for side = [-1 1]
    z = z+1;
    IN = Model.Regions{area}';
    IN = IN(1:800,300:800);
%     if side == -1
%         IN(y>Model.Bregma(2)) = 0;
%     else
%         IN(y<Model.Bregma(2)) = 0;
%     end
    INv = find(IN); %Indices describing position of V1 in the AZI/ELE maps 9bilateral!)
    azi = AZIs(INv);
    ele = ELEs(INv);
%     azirange = (0:1:70).*side;
%     elerange = -40:1:50;
%     [a,e] = meshgrid(azirange,elerange);

    %Where was figure?
    fgx = 40.*side;
    x = fgx+cosd(r).*(fgsz/2);
    y = fgy+sind(r).*(fgsz/2);

    %EXTRACT data based on distance from fig etc
    %Get the distance of each V1 point to teh figure center in spatial
    %co-ordinates
    ar = azi-fgx;
    er = ele-fgy;
    rdist = hypot(ar,er);
    v1idx = sq2brain(INv);
    infig  = v1idx(rdist<(fgsz/2)); 
    figmap(infig) = 1;
    figmap(isnan(cat(2,nan(800,299),FSMit)))=0;
        
 end
