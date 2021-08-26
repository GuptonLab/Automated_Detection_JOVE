%Read in the neurite mask file and skeletonize it
mask_name = 'MaskFiles/VAMP2pHluorin_488_wt_4_mask_file_neur.tif';
csv_file_name = 'DataFiles/VAMP2pHluorin_488_wt_4_fluorescent_traces.csv';


I = imread(mask_name);
I = logical(I);
B = bwskel(I);

%read in the x,y,t file
pos_file = readtable(csv_file_name); %read in the xyt file
x = pos_file.x_pos;
y = pos_file.y_pos;

%remove points not in ROI

%
x2 = 1:size(B,2);
y2 = 1:size(B,1);
[X,Y] = meshgrid(x2,y2);

for i = 1:size(x,1)
    x_t = x(i)
    y_t = y(i)
    X_t=X-x_t;Y_t=Y-y_t;
    D=sqrt(X_t.^2+Y_t.^2);
    D(~B)=inf;%set all 0 pixel to infinity distance so they won't be selected
    [y_ans(i),x_ans(i)]=find(D==min(D(:)),1);%find a/the global minimum
end

radi = repmat(10,size(x_ans,2),1); %constant: radius for circle
temp_image = insertShape(uint8(B),'Circle',[transpose(x_ans) transpose(y_ans) radi], 'LineWidth',4, 'Color','red');

matr=transpose([y_ans;x_ans]);
writematrix(matr,"DataFiles/neurite_exocytic_events.csv");
writematrix(B,"MaskFiles/VAMP2pHluorin_488_wt_4_mask_file_neur.csv");
