function s = shower_tile_plot(inmat)
% A function to plot matrix data with neat lines between everything like
% pcolor...but with all the weird intermediary steps necessary to make sure
% everything's the right way up, and that we don't lose a row/column of
% data. For plotting PLS results, reshape and threshold them before getting
% here. I haven't included syntax for axis ticks, labelling, or titles, so
% add those in the plot script 
%%
nrows = size(inmat,1);
ncols = size(inmat,2);
inmat = flipud(inmat);
inmat(nrows+1,1:ncols) = nan;
inmat(:,ncols+1) = nan;
s = pcolor(inmat);
s.LineWidth = 1.5;
colorbar
