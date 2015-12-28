% NEED APPROPRIAT MATLAB VERSION FOR THIS FUNCTION TO WORK
% WILL NOT WORK WITH MATLAB2010a
% WILL WORK WITH MATLAB2012b
% NOT SURE WHERE THE VERSION CUTOFF IS
%% Example for using writeVideo
% Generate an video object, like a file handle
writerObj = VideoWriter('burgers.avi'); 
% Set the FPS of the movie
set(writerObj, 'FrameRate',10); 
% Open the video object
open(writerObj);

% Prepare the movie
figure;
set(gca,'NextPlot','replaceChildren',...
    'xlim',[0,100],'ylim',[0,5.5]);
% Write the movie
for j = 1:10:1000
    plot(fom.sv(:,j)); %drawnow;
    frame = getframe; 
    % Directly write the frame to video
    writeVideo( writerObj, frame );
end 
% Close the video object
close(writerObj);