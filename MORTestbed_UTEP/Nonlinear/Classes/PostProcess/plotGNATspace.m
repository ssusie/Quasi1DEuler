function  [axOUT] = plotGNATspace(gnatobj,fomobj,romobj,probobj,nType,axIN)
%This function plots the information in gnatobj in GNAT space.
%--------------------------------------------------------------------------
%Inputs:
%-------
%gnatobj - vector of GNAT objects
%fomobj  - FOM object
%romobj  - vector of (or single) ROM objects
%probobj - Problem object
%nType   - positive double (or inf) indicating the norm type to use in
%          computing the error
%axIN    - axes handle to plot the 
%Set up axes
if isempty(axIN)
    figure('Position',[325   272   675   526]);
    ax = axes;
    title('GNAT Space','fontsize',16,'interpreter','latex');
    xlabel('$n_R$','fontsize',14,'interpreter','latex');
    ylabel('$n_J$','fontsize',14,'interpreter','latex');
    zlabel('$n_I$','fontsize',14,'interpreter','latex');
else
    ax = axIN;
end
axOUT = ax;
set(ax,'nextplot','add');

%Define the number of color bins to use and the "color increment"
bins = 10; %11 bin endpoints
inc = 10;
increment = 0;
%Set the colormap for the GNAT points (in "GNAT space", a subset of
%Z^3)
colormap(jet(bins*inc+1));
cmp = colormap;

%Setup array to hold error values for the various values of
%nR, nJ, nI
err = zeros(length(gnatobj),1);
%Loop over all gnat objects
for i = 1:length(gnatobj)
    increment = increment + 1;
    if gnatobj(i).nI < gnatobj(i).nR || gnatobj(i).nI < gnatobj(i).nJ
        %If this is an illegal combination of nR, nJ,
        %and nI (i.e. a point of Z^3 that is not
        %contained in "GNAT space", continue to the
        %next point)
        continue;
    else
        %If this is a point in "GNAT space", run the
        %GNAT model for this combination of nR, nJ, nI,
        %store the error
        if length(romobj) == 1
            svGNAT = gnatobj(i).reconstructFullState(probobj,romobj);
        elseif length(romobj) == length(gnatobj)
            svGNAT = gnatobj(i).reconstructFullState(probobj,romobj(i));
        else
            error('ROM object must be the same length as the GNAT object or it must be of length 1');
        end
        temp = computeRelativeError(svGNAT,fomobj.sv,nType);
        err(i,1) = max(temp);
        %Determine the appropriate color for this point
        %in "GNAT space" based on the error
        cstr = determineColor(err(i,1),cmp,inc);
        %Plot the point with the appropriate color
        %(indicating the error) in GNAT space
        plot3(gnatobj(i).nR,gnatobj(i).nJ,gnatobj(i).nI,...
            'o','MarkerSize',12,'MarkerFaceColor',cstr,'MarkerEdgeColor',cstr);
    end
end

function  [cstr] = determineColor(err,cmp,inc)
%This function determines the appropriate color of the the error specified
%in the err variable based on the colormap array cmp and the color
%increment inc.
%--------------------------------------------------------------------------
%Inputs:
%-------
%err  - scalar specifying the level of error to determine a color for
%cmp  - 3 column matrix where each row specifies a color in the colormap
%inc  - a scalar specifying the intensity difference in adjacent color bins
%
%Outputs:
%--------
%cstr - a 1 x 3 vector specifying the RGB values corresponding to err
%--------------------------------------------------------------------------

%Adjust error values into a percentage
err = err*100;

if err > 40 || isnan(err)
    %If the error is greater than 40% or is not a number, use the last row
    %of cmp
    cstr = cmp(end,:);
elseif err > 30
    %If the error is greater than 30% and less than 40%, use appropriate
    %row of cmp
    cstr = cmp(1+inc*9,:);
elseif err > 20
    %If the error is greater than 20% and less than 30%, use appropriate
    %row of cmp
    cstr = cmp(1+inc*8,:);
elseif err > 10
    %If the error is greater than 10% and less than 20%, use appropriate
    %row of cmp
    cstr = cmp(1+inc*7,:);
elseif err > 5
    %If the error is greater than 5% and less than 10%, use appropriate
    %row of cmp
    cstr = cmp(1+inc*6,:);
elseif err > 4
    %If the error is greater than 4% and less than 5%, use appropriate
    %row of cmp
    cstr = cmp(1+inc*5,:);
elseif err > 3
    %If the error is greater than 3% and less than 4%, use appropriate
    %row of cmp
    cstr = cmp(1+inc*4,:);
elseif err > 2
    %If the error is greater than 2% and less than 3%, use appropriate
    %row of cmp
    cstr = cmp(1+inc*3,:);
elseif err > 1
    %If the error is greater than 1% and less than 2%, use appropriate
    %row of cmp
    cstr = cmp(1+inc*2,:);
elseif err > 0.5
    %If the error is greater than 0.5% and less than 1%, use appropriate
    %row of cmp
    cstr = cmp(1+inc,:);
else
    %If the error is less than 0.5%, use appropriate row of cmp
    cstr = cmp(1,:);
end