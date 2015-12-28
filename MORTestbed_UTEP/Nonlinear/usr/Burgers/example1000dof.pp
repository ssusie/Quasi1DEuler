/* Notice in the model, modelAux, and modelFom fields, I have specified the key words
   fom, rom, tpwl, and gnat with a number appended to them.  This tells the PostProcess.m function
   which fom/rom/gnat/tpwl object to extract the requested data from.  Recall the last line of the workflow.m
   function where we pass {fom(1)}, {rom}, {gnat}, and {tpwl} into PostProcess.  This is trivial since all of
   these were use 1 x 1 arrays, but in general, we need to construct these cell array with many such objects,
   i.e. {fom(1), fom(2), ..., fom(N)}, etc. then in the .pp file, fom1 corresponds to the fom object
   placed in the first component of the cell array (fom(1)) in this case.  Identical explanation holds for
   rom, gnat, and tpwl.  See MEMS for more complex plotting functions. */
VAR:
   myFavLegend = {{'Location','SouthWest'}; {1, 'FOM'}; {2, 'PGROM'}; {3,'GNAT'}; {4,'TPWL'}}
   numForID4 = 3

FIGURE:
    id = 1
    axesLayout = [2,2]
    setfig = {'menubar','figure'}

    AXES:
      id = 1
      subplotNum = [1] 
      xData = {'speedup'} 
      yData = {'AvgRerr%'} 
      zData = [] 
      scaleData = [1,1,1]
      xlabel = {'Speedup','fontsize',14,'interpreter','latex'} 
      ylabel = {'Maximum Relative Error (2-norm)','fontsize',14,'interpreter','latex'} 
      zlabel = {} 
      title = {'Pareto Space','fontsize',16,'interpreter','latex'} 
      plotspec = {'ko','MarkerFaceColor','k','MarkerSize',8;...
                  'bs','MarkerFaceColor','b','MarkerSize',8;...
                  'gp','MarkerFaceColor','g','MarkerSize',8;...
                  'md','MarkerFaceColor','m','MarkerSize',8}
      model = {'fom1';'rom1';'gnat1';'tpwl1'} 
      modelAux = {'';'';'rom1';''} 
      modelFom = {'fom1'} 
      pType = {'plot'} 
      normType = 2 
      numPlotObj = 4 
      legend =  myFavLegend  
    /*  connWithCurve = {[1, 2, 3]; [1, 3, 4]; [1, 2, 4]} 
      connCurveSpec = {'k-','linewidth',2;...
                       'b-','linewidth',2;...
                       'm-','linewidth',2}*/

    AXES:
      id = 2
      subplotNum = 2
      xData = {'xdomain1'} 
      yData = {'clusCenter1@b=1';'clusCenter1@b=2';'pod1@n=1,b=2'} 
      zData = [] 
      scaleData = [1,1,1]
      xlabel = {'Position','fontsize',14,'interpreter','latex'} 
      ylabel = {'Conserved Quantity','fontsize',14,'interpreter','latex'} 
      zlabel = {} 
      title = {'1D Burger''s Equation','fontsize',16,'interpreter','latex'} 
      plotspec = {'bo-','linewidth',2,'Markersize',5,'MarkerFaceColor','b';...
                  'r-' ,'linewidth',2,'','','','';...
                  'm-' ,'linewidth',2,'','','',''} 
      model = {'rom1'} 
      modelAux = {'rom1'} 
      modelFom = {'fom1'} 
      pType = {'plot'} 
      normType = 2 
      numPlotObj = 3 
      legend = []  
      connWithCurve = [] 
      connCurveSpec = []

    AXES:
      id = 3
      subplotNum = [3] 
      xData = {'time'} 
      yData = {'sv1@x=40'} 
      zData = [] 
      scaleData = [1,1,1]
      xlabel = {'time (sec)','fontsize',14,'interpreter','latex'} 
      ylabel = {'Conserved Quantity','fontsize',14,'interpreter','latex'} 
      zlabel = {} 
      title = {'1D Burger''s Equation','fontsize',16,'interpreter','latex'} 
      plotspec = {'b-'}
      model = {'gnat1'} 
      modelAux = {'rom1'} 
      modelFom = {'fom1'} 
      pType = {'plot'} 
      normType = 2 
      numPlotObj = 1 
      legend = []  
      connWithCurve = [] 
      connCurveSpec = []  

    AXES:
      id = 4
      subplotNum = [4] 
      xData = {'xdomain1'} 
      yData = {'sv1@t=2.5';'sv1@t=20';'sv1@t=45';...
               'sv1@t=2.5';'sv1@t=20';'sv1@t=45';...
               'sv1@t=2.5';'sv1@t=20';'sv1@t=45';...
               'sv1@t=2.5';'sv1@t=20';'sv1@t=45'} 
      zData = [] 
      scaleData = [1,1,1]
      xlabel = {'Position','fontsize',14,'interpreter','latex'} 
      ylabel = {'Conserved Quantity','fontsize',14,'interpreter','latex'} 
      zlabel = {} 
      title = {'1D Burger''s Equation','fontsize',16,'interpreter','latex'} 
      plotspec = [repmat({'ko-','linewidth',2},numForID4,1); ...
                  repmat({'b-','linewidth',2},numForID4,1); ...
                  repmat({'g-','linewidth',2},numForID4,1);...
                  repmat({'m-','linewidth',2},numForID4,1)] 
      model = [repmat({'fom1'},numForID4,1);...
               repmat({'rom1'},numForID4,1); ...
               repmat({'gnat1'},numForID4,1);...
               repmat({'tpwl1'},numForID4,1)] 
      modelAux = [repmat({''},2*numForID4,1);repmat({'rom1'},numForID4,1);repmat({''},numForID4,1)] 
      modelFom = {'fom1'} 
      pType = {'plot'} 
      normType = 2 
      numPlotObj = 4*numForID4 
      legend = []  
      connWithCurve = [] 
      connCurveSpec = []


FIGURE:
    id = 2 
    axesLayout = [1,1] 

    AXES:
      id = 1
      subplotNum = [1] 
      xData = {'timestep'} 
      yData = {'prop.LocBasisHist'} 
      zData = [] 
      scaleData = [0.05,1,1]
      xlabel = {'time (sec)','fontsize',14,'interpreter','latex'} 
      ylabel = {'Number of Newton Iterations','fontsize',14,'interpreter','latex'} 
      zlabel = {} 
      title = {'1D Burger''s Equation: Newton Iteration Time History','fontsize',16,'interpreter','latex'} 
      plotspec = {'bo-';'go-'} 
      model = {'rom1';'gnat1'} 
      modelAux = {'';'rom1'} 
      modelFom = {'fom1'} 
      pType = {'plot'} 
      normType = 2 
      numPlotObj = 2 
      legend = []  
      connWithCurve = [] 
      connCurveSpec = []


FIGURE:
    id = 3 
    axesLayout = [1,2] 

    TABLE:
       id = 1 
       subtableNum = 1 
       normType = 1
       numRows = 2 
       numCols = 4 
       rowData = {} 
       colData = {'sv1@t=10','sv1@t=10','sv1@t=10','sv1@t=10'} 
       elemData = {} 
       title = 'Nonlinear Model Reduction Comparison: L2 Relative Error' 
       rowNames = {} 
       colNames = {'FOM','Petrov-Galerkin POD','GNAT','TPWL'} 
       model = {'fom1','rom1','gnat1','tpwl1'} 
       modelAux = {'','','rom1',''} 
       modelFom = {'fom1'}

    TABLE:
       id = 2 
       subtableNum = 2 
       normType = 1
       numRows = 1 
       numCols = 4 
       rowData = {} 
       colData = {'prop.nY','prop.nR','prop.nJ','prop.nI'} 
       elemData = {} 
       title = 'GNAT Parameters' 
       rowNames = {} 
       colNames = {'nY','nR','nJ','nI'} 
       model = {'gnat1'} 
       modelAux = {'rom1'} 
       modelFom = {'fom1'} 