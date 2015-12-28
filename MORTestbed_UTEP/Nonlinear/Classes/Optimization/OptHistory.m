classdef OptHistory < handle
   
    properties (GetAccess=public,SetAccess=private)
       objective = [];
       constraintviolation = [];
       iterates = [];
   end
   
   methods
       function [obj] = OptHistory()
           return;
       end
       
       function [] = incConstrViolate(obj,val)
           obj.constraintviolation = [obj.constraintviolation;val];
       end
       
       function [] = incObj(obj,val)
           obj.objective = [obj.objective;val];
       end
       
       function [] = incIter(obj,val)
           obj.iterates = [obj.iterates,val];
       end
       
       function [] = reshapeIter(obj)
           [a,b] = size(obj.iterates);
           
           if (b == 1)
               niter = length(obj.objective);
               nvar  = a/niter;
               obj.iterates = reshape(obj.iterates,nvar,niter);
           else
               obj.iterates = reshape(obj.iterates,a*b,1);
           end
       end
   end
   
   methods (Static)
       function [globalOpt,globalOptIndex] = globalOpt(history,constrviolationTol)
           globalOptIndex=nan;
           globalOpt = inf;
           for i = 1:length(history)
               if (history(i).constraintviolation(end) < constrviolationTol) && (history(i).objective(end) < globalOpt)
                   globalOpt = history(i).objective(end);
                   globalOptIndex=i;
               end
           end
       end
   end
end