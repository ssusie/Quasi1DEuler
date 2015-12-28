classdef OPTprob < handle
    
    properties
        objective;
        
        ce;
        ci;
        
        Ae; be;
        Ai; bi;
        
        xl; xu;
    end
    
    methods
        function  [obj] = OPTprob(objectiveIn,AeIn,beIn,AiIn,biIn,ceIn,ciIn,xlIn,xuIn)
            
            obj.objective = objectiveIn;
            obj.ce = ceIn;
            obj.ci = ciIn;
            obj.Ae = AeIn;
            obj.be = beIn;
            obj.Ai = AiIn;
            obj.bi = biIn;
            obj.xl = xlIn;
            obj.xu = xuIn;
            
        end
        
        function  [L,dL] = Lagrangian()
            
        end
        
        function  [Q,dQ] = QuadraticPenalty()
            
        end
    end
end