classdef  genfunc < handle
   
    properties
        func     = [];
        funcJac  = [];
        funcHess = [];
        
        JacStatus     = 1; %0 = no jacobian, 1 = exact jacobian, 2 = FD jacobian
        HessianStatus = 1; %0 = no hessian,  1 = exact hessian,  2 = FD hessian
        
        eps = 1e-6;
        M;
        N;
    end
    
    methods
        function  [obj] = genfunc(funcIn,jacIn,hessIn,Min,Nin,JacStat,HessStat)
            
            obj.func = funcIn;
            obj.funcJac = jacIn;
            obj.funcHess = hessIn;
            
            obj.M = Min; %number of components of function
            obj.N = Nin; %number of unknowns
            obj.JacStatus = JacStat;
            obj.HessianStatus = HessStat;
            
            
        end
        
        function  [val] = funcVal(obj,x)
            
            val=obj.func(x);
            
        end
        
        function  [jac] = funcJacobian(obj,x)
            
            if obj.JacStatus == 1
                jac=obj.funcJac(x);
            elseif obj.JacStatus == 2
                jac=obj.funcJacFD(x);
            end
                        
        end
        
        function  [hess] = funcHessian(obj,x)
            
            if obj.HessianStatus == 1
                hess = obj.funcHess(x);
            elseif obj.HessianStatus == 2
                hess = obj.funcHessFD(x);
            end
            
        end
        
        function  [jac] = funcJacFD(obj,x)
            
            jac=zeros(obj.M,obj.N);
            e=zeros(obj.N,1);
            for j = 1:obj.N
                e(j)=1;
                jac(:,j) = (0.5/obj.eps)*(obj.funcVal(x+obj.eps*e)-obj.funcVal(x-obj.eps*e));
                e(j)=0;
            end
            
        end
        
        function  [hess] = funcHessFD(obj,x)
            
            if obj.M > 1
                error('Hessian not implemented for vector-valued functions');
            end
            
            hess=zeros(obj.N,obj.N);
            e=zeros(obj.N,1);
            for j = 1:obj.N
                e(j)=1;
                hess(:,j) = (0.5/obj.eps)*(obj.funcJacobian(x+obj.eps*e)-obj.funcJacobian(x-obj.eps*e))';
                e(j)=0;
            end
            
        end
    end
end