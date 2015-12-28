classdef Contrived < handle
    
    properties (SetAccess = public, GetAccess = public)
        config;
        staticFlag = false;
        ic;
        B;
        G;
        C;
        D;
        
        p;
        param;
    end
    
    methods
        function  [obj] = Contrived(cfgobj,oldobj)
            %This is the constructor for the 2D contrived example.  If
            %there is 1 input, the constructor will create the class
            %instance from the variables defined in this functino.  If
            %there are 2 inputs, the first input will be unused and the
            %output will be a deep copy of the oldobj.
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %cfgobj - CONFIG object
            %oldobj - OneDBurgers object whose properties will be copied
            %         into the new object.
            %
            %Outputs:
            %--------
            %obj - OneDBurgers object
            %--------------------------------------------------------------
            
            if nargin == 0
                %Keep everything at defaults if no input arguments
                return;
            end
            
            % Construct a new object based on a deep copy of an old object
            % of this class by copying properties over.
            if nargin == 2 && isa(oldobj,'Contrived')
                props = properties(oldobj);
                for i = 1:length(props)
                    % Use Dynamic Expressions to copy the required property.
                    % For more info on usage of Dynamic Expressions, refer to
                    % the section "Creating Field Names Dynamically" in:
                    % web([docroot '/techdoc/matlab_prog/br04bw6-38.html#br1v5a9-1'])
                    obj.(props{i}) = oldobj.(props{i});
                end
                return;
            end

            %Set configuration
            obj.config = cfgobj;
            
            obj.param = obj.config.param;
            
            %Determine ic term vector from the ic function and the
            %coordinate vector
            if isa(cfgobj.icFunc,'function_handle')
                obj.ic = cfgobj.icFunc();
            else
                obj.ic = feval(cfgobj.icFunc);
            end
        end %Done
        
        function  [R,J] = ResJac(obj,U,t)
            %This function calcuates the residual and jacobian of Burger's Equation
            %--------------------------------------------------------------
            %Inputs:
            %-------
            %obj - OneDBurgers object
            %U   - N x 1 vector containing the "solution" (or guess or
            %      approximation) at the current time
            %t   - scalar indicating the time
            %
            %Outputs:
            %--------
            %R   - N x 1 vector containing the nonlinear residual of the solution
            %      U at the current time t
            %J   - N x N vector containing the nonlinear jacobian of the solution
            %      U at the current time t
            %--------------------------------------------------------------
            
            switch obj.param{1}
                case 'local-rom'
                    R(1,1) = 1/(U(1)^2+U(2)^2);
                    R(2,1) = -sin(U(1))/(U(1)^2+U(2)^2);
                    
                    J(1,1) = -2*U(1)/(U(1)^2+U(2)^2);
                    J(1,2) = -2*U(2)/(U(1)^2+U(2)^2);
                    
                    J(2,1) = -cos(U(1))/(U(1)^2+U(2)^2) + 2*U(1)*sin(U(1))/(U(1)^2+U(2)^2);
                    J(2,2) = 2*U(2)*sin(U(1))/(U(1)^2+U(2)^2);
                case 'param-ic'
                    R(1,1)=-U(1)*(U(1).^2+U(2).^2);
                    R(2,1)=-U(2)*(U(1).^2+U(2).^2);
                    
                    J(1,1)=-2*U(1).^2-(U(1).^2+U(2).^2);
                    J(1,2)=-U(1)*U(2);
                    J(2,1)=-U(1)*U(2);
                    J(2,2)=-2*U(2).^2-(U(1).^2+U(2).^2);

%                     R(1,1) = 1 + U(1).^2.*U(2);
%                     R(2,1) = U(2);
%                     
%                     J(1,1) = 2*U(1)*U(2);
%                     J(1,2) = U(1).^2;
%                     
%                     J(2,1) = 0;
%                     J(2,2) = 1;
            end
        end %Done
        
        function  [] = setPrevSV(obj,prevsv)
            return;
        end
        
        function  [] = setProperty(obj,prop,val)
            obj.(prop)=val;
        end
    end
end