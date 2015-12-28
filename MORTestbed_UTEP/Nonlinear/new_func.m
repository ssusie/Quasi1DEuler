function  [] = new_func(obj)
            
            itnump1 = obj.cTimeIter + 1;
                              
            %%%%%%% defining constriants %%%%%%
            U=reshape(obj.sv(:,obj.cTimeIter),3,obj.nVol); 
            rho(1)=U(1,1); u(1)=U(2,1); P(1)=U(3,1);
            rho(2)=U(1,end); u(2)=U(2,end); P(2)=U(3,end);
            e(1)= P(1)/(obj.gamma-1)+rho(1)*u(1)^2/2;
            e(2)= P(end)/(obj.gamma-1)+rho(end)*u(end)^2/2;
           
            Um1=reshape(obj.sv(:,obj.cTimeIter-1),3,obj.nVol);
            rhom1(1)=Um1(1,1); um1(1)=Um1(2,1); pm1(1)=Um1(3,1);
            rhom1(2)=Um1(1,end); um1(2)=Um1(2,end); pm1(2)=Um1(3,end);
            em1(1) = Pm1(1)/(obj.gamma-1)+rhom1(1)*um1(1)^2/2;
            em1(2) = Pm1(end)/(obj.gamma-1)+rhom1(end)*um1(end)^2/2;            
            
            %Flux = [rho.*u; rho.*u.*u+P; (e+P).*u];
            flx_a(1,1)=rho(1)*u(1); flx_a(2,1)=rho(1)*u(1)*u(1)+P(1); 
            flx_a(3,1)= (e(1)+P(1))*u(1);
            flx_b(1,1)=rho(2)*u(2); flx_b(2,1)=rho(2)*u(2)*u(2)+P(2); 
            flx_b(3,1)= (e(2)+P(2))*u(2);
                        
            flx_am1(1,1)=rhom1(1)*um1(1); flx_am1(2,1)=rhom1(1)*um1(1)*um1(1)+Pm1(1); 
            flx_am1(3,1)=(em1(1)+Pm1(1))*um1(1);
            flx_bm1(1,1)=rhom1(2)*um1(2); flx_bm1(2,1)=rhom1(2)*um1(2)*um1(2)+Pm1(2); 
            flx_bm1(3,1)=(em1(2)+Pm1(2))*um1(2);
                        
            vol=obj.mesh.node(2:end,1)-obj.mesh.node(1:end-1,1); %column vector
            U(2,1)=rho(1)*u(1); %first column is primitive variables
            U(3,1)=e(1);
            U(2,end)=rho(end)*u(end); %last column is primitive variables
            U(3,end)=e(end);
            
            cnstr(1)=(U(1,1:end-1)-Um1(1,1:end))*vol-obj.time.dt/2*(flx_a(1)+flx_am1(1)-flx_b(1)-flx_bm1(1));
            cnstr(2)=(U(2,1:end-1)-Um1(2,1:end))*vol-obj.time.dt/2*(flx_a(2)+flx_am1(2)-flx_b(2)-flx_bm1(2));
            cnstr(3)=(U(3,1:end-1)-Um1(3,1:end))*vol-obj.time.dt/2*(flx_a(3)+flx_am1(3)-flx_b(3)-flx_bm1(3));
                      
            
            %Determine the residual and jacobian based on initial guess
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1-1),obj.sv(:,itnump1-1),t);
            % determine Hessian of constraint and residuale HR, HC
            
            
            
            %Determine convergence criteria to use
            res0 = norm(R,2);
            [tol,tolIt] = determineConvergeCriterion(obj,res0);
            
            indexAdj=1;
            if obj.printLevel > 1 
                fprintf('---- Newton Step    # --------- ||R|| ------------ ||W''*R|| ----------- ||du|| ------------- tol -------------- tol_du ------ snap_coll ----------\n');
            end
            
            %initial guesses for lambda and solution
            num_cnstr=3;
            l=3*ones(num_cnstr,1); %size of constraint 
            x=obj.sv(:,obj,cTimeIter);
            
                     
            for i_N = 1:obj.newt.maxIter
               %obj.nY, obj.nY(basisNum)
                H=[HR+HC'*l, dcnstr'; dcnstr, zeros(num_cnst)];
                h=[J+dcnstr'*l;dcnstr];
                sol=H\h;
                dx=sol(1,end-num_cnstr); dl=sol(end-num_cnstr:end);
                
                %If nan or complex encountered -> kill simulation
                if sum(isnan(dx)) || ~isreal(dx)
                    obj.killflag=true;
                    return;
                end
                
                x=x+dx; l=l+dl;
                p=obj.phi*dx; 
                obj.sv(:,itnump1)=obj.sv(:,itnump1-1)+p;
                
                %check criterion
                
               
                %update HR, HC,dcnstr,J
                
                               

            end
                
            
%             obj.newt.iter(obj.cTimeIter) = i_N;
%             
%             if obj.printLevel == 1.5
%                 fprintf('---- Newton Step %4i ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ------- %1i --------------\n',i_N,norm(R,2),conv,norm(p,2),tol,tolIt,obj.saveNL);
%             end
%             
%             %If the maxiumum newton iterations are reached, warn the user
%             if ~obj.newt.quiet && (i_N==obj.newt.maxIter)
%                 disp(['*** Warning; Newton solver reached max number of iterations before convergence : res = ',...
%                     num2str(norm(R,2)/tol),'  eps = ',num2str(obj.newt.eps), '***']);
%             end
end
