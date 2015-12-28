function  [] = opt_func(obj)
            
            itnump1 = obj.cTimeIter + 1;
                              
            %%%%%%% defining constriants %%%%%%
            U=reshape(obj.sv(:,obj.cTimeIter),3,obj.nVol); 
            rho(1)=U(1,1); u(1)=U(2,1); P(1)=U(3,1);
            rho(2)=U(1,end); u(2)=U(2,end); P(2)=U(3,end);
            e(1)= P(1)/(obj.gamma-1)+rho(1)*u(1)^2/2;
            e(2)= P(end)/(obj.gamma-1)+rho(end)*u(end)^2/2;
            
%             c   = sqrt(obj.gamma*P(1)/rho(1)); %Speed of sound
           
            %Flux = [rho.*u; rho.*u.*u+P; (e+P).*u];
            flx_a(1,1)=rho(1)*u(1); flx_a(2,1)=rho(1)*u(1)*u(1)+P(1); 
            flx_a(3,1)= (e(1)+P(1))*u(1);
            flx_b(1,1)=rho(2)*u(2); flx_b(2,1)=rho(2)*u(2)*u(2)+P(2); 
            flx_b(3,1)= (e(2)+P(2))*u(2);
            
            Phi1=obj.phi(1:3:end,:); %basis for the first conserved quantity rho
            Phi2=obj.phi(2:3:end,:); %basis for the second conserved quantity rho*u
            Phi3=obj.phi(3:3:end,:); %basis for the thirs conserved quantity e
                        
            vol=(obj.mesh.node(2:end,1)-obj.mesh.node(1:end-1,1))'; %row vector
            U(2,1)=rho(1)*u(1); %first column is primitive variables
            U(3,1)=e(1);
            U(2,end)=rho(end)*u(end); %last column is primitive variables
            U(3,end)=e(end);
            
            w0 = obj.phi'* svobj.sv(:,obj.cTimeIter); %solution at previous time step
            w_guess = obj.phi'* svobj.sv(:,obj.cTimeIter); %initial guess for generalized coordinates for Newton iteration
            l_guess = [0;0;0]; % guess for lambda in optimization problem for Newton interation
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            
            if obj.printLevel > 1 
                fprintf('---- Newton Step    # --------- ||R|| ------------ ||du|| ------------- \n');
            end
            for k=1:obj.newt.maxIter
                
                sv_guess= obj.phi*w_guess;
                
                Ug=reshape(sv_guess,3,obj.nVol); 
                rhog(1)=Ug(1,1); ug(1)=Ug(2,1); Pg(1)=Ug(3,1);
                rhog(2)=Ug(1,end); ug(2)=Ug(2,end); Pg(2)=Ug(3,end);
                eg(1)= Pg(1)/(obj.gamma-1)+rhog(1)*ug(1)^2/2;
                eg(2)= Pg(end)/(obj.gamma-1)+rhog(end)*ug(end)^2/2;
                cg(1) = sqrt(obj.gamma*Pg(1)/rhog(1)); %Speed of sound at boundary a
                cg(1) = sqrt(obj.gamma*Pg(2)/rhog(2)); %Speed sound at boundary b
                flxGuess_a(1,1) = rhog(1)*ug(1); flxGuess_a(2,1)=rhog(1)*ug(1)*ug(1)+Pg(1); 
                flxGuess_a(3,1) = (eg(1)+Pg(1))*ug(1);
                flxGuess_b(1,1) = rhog(2)*ug(2); flxGuess_b(2,1)=rhog(2)*ug(2)*ug(2)+Pg(2); 
                flxGuess_b(3,1) = (eg(2)+Pg(2))*ug(2);                       
            
                constr(1)=vol*Phi1*(w_guess-w0)-obj.time.dt/2*(flxGuess_a(1)+flx_a(1)-flxGuess_b(1)-flx_b(1));
                constr(2)=vol*Phi2*(w_guess-w0)-obj.time.dt/2*(flxGuess_a(2)+flx_a(2)-flxGuess_b(2)-flx_b(2));
                constr(3)=vol*Phi1*(w_guess-w0)-obj.time.dt/2*(flxGuess_a(3)+flx_a(3)-flxGuess_b(3)-flx_b(3));
                
                [Sg_a,Sgi_a,Cg_a,Cgi_a,Lam_a]=eulerJac(obj,rhog(1), ug(1),cg(1);
                [Sg_b,Sgi_b,Cg_b,Cgi_b,Lam_b]=eulerJac(obj,rhog(2), ug(2),cg(2);
                dF_a=Sgi_a*Cgi_a*Lam_a*Cg_a*Sg_a;
                dF_b=Sgi_b*Cgi_b*Lam_b*Cg_b*Sg_b;
                
                Psi_a= [Phi1(1,:);Phi2(1,:); Phi3(1,:)];
                Psi_b= [Phi1(end,:);Phi2(end,:); Phi(end,:)];
                Jconstr=[ vol*Phi1+dF_b(1,:)*Psi_b-dF_a(1,:)*Psi_a;
                          vol*Phi2+dF_b(2,:)*Psi_b-dF_a(2,:)*Psi_a;
                          vol*Phi3+dF_b(3,:)*Psi_b-dF_a(3,:)*Psi_a; ];
                t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;      
                [Res,Jres]=obj.TimeScheme.TimeIntNLFunc(obj.prob,sv_guess,sv_guess,t);
                Hres=Jres'*Jres;
%                 Hconstr=Jconstr'*Jconstr;
                H = [Hres, Jconstr'; Jconstr, zeros(3)];
                h = [Jres*Res+Jconstr'*l_guess;constr];
                delta=H\h;
                w_guess=w_guess+delta(1:end-3);
                l_guess=l_guess+delta(end-2:end);
                
                if k==1
                     res0 = norm(Res,2);
                     [tol,tolIt] = determineConvergeCriterion(obj,res0);
                end
%                 
%                 if checkConverge(obj,conv,du,tol,tolIt)
%                     break;
%                 end
%                     
                if norm(delta,2)<10^(-5)
                    break;
                    
                end
                if obj.printLevel > 2
                    fprintf('---- Newton Step %4i ----- %10.7e ----- %10.7e -----\n',...
                        k,norm(Res,2),norm(delta,2));
                end
            end
                      
            
            
end
