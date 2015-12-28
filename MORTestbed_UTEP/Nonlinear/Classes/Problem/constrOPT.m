

function  [] = constrOPT(obj)
            
            itnump1 = obj.cTimeIter + 1;
           
            t = obj.time.T(1) + obj.time.dt*obj.cTimeIter;
            [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1-1),obj.sv(:,itnump1-1),t);
            %Determine convergence criteria to use
            res0 = norm(R,2);
            [tol,tolIt] = determineConvergeCriterion(obj,res0);
            
            indexAdj=1;
            if obj.printLevel > 1 
                fprintf('---- Newton Step    # --------- ||R|| ------------ ||W''*R|| ----------- ||du|| ------------- tol -------------- tol_du ------ snap_coll ----------\n');
            end
            
            %constraints
            [rho,u,P,c,dcdU] = conservativeToPrimitive(obj,obj.sv(:, itnump1-1))
            [R,J,Cp,Cpp,dCp,dCpp] = ResJac(obj,U,~)
            [roeF,droeF] = roeFlux(obj,rho,u,P,c,e,dc)
            [R1,J1,Cpp,dCpp,dR1dp] = fullImplicitInletBC(obj,rho,u,P,c,dc)
            
            for i_N = 1:obj.newt.maxIter
                %Solve for the search direction using Newton or Gauss-Newton
                if lower(obj.Rtype(1)) == 'g'
                    newR = obj.phi'*R;
                    du = -((obj.phi'*J*obj.phi)\newR);
                    conv = norm(newR,2);
                    %conv(i_N) = norm(newR,2);
                elseif strcmpi(obj.Rtype,'Petrov-Galerkin')
                    newJ = J*obj.phi;
                    if obj.augment
                        if i_N == 1
                            S_k = zeros(obj.nY);
                            du = -(newJ\R);
                        else
                            [Ro,Jo] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1)-p,obj.sv(:,itnump1-1),t);
                            S_k = obj.GN_HessianMod(S_k,du,newJ,R,Jo*obj.phi,Ro);
                            du = -((newJ'*newJ + S_k)\(newJ'*R));
                        end
                    else
                        du = -(newJ\R);
                    end
                    
                    conv = norm(newJ'*R,2);
                elseif strcmpi(obj.Rtype,'Petrov-Galerkin-Regularized')
%                     mu    = 0.1;
%                     sigma = 1;
                    
                    newJ = [ J*obj.phi; obj.mu*eye(obj.nY)] ;
                    newR = [ obj.sigma*R ; zeros(obj.nY,1) ];
                    du = -(newJ\newR);
                    conv = norm(newJ'*newR,2);
                    %conv(i_N) = norm(newJ'*R,2);
                end
                
                if i_N == 1, [tol,tolIt] = determineConvergeCriterion(obj,conv); end;
                
                %If nan or complex encountered -> kill simulation
                if sum(isnan(du)) || ~isreal(du)
                    obj.killflag=true;
                    return;
                end
               
                %Linesearch if necessary
                if obj.newt.linesrch.status
                    %alpha = linesearch(@(alpha) obj.TimeScheme.TimeIntNLFuncDir(obj.prob,obj.sv(:,itnump1-indexAdj),obj.sv(:,itnump1-1),t,alpha,du,obj),obj.newt.linesrch.prop);                
                    alpha = linesrchBackNewton(@(alpha) obj.TimeScheme.TimeIntNLFuncDir(obj.prob,obj.sv(:,itnump1-indexAdj),obj.sv(:,itnump1-1),t,alpha,du,obj));
                    du = alpha*du;
                end
                %Reconstructed search direction
                p = obj.phi*du;
                
                obj.sv(:,itnump1)=obj.sv(:,itnump1-indexAdj)+p;
                
                if obj.printLevel > 2
                    fprintf('---- Newton Step %4i ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ------- %1i --------------\n',i_N,norm(R,2),conv,norm(du,2),tol,tolIt,obj.saveNL);
                end
            
                %WHY DO WE SKIP THE FIRST ONE?
                if i_N > 1 || obj.newt.maxIter == 1
                    %Write nonlinear snapshots
                    obj.writeNonlinearSnapshot(R,J,p,1);
                end
                
                indexAdj=0;
                if checkConverge(obj,conv,du,tol,tolIt)
                    break;
                end
                
                %Compute residual and jacobian with updated vector for next iteration
                [R,J] = obj.TimeScheme.TimeIntNLFunc(obj.prob,obj.sv(:,itnump1),obj.sv(:,itnump1-1),t);
            end
                
            %Write the last nonlinear snapshot.  We note that we are using
            %J^(k+1)*0 if using Snapshot 1.5 or 2 (we don't have a
            %search direction p^(k+1) because we converged!).  It shouldn't
            %matter much because ||p|| should be small since we converged,
            %i.e. the next step shouldn't take us away from the solution!
            %obj.writeNonlinearSnapshot(R,J,zeros(size(p)),1);
            
            %Store the number of newton iterations at this time step
            obj.newt.iter(obj.cTimeIter) = i_N;
            
            if obj.printLevel == 1.5
                fprintf('---- Newton Step %4i ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ----- %10.7e ------- %1i --------------\n',i_N,norm(R,2),conv,norm(p,2),tol,tolIt,obj.saveNL);
            end
            
            %If the maxiumum newton iterations are reached, warn the user
            if ~obj.newt.quiet && (i_N==obj.newt.maxIter)
                disp(['*** Warning; Newton solver reached max number of iterations before convergence : res = ',...
                    num2str(norm(R,2)/tol),'  eps = ',num2str(obj.newt.eps), '***']);
            end
        end %Done

