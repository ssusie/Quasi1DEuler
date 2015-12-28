classdef nlHF
    
    properties
        config;
        ic;
        B;
        C;
        D;
        G;
        dx;
        HFparam;
    end
    
    methods
        function  [obj] = nlHF(cfgobj)
            obj.config = cfgobj;
            
            [Coord,obj.HFparam,obj.C] = feval(obj.config.altFile,obj.config.nNodes);

        end
        
        function [Minv_hat,M_bar,F_hat,G_bar,G_hat, T_dbc_Prev, T_dbc] = evalFullGlobalMAT(obj,t,dt)
            %This function evaluates the full global matrices (mass, stiffness, etc)
            %for the models that require them (i.e. FOM, GROM, PGROM, TPWL)
            %--------------------------------------------------------------------------
            %Inputs:
            %-------
            %t         -
            %dt        -
            %param     -
            %
            %Outputs:
            %--------
            %GlobalMAT -
            %DBCvals   -
            %--------------------------------------------------------------------------
            
            %Extract the appropriate parameters
            LM = obj.HFparam.LM;
            
            nip = obj.HFparam.quad{1,1};
            Z = obj.HFparam.quad{2,1};
            W = obj.HFparam.quad{3,1};
            
            BCCode = obj.HFparam.BCs{1,1};
            DBCfunc = obj.HFparam.BCs{4,1};
            
            kfunc = obj.HFparam.MATprop{1,1};
            rcfunc = obj.HFparam.MATprop{2,1};
            betafunc = obj.HFparam.MATprop{3,1};
            sefunc = obj.HFparam.MATprop{4,1};
            Tinffunc = obj.HFparam.MATprop{5,1};
            qfunc = obj.HFparam.MATprop{6,1};
            
            nel = obj.HFparam.mesh{1,1};
            nen = obj.HFparam.mesh{2,1};
            nnp = obj.HFparam.mesh{3,1};
            xhat = obj.HFparam.mesh{4,1};
            Belem = obj.HFparam.mesh{5,1};
            
            L = obj.HFparam.domain{1,1};
            rfunc = obj.HFparam.domain{2,1};
            
            numel = nel(1)*nel(2);
            
            M = zeros(nnp,nnp);
            K = zeros(nnp,nnp);
            H = zeros(nnp,nnp);
            F = zeros(nnp,1);
            
            for iel = 1:numel
                Me = zeros(nen,nen);
                Ke = zeros(nen,nen);
                Fe = zeros(nen,1);
                
                xe = xhat([LM(1,iel);LM(2,iel);LM(3,iel);LM(4,iel)],1);
                ye = xhat([LM(1,iel);LM(2,iel);LM(3,iel);LM(4,iel)],2);
                
                for i = 1:nip
                    
                    for j = 1:nip
                        [Ne,~,J,Bmat] = CalcInterpFunctions(Z(i),Z(j),xe,ye,nen);
                        
                        x = Ne*xe;
                        y = Ne*ye;
                        
                        k = kfunc(x,y);
                        r = rfunc(x,y,t); %body force at time = 0
                        rc = rcfunc(x,y);
                        
                        Me = Me + Ne'*rc*Ne*det(J)*W(i)*W(j);
                        Ke = Ke + Bmat'*k*Bmat*det(J)*W(i)*W(j);
                        Fe = Fe + Ne'*r*det(J)*W(i)*W(j);
                    end
                end
                M(LM(:,iel),LM(:,iel)) = M(LM(:,iel),LM(:,iel)) + Me;
                K(LM(:,iel),LM(:,iel)) = K(LM(:,iel),LM(:,iel)) + Ke;
                F(LM(:,iel),1) = F(LM(:,iel),1) + Fe;
            end
            
            M = sparse(M);
            K = sparse(K);
            
            Minv = diag(1./diag(M));
            T_dbc_Prev = [];
            T_dbc = [];
            DBCnodesI = [];
            for iel = Belem(:,1)'
                xe = xhat([LM(1,iel);LM(2,iel);LM(3,iel);LM(4,iel)],1);
                ye = xhat([LM(1,iel);LM(2,iel);LM(3,iel);LM(4,iel)],2);
                Fe_NBC = zeros(nen,1);
                He = zeros(nen,nen);
                ind = (Belem(:,1) == iel);
                for p = 1:sum(~isnan(Belem(ind,2:3)))
                    bndy = Belem(ind,p+1);
                    for q = 1:nip
                        [Ne,J] = CalcBNDInterpFunctions(xe,ye,nen,Z(q),bndy);
                        x = Ne*xe;
                        y = Ne*ye;
                        
                        beta = betafunc(x,y);
                        Tinf = Tinffunc(x,y);
                        qbarI = qfunc(x,y);
                        qbar = qbarI(bndy);
                        se   = sefunc(x,y);
                        
                        [fe,he] = evalNBC(Ne,bndy,BCCode,nen,beta,Tinf,qbar,se);
                        Fe_NBC = Fe_NBC + W(p)*fe*J;
                        He = He + W(p)*he*J;
                    end
                    [Ti_Prev,Ti,LN] = evalDBC(bndy,BCCode,xe,ye,DBCfunc,nen,L,t,dt);
                    T_dbc_Prev(LM(LN,iel),1) = Ti_Prev; clear Ti_Prev;
                    T_dbc(LM(LN,iel),1) = Ti; clear Ti;
                    DBCnodesI = [DBCnodesI;LM(LN,iel)]; %#ok<AGROW>
                end
                H(LM(:,iel),LM(:,iel)) = H(LM(:,iel),LM(:,iel)) + He;
                F(LM(:,iel),1) = F(LM(:,iel),1) + Fe_NBC;
            end
            
            F = sparse(F);
            T_dbc_Prev = T_dbc_Prev(unique(DBCnodesI),1);
            T_dbc = T_dbc(unique(DBCnodesI),1);
            
            %Remove DBCs
            Minv_hat = Minv(~DBCnodes,~DBCnodes);
            G_hat = K(~DBCnodes,~DBCnodes) + H(~DBCnodes,~DBCnodes);
            F_hat = F(~DBCnodes,1);
            
            M_bar = M(~DBCnodes,DBCnodes);
            G_bar = K(~DBCnodes,DBCnodes) + H(~DBCnodes,DBCnodes);
        end
        
        function  [R,J] =  ResJac(obj,U,t)
            %This function calcuates the residual and jacobian of the nonlinear heat
            %flow problem for the all models that require the full residual and
            %jacobian (i.e. FOM, GROM, PGROM, TPWL).
            %--------------------------------------------------------------------------
            %Inputs:
            %-------
            %GlobalMAT -
            %DBCvals   -
            %param     -
            %T_hat     -
            %t         -
            %dt        -
            %PHASEflag -
            %
            %Outputs:
            %--------
            %R       -
            %J       -
            %--------------------------------------------------------------------------
            
            LM = obj.HFparam.LM;
            
            nip = obj.HFparam.quad{1,1};
            Z = obj.HFparam.quad{2,1};
            W = obj.HFparam.quad{3,1};
            
            DBCnodes = obj.HFparam.BCs{3,1};
            Relem = obj.HFparam.BCs{5,1};
            Rbndy = obj.HFparam.BCs{6,1};
           
            sefunc = obj.HFparam.MATprop{4,1};
            
            nen = obj.HFparam.mesh{2,1};
            nnp = obj.HFparam.mesh{3,1};
            xhat = obj.HFparam.mesh{4,1};
            
            [Minv_hat,M_bar,F_hat,G_bar,G_hat, T_dbc_Prev, T_dbc] = evalFullGlobalMAT(obj,t,dt);
            
            
            F_rad = zeros(nnp,1);
            dF_rad = zeros(nnp,nnp);
            
            T(~DBCnodes,1) = T_hat;
            T(DBCnodes,1) = T_dbc;
            
            ind = 1;
            for iel = Relem(:)'
                f_rad = zeros(nen,1);
                df_rad = zeros(nen,nen);
                %     indT = Rnode(ind:ind+1,1);
                Ti = T(LM(:,iel),1);
                
                ind = ind + 1;
                xe = xhat([LM(1,iel);LM(2,iel);LM(3,iel);LM(4,iel)],1);
                ye = xhat([LM(1,iel);LM(2,iel);LM(3,iel);LM(4,iel)],2);
                for q = 1:nip
                    [Ne,jac] = CalcBNDInterpFunctions(xe,ye,nen,Z(q),Rbndy);
                    x = Ne*xe;
                    y = Ne*ye;
                    se   = sefunc(x,y);
                    f_rad = f_rad + W(q)*Ne'*se*(Ne*Ti)^4*jac;
                    df_rad = df_rad + 4*W(q)*se*(Ne*Ti)^3*(Ne'*Ne)*jac;
                end
                F_rad(LM(:,iel),1) = F_rad(LM(:,iel),1) + f_rad;
                dF_rad(LM(:,iel),LM(:,iel)) = dF_rad(LM(:,iel),LM(:,iel)) + df_rad;
            end
            
            R = Minv_hat*(F_hat - M_bar*(T_dbc-T_dbc_Prev)/dt - G_bar*T_dbc - G_hat*T_hat - F_rad(~DBCnodes,1));
            J = Minv_hat*(-G_hat - dF_rad(~DBCnodes,~DBCnodes));
        end
        
        
        
        
        
        function  [R,J] =  ResJacGNAT(obj)
        end
    end
    
end

