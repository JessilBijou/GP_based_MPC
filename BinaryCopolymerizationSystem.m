classdef BinaryCopolymerizationSystem
   
    
    properties
        k;          % Kinetic constants
        monomers;   % Names of the monomers
        mw;         % Mass weights
        r;          % Reactivity ratios
        rho;        % Densities
    end
    
    methods
        function obj = BinaryCopolymerizationSystem(k,monomers,mw,r,rho)
              
            obj.k = k;
            obj.monomers = monomers;
            obj.mw = mw;
            obj.r = r;
            obj.rho = rho;
            
        end
        
        function [kp_star,kfs_star,kfm_star,kt_star,r] = pseudo_kinetic_constants(obj,T,A,B,wp)
            
            
            % Evaluating the kinetic constants at the current temperature:
            kp = [obj.k.p.AA(T) obj.k.p.BB(T)];     %[L/mol/s]
            
            % If depropagation is active we need to re evaluate a few 
            % parameters:
            if isfield(obj.k,'dp')
                % Evaluate effective propagation constants:
                kdp = [obj.k.dp.AA(T,wp) obj.k.dp.BB(T,wp)];  %[1/s]
                kp =  [kp(1)-kdp(1)/A kp(2)-kdp(2)/B];        %[L/mol/s]  
            end
            
            r = [obj.r.A(T) obj.r.B(T)];                  %[-]
            kxp = kp./r;                            %[L/mol/s]
            kfs = [obj.k.fs.AA(T) obj.k.fs.BB(T)];  %[L/mol/s]
            kfm = [obj.k.fm.AA(T) obj.k.fm.BB(T)];  %[L/mol/s]
            kxfm = kfm./r;                          %[L/mol/s]
            kt =  [obj.k.t.AA(T) obj.k.t.BB(T)];    %[L/mol/s]
      
            % Pseudo-Kinetic constants:             
            p_A = kxp(2)*A/(kxp(2)*A+kxp(1)*B);         %[-]
            p_B = kxp(1)*B/(kxp(2)*A+kxp(1)*B);         %[-]
            X_A = A/(A+B);                              %[-]
            X_B = B/(A+B);                              %[-]
            kp_star.A = (kp(1)*p_A+kxp(2)*p_B);         %[L/mol/s]      
            kp_star.B = (kp(2)*p_B+kxp(1)*p_A);         %[L/mol/s]
            kp_star.tot = kp_star.A*X_A+kp_star.B*X_B;  %[L/mol/s]
            kt_star.tot = exp(p_A*log(kt(1))+p_B*log(kt(2)));     %[L/mol/s]
            kfs_star = kfs(1)*p_A+kfs(2)*p_B;           %[L/mol/s]
            kfm_star = (kfm(1)*X_A+kxfm(1)*X_B)*p_A+(kfm(2)*X_B+kxfm(2)*X_A)*p_B; %[L/mol/s]
            
            % Kinetic constants of bimolecular terminations:
            ktc = kt_star.tot*[1-obj.k.t.a, 1-obj.k.t.b];
            kxtc = kt_star.tot*(1-obj.k.t.c);
            ktd = kt_star.tot*[obj.k.t.a, obj.k.t.b];
            kxtd = kt_star.tot*obj.k.t.c;
            kt_star.c = ktc(1)*p_A^2+ktc(2)*p_B^2+2*kxtc*p_A+p_B;
            kt_star.d = ktd(1)*p_A^2+ktd(2)*p_B^2+2*kxtd*p_A+p_B;
            
        end
    end
end

