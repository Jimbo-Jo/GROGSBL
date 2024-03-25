function [Pm,search_area] = GROGSBL(Y,dd,lambda,search_area,etc)

[~,~,D_svd] = svd(Y,'econ');
Y = Y*D_svd(:,1:etc);
[M,T]=size(Y);

A_theta = exp(-1i*2*pi*(0:M-1)'*(dd/lambda)*cosd(search_area));
B_theta = 1i*2*pi*(0:M-1)'*(dd/lambda)*sind(search_area).*A_theta;
BHB = B_theta'*B_theta;

%parameters
maxiter = 500;
maxI1 = 5;
tol = 1e-2;
minigap = 1;

%initialization
beta=(M*T)/(0.01*norm(Y, 'fro')^2);
Gamma=mean(abs(A_theta'*Y),2);

converged = false;
iter = 0;

while ~converged
   if iter<maxI1
      iter = iter + 1;
      
      %calculate mu and Sigma
      sigmaY = eye(M)/beta + A_theta*diag(Gamma)*A_theta';
      mu = diag(Gamma)*A_theta'/sigmaY*Y;
      sigma = diag(Gamma) - diag(Gamma)*A_theta'/sigmaY*A_theta*diag(Gamma);
      
      %calculate gamma
      Gamma = diag(abs(mu*mu'/T))+real(diag(sigma));
      
      %update beta
      [~,locs] = findpeaks( abs(diag(mu*mu')), 'SortStr','descend','Npeaks',etc);
      A_active = A_theta(:,sort(locs));
      Proj = A_active*pinv(A_active);
      beta  = (M-etc)/(trace( (eye(M)-Proj)*(Y*Y'/T) ));
      
      %grid refine
      [~, index_temp] = sort(Gamma, 'descend');
      [index,~] = sort(index_temp(1:etc));
      t=0;
      for j = 1:length(index)
          ii = index(j)+t;
          if ii==1
              search_area = [search_area(1) (search_area(1)+search_area(2))/2 search_area(2:end)];
              Gamma = [Gamma(1)/2;Gamma(1)/2;Gamma(2:end)];
          elseif ii == length(search_area)
               search_area = [search_area(1:end-1) (search_area(end-1)+search_area(end))/2 search_area(end)];
               Gamma = [Gamma(1:end-1);Gamma(end)/2;Gamma(end)/2];
          else
               midright = (search_area(ii+1)+search_area(ii))/2;
               midleft = (search_area(ii-1)+search_area(ii))/2;
               if (search_area(ii+1)-search_area(ii))>2*minigap || (search_area(ii)-search_area(ii-1))>2*minigap                   
                   t=t+1;
                   P = real(conj(BHB(ii,ii)).*(mu(ii,:)*mu(ii,:)'+T*sigma(ii,ii)));
                   v = 0;
                   for tt = 1:T
                       v = v+real(conj(mu(ii,tt)).*(B_theta(:,ii)'*(Y(:,tt)-A_theta*mu(:,tt))));
                   end
                   v = v-T*real(diag(B_theta(:,ii)'*A_theta*sigma(:,ii)));
                   beta0 = v/P;
                   
                   if beta0 >= 0
                       search_area = [search_area(1:ii) midright search_area(ii+1:end)];
                       Gamma = [Gamma(1:ii-1);Gamma(ii)/2;Gamma(ii)/2;Gamma(ii+1:end)];
                   else
                       search_area = [search_area(1:ii-1) midleft search_area(ii:end)];
                       Gamma = [Gamma(1:ii-1);Gamma(ii)/2;Gamma(ii)/2;Gamma(ii+1:end)];                     
                   end             
               end
          end
          
          A_theta=exp(-1i*2*pi*(0:M-1)'*(dd/lambda)*cosd(search_area));
          B_theta = 1i*2*pi*(0:M-1)'*(dd/lambda)*sind(search_area).*A_theta;
          BHB = B_theta'*B_theta;
          sigmaY = eye(M)/beta + A_theta*diag(Gamma)*A_theta';
          mu = diag(Gamma)*A_theta'/sigmaY*Y;
          sigma = diag(Gamma) - diag(Gamma)*A_theta'/sigmaY*A_theta*diag(Gamma);
          
      end
   else
   
      iter = iter + 1;
      Gamma_last = Gamma;
   
      %calculate mu and Sigma
      sigmaY = eye(M)/beta + A_theta*diag(Gamma)*A_theta';
      mu = diag(Gamma)*A_theta'/sigmaY*Y;
      sigmann = Gamma - Gamma.^2.*diag(A_theta'/sigmaY*A_theta);
      
      %calculate gamma
      Gamma = diag(abs(mu*mu'/T))+real(sigmann);
      
      %update beta
      [~,locs] = findpeaks( abs(diag(mu*mu')), 'SortStr','descend','Npeaks',etc);
      A_active = A_theta(:,sort(locs));
      Proj = A_active*pinv(A_active);
      beta  = (M-etc)/(trace( (eye(M)-Proj)*(Y*Y'/T) ));

    % stopping criteria
      erro=norm(Gamma - Gamma_last)/norm(Gamma_last);
      if erro < tol || iter >= maxiter
          converged = true;
      end
   end   
end

%offgrid DOA 
      [~, index_temp] = sort(Gamma, 'descend');
      index = index_temp(1:etc);
      
      for j = 1:length(index)
          ii = index(j);   
          if ii ~= 1&&ii ~=length(search_area)
              if Gamma(ii+1)>=Gamma(ii-1)
                  Aii = [A_theta(:,ii) A_theta(:,ii+1)];
                  Gammaii = diag([Gamma(ii);Gamma(ii+1)]);
              else
                  Aii = [A_theta(:,ii-1) A_theta(:,ii)];
                  Gammaii = diag([Gamma(ii-1);Gamma(ii)]);                 
              end
              sigmaYii = sigmaY - Aii*Gammaii*Aii';
              
              search_gap = 0.1;
              search_temp = search_area(ii-1):search_gap:search_area(ii+1);
              A_temp = exp(-1i*2*pi*(0:M-1)'*(dd/lambda)*cosd(search_temp));
              Zm = diag(A_temp'/sigmaYii*A_temp);
              qm_temp = abs(A_temp'/sigmaYii*Y).^2;
              qm = sum(qm_temp,2);
              
              Gammam = (qm-T*Zm)./(T*Zm.^2);
              L = -T*log(1+Gammam.*Zm) + qm./(1./Gammam+Zm);
              
              [~,indexmax] = max(L);
              search_area(ii) =  search_temp(indexmax);
          end
      end
      
Pm=mean(mu.*conj(mu),2);
end

