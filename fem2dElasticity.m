function fem2dElasticity(nvm)

    
 global lambda mu
 global nconn nv 
%--------Five point Guassian quadrature  ------------     

intx=[1/3 (6+sqrt(15))/21 (9-2*sqrt(15))/21 (6+sqrt(15))/21 (6-sqrt(15))/21 (9+2*sqrt(15))/21 (6-sqrt(15))/21];

inty=[1/3 (6+sqrt(15))/21 (6+sqrt(15))/21 (9-2*sqrt(15))/21 (6-sqrt(15))/21 (6-sqrt(15))/21 (9+2*sqrt(15))/21];

intw=[9/80 (155+sqrt(15))/2400 (155+sqrt(15))/2400 (155+sqrt(15))/2400 (155-sqrt(15))/2400 (155-sqrt(15))/2400 (155-sqrt(15))/2400];
%======================================================
max_it=1000;
tol=10^(-6);

% necessary constants
     lambda = 1;
     mu = 1;

% computational domain
     
     rmb=0;
     rme=1;
     rnb=0;
     rne=1;

% total number of nodes
     
     nvn = nvm;
     nvm1=nvm+1;
     nvn1=nvn+1;
     nv=nvm1*nvn1;
     
% total number of elements

     ne = nvm*nvn*2;
     
% get the mesh in each subdomain     
     
lijtk=(reshape(linspace(1,nv,nv),nvn1,nvm1))';
    

%  to get the local connectivity matrix


      ii=1;
      for i=1:nvm
          ip=i+1;
          for  j=1:nvn
               jp=j+1;
               nconn(1,ii)  =lijtk(i,j);
               nconn(2,ii)  =lijtk(ip,j);
               nconn(3,ii)  =lijtk(ip,jp);
               nconn(1,ii+1)=lijtk(i,j);
               nconn(2,ii+1)=lijtk(i,jp);
               nconn(3,ii+1)=lijtk(ip,jp);
               ii=ii+2;      
          end
      end


% initialize the stiffness matrix and load vector

      aK=sparse(2*nv,2*nv);
      b=zeros(2*nv,1);
	
% to obtain the coordinates

     for i=1:nvm1
         x(2,lijtk(i,:))=linspace(rnb ,rne ,nvn1);
     end
     for j=1:nvn1
         x(1,lijtk(:,j))=linspace(rmb ,rme ,nvm1);
     end
     

     
% calculate over the elements

     for k=1:ne
          n1 =nconn(1,k);
          n2 =nconn(2,k);
          n3 =nconn(3,k);
          n1v =n1+nv;
          n2v =n2+nv;
          n3v =n3+nv;
 	  
          % [k n1 n2 n3 n1v n2v n3v]
          x1=x(1,n1);
          y1=x(2,n1);
          x2=x(1,n2);
          y2=x(2,n2);
          x3=x(1,n3);
          y3=x(2,n3);

          b11=x2-x1;
          b12=x3-x1;
          b21=y2-y1;
          b22=y3-y1;
         
          detb=b11*b22-b12*b21;
          adetb=abs(detb);
          mydet(k)=adetb;
          d11=b22/detb;
          d12=-b21/detb;
          d21=-b12/detb;
          d22=b11/detb;
          w1x=-(d11+d12);
          w1y=-(d21+d22);
          w2x=d11;
          w2y=d21;
          w3x=d12;
          w3y=d22;
   
    %row 1
     aK(n1,n1)=aK(n1,n1)+adetb/2 * (lambda * w1x*w1x + 2*mu * (w1x*w1x + w1y/2*w1y/2*2));
     aK(n1,n2)=aK(n1,n2)+adetb/2 * (lambda * w1x*w2x + 2*mu * (w1x*w2x + w1y/2*w2y/2*2));
     aK(n1,n3)=aK(n1,n3)+adetb/2 * (lambda * w1x*w3x + 2*mu * (w1x*w3x + w1y/2*w3y/2*2));
     aK(n1,n1v)=aK(n1,n1v)+adetb/2 * (lambda * w1x*w1y + 2*mu * (w1y/2*w1x/2*2));
     aK(n1,n2v)=aK(n1,n2v)+adetb/2 * (lambda * w1x*w2y + 2*mu * (w1y/2*w2x/2*2));
     aK(n1,n3v)=aK(n1,n3v)+adetb/2 * (lambda * w1x*w3y + 2*mu * (w1y/2*w3x/2*2));
     %row 2
     aK(n2,n1)=aK(n2,n1)+adetb/2 * (lambda * w2x*w1x + 2*mu * (w2x*w1x + w2y/2*w1y/2*2));
     aK(n2,n2)=aK(n2,n2)+adetb/2 * (lambda * w2x*w2x + 2*mu * (w2x*w2x + w2y/2*w2y/2*2));
     aK(n2,n3)=aK(n2,n3)+adetb/2 * (lambda * w2x*w3x + 2*mu * (w2x*w3x + w2y/2*w3y/2*2));
     aK(n2,n1v)=aK(n2,n1v)+adetb/2 * (lambda * w2x*w1y + 2*mu * (w2y/2*w1x/2*2));
     aK(n2,n2v)=aK(n2,n2v)+adetb/2 * (lambda * w2x*w2y + 2*mu * (w2y/2*w2x/2*2));
     aK(n2,n3v)=aK(n2,n3v)+adetb/2 * (lambda * w2x*w3y + 2*mu * (w2y/2*w3x/2*2));
     %row 3
     aK(n3,n1)=aK(n3,n1)+adetb/2 * (lambda * w3x*w1x + 2*mu * (w3x*w1x + w3y/2*w1y/2*2));
     aK(n3,n2)=aK(n3,n2)+adetb/2 * (lambda * w3x*w2x + 2*mu * (w3x*w2x + w3y/2*w2y/2*2));
     aK(n3,n3)=aK(n3,n3)+adetb/2 * (lambda * w3x*w3x + 2*mu * (w3x*w3x + w3y/2*w3y/2*2));
     aK(n3,n1v)=aK(n3,n1v)+adetb/2 * (lambda * w3x*w1y + 2*mu * (w3y/2*w1x/2*2));
     aK(n3,n2v)=aK(n3,n2v)+adetb/2 * (lambda * w3x*w2y + 2*mu * (w3y/2*w2x/2*2));
     aK(n3,n3v)=aK(n3,n3v)+adetb/2 * (lambda * w3x*w3y + 2*mu * (w3y/2*w3x/2*2));
     %row 4
     aK(n1v,n1)=aK(n1v,n1)+adetb/2 * (lambda * w1y*w1x + 2*mu * (w1x/2*w1y/2*2));
     aK(n1v,n2)=aK(n1v,n2)+adetb/2 * (lambda * w1y*w2x + 2*mu * (w1x/2*w2y/2*2));
     aK(n1v,n3)=aK(n1v,n3)+adetb/2 * (lambda * w1y*w3x + 2*mu * (w1x/2*w3y/2*2));
     aK(n1v,n1v)=aK(n1v,n1v)+adetb/2 * (lambda * w1y*w1y + 2*mu * (w1x/2*w1x/2*2 + w1y*w1y));
     aK(n1v,n2v)=aK(n1v,n2v)+adetb/2 * (lambda * w1y*w2y + 2*mu * (w1x/2*w2x/2*2 + w1y*w2y));
     aK(n1v,n3v)=aK(n1v,n3v)+adetb/2 * (lambda * w1y*w3y + 2*mu * (w1x/2*w3x/2*2 + w1y*w3y));
     %row 5
     aK(n2v,n1)=aK(n2v,n1)+adetb/2 * (lambda * w2y*w1x + 2*mu * (w2x/2*w1y/2*2));
     aK(n2v,n2)=aK(n2v,n2)+adetb/2 * (lambda * w2y*w2x + 2*mu * (w2x/2*w2y/2*2));
     aK(n2v,n3)=aK(n2v,n3)+adetb/2 * (lambda * w2y*w3x + 2*mu * (w2x/2*w3y/2*2));
     aK(n2v,n1v)=aK(n2v,n1v)+adetb/2 * (lambda * w2y*w1y + 2*mu * (w2x/2*w1x/2*2 + w2y*w1y));
     aK(n2v,n2v)=aK(n2v,n2v)+adetb/2 * (lambda * w2y*w2y + 2*mu * (w2x/2*w2x/2*2 + w2y*w2y));
     aK(n2v,n3v)=aK(n2v,n3v)+adetb/2 * (lambda * w2y*w3y + 2*mu * (w2x/2*w3x/2*2 + w2y*w3y)); 
     %row 6
     aK(n3v,n1)=aK(n3v,n1)+adetb/2 * (lambda * w3y*w1x + 2*mu * (w3x/2*w1y/2*2));
     aK(n3v,n2)=aK(n3v,n2)+adetb/2 * (lambda * w3y*w2x + 2*mu * (w3x/2*w2y/2*2));
     aK(n3v,n3)=aK(n3v,n3)+adetb/2 * (lambda * w3y*w3x + 2*mu * (w3x/2*w3y/2*2));
     aK(n3v,n1v)=aK(n3v,n1v)+adetb/2 * (lambda * w3y*w1y + 2*mu * (w3x/2*w1x/2*2 + w3y*w1y));
     aK(n3v,n2v)=aK(n3v,n2v)+adetb/2 * (lambda * w3y*w2y + 2*mu * (w3x/2*w2x/2*2 + w3y*w2y));
     aK(n3v,n3v)=aK(n3v,n3v)+adetb/2 * (lambda * w3y*w3y + 2*mu * (w3x/2*w3x/2*2 + w3y*w3y));


      ointx=x1+b11*intx+b12*inty;
      ointy=y1+b21*intx+b22*inty;
      int1=0;
      int2=0;
      int3=0;
      int1v=0;
      int2v=0;
      int3v=0;
      for i=1:7
	  xxx=ointx(i);
	  yyy=ointy(i);
	  ff1=f1(xxx,yyy);
          int1=int1+ff1*(1-intx(i)-inty(i))*intw(i); 
          int2=int2+ff1*(  intx(i)        )*intw(i);   
          int3=int3+ff1*(          inty(i))*intw(i);
      ff2=f2(xxx,yyy);
          int1v=int1v+ff2*(1-intx(i)-inty(i))*intw(i);
          int2v=int2v+ff2*(  intx(i)        )*intw(i); 
          int3v=int3v+ff2*(          inty(i))*intw(i);
      end
      b(n1)=b(n1)+adetb*int1;
      b(n2)=b(n2)+adetb*int2;
      b(n3)=b(n3)+adetb*int3;
      b(n1v)=b(n1v)+adetb*int1v;
      b(n2v)=b(n2v)+adetb*int2v;
      b(n3v)=b(n3v)+adetb*int3v;                                                    
     end

  % to obtain the true solution
  
  %  utrue1=zeros(2*nv,2);
  % utrue2=zeros(2*nv,1);
  utrue1=zeros(nv,1);
  utrue2=zeros(nv,1);
  for iit= 1: nv
       xxx=x(1,iit);
       yyy=x(2,iit);
       futruexy=futrue(xxx,yyy);
       utrue1(iit)=futruexy(1);
       utrue2(iit)=futruexy(2);
  end
  %  utrue = [utrue1, utrue2];

% boundary treatment --------------------
  
      id=[];
      id=union(id,lijtk(1,:));
      id=union(id,lijtk(nvm1,:));
      %id= union(id,lijtk(:,1));
      %id = union(id,lijtk(:,nvn1));
      
      % zero dirichlet for u at x=0 and x=1, normal =0 for y=0 and y=1
      aK(id,:)=0;
      aK(:,id)=0;
      dd=size(id);
      iidd=max(dd(1),dd(2));
      for ddd=1:iidd
        iidd=id(ddd);
        aK(iidd,iidd)=1;
      end
      b(id)=0;	 
 
      id=[];
      %      id=union(id,lijtk(1,:));
      %id=union(id,lijtk(nvm1,:));
      id= union(id,lijtk(:,1));
      id = union(id,lijtk(:,nvn1));
      id=id+nv;
      % zero dirichlet for u at x=0 and x=1, normal =0 for y=0 and y=1
      aK(id,:)=0;
      aK(:,id)=0;
      dd=size(id);
      iidd=max(dd(1),dd(2));
      for ddd=1:iidd
        iidd=id(ddd);
        aK(iidd,iidd)=1;
      end
      b(id)=0;	 
      
   
      
     u0=aK\b;

     
     u=u0(1:nv);
     v=u0(nv+1:2*nv);
     
  
     % compute the error
     
     %     size(utrue1)
     %size(u)

      resl2u=norm(utrue1-u)/nvm;
     errmaxu=max(abs(utrue1-u));
        
     resl2v=norm(utrue2-v)/nvm;
     errmaxv=max(abs(utrue2-v));
     
     fprintf('L2u = %12.15f\n',resl2u); 
     fprintf('maxu = %12.15f\n',errmaxu); 

     fprintf('L2v = %12.15f\n',resl2v); 
     fprintf('maxv = %12.15f\n',errmaxv); 
     
stop
resl2=0;
resh1=0;

       for k=1:ne
           n1=nconn(1,k);
           n2=nconn(2,k);
           n3=nconn(3,k);
           %n1v=n1+nv;
           %n2v=n2+nv;
           %n3v=n3+nv;
           x1=x(1,n1);
           y1=x(2,n1);
           x2=x(1,n2);
           y2=x(2,n2);
           x3=x(1,n3);
           y3=x(2,n3);
           x1v=x(1,n1);
           y1v=x(2,n1);
           x2v=x(1,n2);
           y2v=x(2,n2);
           x3v=x(1,n3);
           y3v=x(2,n3);
           b11=x2-x1;
           b12=x3-x1;
           b21=y2-y1;
           b22=y3-y1;

        ointx=x1+b11*intx+b12*inty;
        ointy=y1+b21*intx+b22*inty;

A=[x1 y1 1;x2 y2 1;x3 y3 1;x1v y1v 1;x2v y2v 1;x3v y3v 1];
u1 =A\[u(n1); 0    ; 0     ; 0      ; 0      ; 0     ];
u2 =A\[0    ; u(n2); 0     ; 0      ; 0      ; 0     ];
u3 =A\[0    ; 0    ; u(n3) ; 0      ; 0      ; 0     ];
u1v=A\[0    ; 0    ; 0     ; u(n1v) ; 0      ; 0     ];
u2v=A\[0    ; 0    ; 0     ; 0      ; u(n2v) ; 0     ];
u3v=A\[0    ; 0    ; 0     ; 0      ; 0      ; u(n3v)];

uh=[ointx' ointy' [1; 1; 1; 1; 1; 1; 1]]*(A\[u(n1); u(n2); u(n3); u(n1v); u(n2v); u(n3v)]);
uh=uh';
uh1= u1(1)+u2(1)+u3(1)+u1v(1)+u2v(1)+u3v(1);
uh2= u1(2)+u2(2)+u3(2)+u1v(2)+u2v(2)+u3v(2);


        ll2=0;
        hh1=0;
       for i=1:7
       ll2 = ll2+(futrue(ointx(i),ointy(i))-uh(i)).^2*intw(i);
       hh1 = hh1+((fux(ointx(i),ointy(i))-uh1).^2+(fuy(ointx(i),ointy(i))-uh2).^2)*intw(i);
       

       end            
        ll2sum = ll2(1) + ll2(2);
        hh1sum = hh1(1)+hh1(2);
    
resl2= resl2+mydet(k)*ll2sum;
resh1= resh1+mydet(k)*hh1sum;

end

resl2=sqrt(resl2);
resh1=sqrt(resh1);
err=[resl2 resh1]

function result=f1(x,y)
    global lambda mu
    %    mu =1;
    %  result=-(-pi*pi*mu*sin(pi*x)*0.01 + pi*pi*mu/2*sin(pi*x)*0.01);
  
    result=(lambda+mu)*pi^2*sin(pi*x);
return;

function result=f2(x,y)
global lambda mu
%    mu =1;
%    result = -pi*pi*pi*mu/2*y*cos(pi*x)*0.01;
    result=(lambda+mu)*pi^2*sin(pi*y);

return;

function result=futrue(x,y)
 
%      result=[0.01*sin(pi*x) 0.01*(-pi*y*cos(pi*x))];
 
    result=[sin(pi*x) sin(pi*y)];
      return
 
      

function result=fux(x,y)

 result = [0.01*pi*cos(pi*x) 0.01*pi*pi*y*sin(pi*x)];
 
 return

 


function result=fuy(x,y)

 result = [0 -0.01*pi*cos(pi*x)];
 
 return

