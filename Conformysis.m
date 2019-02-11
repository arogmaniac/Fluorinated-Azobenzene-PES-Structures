function [Rotationfile] = Conformysis(A,atomsymb,p1,p2,p3,p4,...
    p12s,p12e,p34s,p34e,ang,moleculename)
%Code for the development of conformational analysis of N-C bond rotation
%in Azobenzenes for both phenyl moieties. Can be extended to other
%molecules as well.

%Variables
%A=[x,y,z] set of cartesian coordinates for complete molecule
%atomsymb= {'symb') set of atom symbols corresponding to each x,y,z locale
%p1,p2=(a,b,c) and (d,e,f) N-C bond xyz coordinates (moiety 1)
%p3,p4=(a,b,c) and (d,e,f) N-C bond xyz coordinates (moiety 2)
%..start and ..end are moiety subsets of A (moieties must be grouped in A)
%"ang" is the rotation angle range and increments about axis p#,p#
%"name" is the 

[i,j]=size(A);       %Assign i and j to rows and columns of A
subset12=A(p12s:p12e,:);
subsetloop12=subset12';
g=length(subsetloop12);
B=ones(1,g);
C=[subsetloop12;B];
subset12b=A(1:p12s-1,:);
subset12a=A(p12e+1:i,:);
DDD=cell2mat(atomsymb);
[iii,jjj]=size(atomsymb);
for zzz=1:iii;
    CCC=abs(DDD(1:zzz,:));
end
FFF=CCC;
for z=[p1];
    a=A(z,1);
    b=A(z,2);
    c=A(z,3);
    for o=[p2];
        ncaxis=[A(o,1)-A(z,1),A(o,2)-A(z,2),A(o,3)-A(z,3)];
        unitaxis=ncaxis./norm(ncaxis);
    end
    u=unitaxis(1,1);
    v=unitaxis(1,2);
    w=unitaxis(1,3);
    [d,q]=size(ang);
    for l=1:q;
        t=ang(1,l)/180*pi;
        Rmat=[(u^2)+(((v^2)+(w^2))*cos(t)) , (u*v*(1-cos(t)))-(w*sin(t)) , (u*w*(1-cos(t)))+(v*sin(t)) , (((a*((v^2)+(w^2)))-(u*((b*v)+(c*w))))*(1-cos(t)))+(((b*w)-(c*v))*sin(t));...
              (u*v*(1-cos(t)))+(w*sin(t)) , (v^2)+(((u^2)+(w^2))*cos(t)) , (v*w*(1-cos(t)))-(u*sin(t)) , (((b*((u^2)+(w^2)))-(v*((a*u)+(c*w))))*(1-cos(t)))+(((c*u)-(a*w))*sin(t));...
              (u*w*(1-cos(t)))-(v*sin(t)) , (v*w*(1-cos(t)))+(u*sin(t)) , (w^2)+(((u^2)+(v^2))*cos(t)) , (((c*((u^2)+(v^2)))-(w*((a*u)+(b*v))))*(1-cos(t)))+(((a*v)-(b*u))*sin(t));...
              0,0,0,1];
        Newcoords=Rmat*C;
        xyzform=Newcoords';
        backtoA=xyzform(:,1:3);
        Areconstruct=[subset12b;backtoA;subset12a];
        dg=ang(1,l);
        
        [ii,jj]=size(Areconstruct);
        subset34=Areconstruct(p34s:p34e,:);
        subsetloop34=subset34';
        gg=length(subsetloop34);
        BB=ones(1,gg);
        CC=[subsetloop34;BB];
        subset34b=Areconstruct(1:p34s-1,:);
        subset34a=Areconstruct(p34e+1:ii,:);

        for zz=[p3];
            aa=Areconstruct(zz,1);
            bb=Areconstruct(zz,2);
            cc=Areconstruct(zz,3);
            for oo=[p4];
                ncaxis1=[Areconstruct(oo,1)-Areconstruct(zz,1),Areconstruct(oo,2)-Areconstruct(zz,2),Areconstruct(oo,3)-Areconstruct(zz,3)];
                unitaxis1=ncaxis1./norm(ncaxis1);
            end
            uu=unitaxis1(1,1);
            vv=unitaxis1(1,2);
            ww=unitaxis1(1,3);
            [dd,qq]=size(ang);
            for ll=1:qq;
                tt=ang(1,ll)/180*pi;
                Rmat1=[(uu^2)+(((vv^2)+(ww^2))*cos(tt)) , (uu*vv*(1-cos(tt)))-(ww*sin(tt)) , (uu*ww*(1-cos(tt)))+(vv*sin(tt)) , (((aa*((vv^2)+(ww^2)))-(uu*((bb*vv)+(cc*ww))))*(1-cos(tt)))+(((bb*ww)-(cc*vv))*sin(tt));...
                       (uu*vv*(1-cos(tt)))+(ww*sin(tt)) , (vv^2)+(((uu^2)+(ww^2))*cos(tt)) , (vv*ww*(1-cos(tt)))-(uu*sin(tt)) , (((bb*((uu^2)+(ww^2)))-(vv*((aa*uu)+(cc*ww))))*(1-cos(tt)))+(((cc*uu)-(aa*ww))*sin(tt));...
                       (uu*ww*(1-cos(tt)))-(vv*sin(tt)) , (vv*ww*(1-cos(tt)))+(uu*sin(tt)) , (ww^2)+(((uu^2)+(vv^2))*cos(tt)) , (((cc*((uu^2)+(vv^2)))-(ww*((aa*uu)+(bb*vv))))*(1-cos(tt)))+(((aa*vv)-(bb*uu))*sin(tt));...
                        0,0,0,1];
                Newcoords1=Rmat1*CC;
                xyzform1=Newcoords1';
                backtoA1=xyzform1(:,1:3);
                Areconstruct1=[subset34b;backtoA1;subset34a];
                
                
                fg=ang(1,ll);
                name=[moleculename,'-',num2str(dg),'-',num2str(fg),'.com'];
                save(name,'-ascii');
                fid=fopen(name,'w');
                D=Areconstruct1;
                y=name;
                checkstring=[moleculename,'-',num2str(dg),'-',num2str(fg),'.chk'];
                fprintf(fid,'%%chk=%s\r\n',checkstring);
                theorystring=['# b3lyp/6-31++g(d''',',','p''',')']
                fprintf(fid,'%s\tnosymm\r\n\r\n',theorystring)
                randnum=['00000001']
                fprintf(fid,'%s\r\n\r\n',randnum)
                chargespin=['0 1']
                fprintf(fid,'%s\r\n',chargespin)
                rrr=numel(FFF);
                for xxx=1:rrr;
                fprintf(fid,'%s\t%10.5f\t%10.5f\t%10.5f\r\n',FFF(xxx,1),D(xxx,:));
                end
                fprintf(fid,'\r\n')
                fclose(fid);
            end

            end
        end
    end
end
  




