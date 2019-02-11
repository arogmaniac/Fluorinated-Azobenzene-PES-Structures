function [PlotSurf] = ConformysisResultsNoContour(atomcall,ang)
%UNTITLED reads in the results of G09 output files for tabulation and
%plotting the dat

fid=atomcall;
J=ang;
[m,n]=size(J);
A=zeros(1369,3)
for i=1:n;
    for j=1:n;
        name=[atomcall,'-',num2str(J(1,i)),'-',num2str(J(1,j)),'.log'];
        fid=fopen(name);
        string=fscanf(fid,'%s');
        myID='HF=-[0-9]+.[0-9]+';
        
        HFEnergies=regexp(string,myID,'match');
        
        HFEString=cell2mat(HFEnergies);
        
        if HFEString>0;
            HFES=mat2str(HFEString);
            Energy=regexp(HFES,'-\d+.\d+','match');
            string1=Energy{:};
            v(j)=str2num(string1);
        else
            v(j)=0;
        end
        X(j)=J(1,i);
        Y(j)=J(1,j);
        column=(i*37)-(37-j);
        A(column,:)=[X(j);Y(j);v(j)];
        
        fclose(fid);
    end;
    
end;
b=min(A(:,3));
B=A-b;
for i=1:37;
    for j=1:37;
        Row=(37*i)-(37-j);
        PlotSurf(i,j)=A(Row,3);
    end
end
F=PlotSurf.*627.509469;
yu=min(F);
ye=min(yu);
AA(1:37,1:37)=ye;
BB=F-AA;
for i=1:37;
    for j=1:37;
        if BB(i,j)>=50;
        BB(i,j)=50;
        end
    end
end
CC=[BB,BB;BB,BB]
DD=CC(19:55,19:55)
surf(DD,'FaceColor','interp','EdgeColor','none')
axis tight
set(get(gca,'XLabel'),'String','Rotation Angle \theta','FontName','Times New Roman','FontSize',14)
set(get(gca,'YLabel'),'String','Rotation Angle \phi','FontName','Times New Roman','FontSize',14)
set(get(gca,'ZLabel'),'String','Energy (kcal/mol)','FontName','Times New Roman','FontSize',14)
set(gca,'XTick',1:6:37)
set(gca,'XTickLabel',{'-180','-120','-60','0','60','120','180'},'FontSize',12,'FontName','Times New Roman')
set(gca,'YTick',1:6:37)
set(gca,'YTickLabel',{'-180','-120','-60','0','60','120','180'})
set(gca,'ZTick',0:10:30)
set(gca,'ZTickLabel',{'0','10','20','30'})



box on
colormap cool
camlight left; lighting phong

end

