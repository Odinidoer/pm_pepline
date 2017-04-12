clc;close all;clear all;
inf=importdata('to_figure_classify.txt');
da=inf.data;test=inf.textdata;
type1=unique(da(:,1));
type=type1';
lab={'A' 'B' 'C' 'D' 'E' 'F' 'G'};
%color1=[0 128 0;128 0 128;220 20 60;0 0 205;184 134 11;255 255 0;0 255 255;255 192 203;148 0 211;0 128 0;128 0 128;220 20 60;0 0 205;184 134 11;255 255 0;0 255 255;255 192 203;148 0 211]/255;
color1=[255 215 0;255 255 0;0 255 0;0 255 255;148 0 211;255 20 147]/255;
x_max=max(da(:,3));
x_total=x_max/0.6;
width=0.7;y_start=1;
for j=1:length(type);
    tar1=da(find(da(:,1) == type(j)),:);
    test1=test(find(da(:,1) == type(j)));
    rectangle('position',[x_total*0.7,y_start,0.01*x_total,length(tar1(:,2))-1+width],'facecolor',color1(j,:),'Curvature',[1,.5/length(tar1(:,2))]);hold on; %对分类进行标注
    text(x_total*0.74,y_start+length(tar1(:,2))/2,lab{j},'HorizontalAlignment','left','Interpreter','none','fontsize',10);hold on;
    
    for i=1:length(tar1(:,2))
        kk1=text(-0.002,y_start+width/2,test1{i},'HorizontalAlignment','right','Interpreter','none');
        kk2=text(tar1(i,3)+0.02*x_max,y_start+width/2,num2str(tar1(i,2)),'HorizontalAlignment','left','Interpreter','none','fontsize',9); %基因数量
        rectangle('position',[0,y_start,tar1(i,3),width],'facecolor',color1(j,:));hold on;
        y_start=y_start+1;
    end
end
%添加x轴刻度
plot(linspace(0,x_total*0.7,20),zeros(1,20),'k');hold on; %x轴下边线
jianju=0.3;
for i=0:0.05:x_total*0.7
    plot(i*ones(1,20),linspace(-jianju,0,20),'k');hold on;
    text(i,-jianju-0.8,num2str(100*i),'HorizontalAlignment','center');
end
text(x_max*0.4,-3,'Percent of Protein (%)','HorizontalAlignment','left','Interpreter','none','fontsize',10);hold on;
%plot(x_max*ones(1,20),linspace(-jianju,0,20),'k');hold on;
set(gca,'xtick',[],'ytick',[],'xcolor',[254 254 254]/255,'ycolor',[254 254 254]/255);
axis([-0.3*x_total x_total*0.8 -1.2 y_start]);
saveas(gca,'kegg_classification.pdf');close;
