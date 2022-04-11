% Monte-carlo random walk simulation
% Here adsorption of non-ionic surfactants on organoclays in drilling fluid
% investigated by molecular descriptiors and Monte Carlo random walk
% simulations
clc;clear;
% boundary points
p1=[20,20];p2=[20,80];
p3=[80,20];p4=[80,80];
nB=100; % when the particle is out of boundary B, it is variable for laying in the B.
nt=200; % number of steps
nl=50; % number of molecules
ns=3;   % size of step
np=0.5; np2=0.7*np; np3=0.5*np; % Adsorption probability
%variable to store last position of particles
x=zeros(nl);x_left=zeros(1,nl);x_right=zeros(1,nl);
y=zeros(nl);y_down=zeros(1,nl);y_top=zeros(1,nl);

count1=0; count2=0; count3=0;
count11=0; count12=0; count13=0;
% adsorption of molecules with 2 chains
figure;
clf; hold on;
axis([-5 105 -5 105]);
line(p1,p2,'color','red')
line(p2,p4,'color','red')
line(p4,p3,'color','red')
line(p3,p1,'color','red')
for i=20:80
    for j=20:80
        x=i;
        y=j;
        hold on;
        plot(x,y,'r--.', ...
            'MarkerSize',1)
    end
end
P=zeros(nl,1);
for l=1:nl
    % set the initial point
    xy0=randi([0,100],1,2);
    if 20<xy0(1,1) && xy0(1,1)<80
        num=randi([0,100],1);
        xy0(1,1)=num;
    end
    if 20<xy0(1,2) && xy0(1,2)<80
        num=randi([0,100],1);
        xy0(1,2)=num;
    end
    xy=zeros(nt,2);
    xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
    % generate the random size of steps
    deltax=randi([-ns,ns],1,nt);
    deltay=randi([-ns,ns],1,nt);
    for step=1:nt-1
        % Walk in the x direction.
        if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
            p=rand(1);
            P(step)=p;
            if np>p
                xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                break;
            end
        end
        xy(step+1,1)=xy(step,1)+deltax(step);
        % apply periodic boundary condition
        if xy(step+1,1)>=0 && xy(step+1,1)<=nB
            xy(step+1,1)=xy(step+1,1);
        end
        if xy(step+1,1)>nB
            xy(step+1,1)=xy(step+1,1)-nB;
        end
        if xy(step+1,1)<0
            xy(step+1,1)=xy(step+1,1)+nB;
        end
	    % Walk in the y direction.
	    xy(step+1,2)=xy(step,2)+deltay(step);
        if xy(step+1,2)>=0 && xy(step+1,2)<=nB
            xy(step+1,2)=xy(step+1,2);
        end
        if xy(step+1,2)>nB
            xy(step+1,2)=xy(step+1,2)-nB;
        end
        if xy(step+1,2)<0
            xy(step+1,2)=xy(step+1,2)+nB;
        end
    end
    x(l)=xy(nt,1); y(l)=xy(nt,2);
    x_left(l)=x(l)-1;x_right(l)=x(l)+1;
    if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
        count11=count11+1;
    end
    if x_left(l)>20 && x_left(l)<80
        count1=count1+1;
    end
    if x_right(l)>20 && x_right(l)<80
        count1=count1+1;
    end
end
acc1=count1*100/3600;
ac1=count11*100/nl;
hold on
plot(x_left,y,'b.')
plot(x_right,y,'b.')


% adsorption of molecules with 3 chains
figure;
clf; hold on;
axis([-5 105 -5 105]);
line(p1,p2,'color','red')
line(p2,p4,'color','red')
line(p4,p3,'color','red')
line(p3,p1,'color','red')
for i=20:80
    for j=20:80
        x=i;
        y=j;
        hold on;
        plot(x,y,'r--.', ...
            'MarkerSize',1)
    end
end
P=zeros(nl,1);
for l=1:nl
    % set the initial point
    xy0=randi([0,100],1,2);
    if 20<xy0(1,1) && xy0(1,1)<80
        num=randi([0,100],1);
        xy0(1,1)=num;
    end
    if 20<xy0(1,2) && xy0(1,2)<80
        num=randi([0,100],1);
        xy0(1,2)=num;
    end
    xy=zeros(nt,2);
    xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
    % generate the random size of steps
    deltax=randi([-ns,ns],1,nt);
    deltay=randi([-ns,ns],1,nt);
    for step=1:nt-1
        % Walk in the x direction.
        if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
            p=rand(1);
            P(step)=p;
            if np2>p
                xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                break;
            end
        end
        xy(step+1,1)=xy(step,1)+deltax(step);
        % apply periodic boundary condition
        if xy(step+1,1)>=0 && xy(step+1,1)<=nB
            xy(step+1,1)=xy(step+1,1);
        end
        if xy(step+1,1)>nB
            xy(step+1,1)=xy(step+1,1)-nB;
        end
        if xy(step+1,1)<0
            xy(step+1,1)=xy(step+1,1)+nB;
        end
	    % Walk in the y direction.
	    xy(step+1,2)=xy(step,2)+deltay(step);
        if xy(step+1,2)>=0 && xy(step+1,2)<=nB
            xy(step+1,2)=xy(step+1,2);
        end
        if xy(step+1,2)>nB
            xy(step+1,2)=xy(step+1,2)-nB;
        end
        if xy(step+1,2)<0
            xy(step+1,2)=xy(step+1,2)+nB;
        end
    end
    x(l)=xy(nt,1); y(l)=xy(nt,2);
    x_left(l)=x(l)-1;x_right(l)=x(l)+1;
    y_down(l)=y(l)-1;
    if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
        count12=count12+1;
    end
    if x_left(l)>20 && x_left(l)<80
        count2=count2+1;
    end
    if x_right(l)>20 && x_right(l)<80
        count2=count2+1;
    end
    if y_down(l)>20 && y_down(l)<80
        count2=count2+1;
    end
end
acc2=count2*100/3600;
ac2=count12*100/nl;
hold on
plot(x_left,y,'b.')
plot(x_right,y,'b.')
plot(x,y_down,'b.')

% adsorption of molecules with 4 chains
figure;
clf; hold on;
axis([-5 105 -5 105]);
line(p1,p2,'color','red')
line(p2,p4,'color','red')
line(p4,p3,'color','red')
line(p3,p1,'color','red')
for i=20:80
    for j=20:80
        x=i;
        y=j;
        hold on;
        plot(x,y,'r--.', ...
            'MarkerSize',1)
    end
end
P=zeros(nl,1);
for l=1:nl
    % set the initial point
    xy0=randi([0,100],1,2);
    if 20<xy0(1,1) && xy0(1,1)<80
        num=randi([0,100],1);
        xy0(1,1)=num;
    end
    if 20<xy0(1,2) && xy0(1,2)<80
        num=randi([0,100],1);
        xy0(1,2)=num;
    end
    xy=zeros(nt,2);
    xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
    % generate the random size of steps
    deltax=randi([-ns,ns],1,nt);
    deltay=randi([-ns,ns],1,nt);
    for step=1:nt-1
        % Walk in the x direction.
        if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
            p=rand(1);
            P(step)=p;
            if np3>p
                xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                break;
            end
        end
        xy(step+1,1)=xy(step,1)+deltax(step);
        % apply periodic boundary condition
        if xy(step+1,1)>=0 && xy(step+1,1)<=nB
            xy(step+1,1)=xy(step+1,1);
        end
        if xy(step+1,1)>nB
            xy(step+1,1)=xy(step+1,1)-nB;
        end
        if xy(step+1,1)<0
            xy(step+1,1)=xy(step+1,1)+nB;
        end
	    % Walk in the y direction.
	    xy(step+1,2)=xy(step,2)+deltay(step);
        if xy(step+1,2)>=0 && xy(step+1,2)<=nB
            xy(step+1,2)=xy(step+1,2);
        end
        if xy(step+1,2)>nB
            xy(step+1,2)=xy(step+1,2)-nB;
        end
        if xy(step+1,2)<0
            xy(step+1,2)=xy(step+1,2)+nB;
        end
    end
    x(l)=xy(nt,1); y(l)=xy(nt,2);
    x_left(l)=x(l)-1;x_right(l)=x(l)+1;
    y_down(l)=y(l)-1;y_top(l)=y(l)+1;
    if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
        count13=count13+1;
    end
    if x_left(l)>20 && x_left(l)<80
        count3=count3+1;
    end
    if x_right(l)>20 && x_right(l)<80
        count3=count3+1;
    end
    if y_down(l)>20 && y_down(l)<80
        count3=count3+1;
    end
    if y_top(l)>20 && y_top(l)<80
        count3=count3+1;
    end
end
acc3=count3*100/3600;
ac3=count13*100/nl;
hold on
plot(x_left,y,'b.')
plot(x_right,y,'b.')
plot(x,y_down,'b.')
plot(x,y_top,'b.')

% percent depend on the number of steps
accuracy11=zeros(1,5);accuracy12=zeros(1,5);accuracy13=zeros(1,5);
accuracy21=zeros(1,5);accuracy22=zeros(1,5);accuracy23=zeros(1,5);
for iter=1:5
    size=zeros(1,5);
    ns=5;np=0.4;nl=100;np2=0.7*np; np3=0.5*np;
    acc1=zeros(1,5);acc2=zeros(1,5);acc3=zeros(1,5);
    ac1=zeros(1,5);ac2=zeros(1,5);ac3=zeros(1,5);
    for i=1:5
        count1=0; count2=0; count3=0;
        count11=0; count12=0; count13=0;
        nt=(i-1)*50;
        if i==1
            nt=nt+25;
        end
        for l=1:nl
            % set the initial point
            xy0=randi([0,100],1,2);
            if 20<xy0(1,1) && xy0(1,1)<80
                num=randi([0,100],1);
                xy0(1,1)=num;
            end
            if 20<xy0(1,2) && xy0(1,2)<80
                num=randi([0,100],1);
                xy0(1,2)=num;
            end
            xy=zeros(nt,2);
            xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
            % generate the random size of steps
            deltax=randi([-ns,ns],1,nt);
            deltay=randi([-ns,ns],1,nt);
            for step=1:nt-1
                % Walk in the x direction.
                if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
                    p=rand(1);
                    P(step)=p;
                    if np>p
                        xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                        break;
                    end
                end
                xy(step+1,1)=xy(step,1)+deltax(step);
                % apply periodic boundary condition
                if xy(step+1,1)>=0 && xy(step+1,1)<=nB
                    xy(step+1,1)=xy(step+1,1);
                end
                if xy(step+1,1)>nB
                    xy(step+1,1)=xy(step+1,1)-nB;
                end
                if xy(step+1,1)<0
                    xy(step+1,1)=xy(step+1,1)+nB;
                end
	            % Walk in the y direction.
	            xy(step+1,2)=xy(step,2)+deltay(step);
                if xy(step+1,2)>=0 && xy(step+1,2)<=nB
                    xy(step+1,2)=xy(step+1,2);
                end
                if xy(step+1,2)>nB
                    xy(step+1,2)=xy(step+1,2)-nB;
                end
                if xy(step+1,2)<0
                    xy(step+1,2)=xy(step+1,2)+nB;
                end
            end
            x(l)=xy(nt,1); y(l)=xy(nt,2);
            x_left(l)=x(l)-1;x_right(l)=x(l)+1;
            if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
                count11=count11+1;
            end
            if x_left(l)>20 && x_left(l)<80
                count1=count1+1;
            end
            if x_right(l)>20 && x_right(l)<80
                count1=count1+1;
            end
        end
        acc1(i)=count1*100/3600;
        ac1(i)=count11*100/nl;
        size(i)=nt;
            
        for l=1:nl
            % set the initial point
            xy0=randi([0,100],1,2);
            if 20<xy0(1,1) && xy0(1,1)<80
                num=randi([0,100],1);
                xy0(1,1)=num;
            end
            if 20<xy0(1,2) && xy0(1,2)<80
                num=randi([0,100],1);
                xy0(1,2)=num;
            end
            xy=zeros(nt,2);
            xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
            % generate the random size of steps
            deltax=randi([-ns,ns],1,nt);
            deltay=randi([-ns,ns],1,nt);
            for step=1:nt-1
                % Walk in the x direction.
                if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
                    p=rand(1);
                    P(step)=p;
                    if np2>p
                        xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                        break;
                    end
                end
                xy(step+1,1)=xy(step,1)+deltax(step);
                % apply periodic boundary condition
                if xy(step+1,1)>=0 && xy(step+1,1)<=nB
                    xy(step+1,1)=xy(step+1,1);
                end
                if xy(step+1,1)>nB
                    xy(step+1,1)=xy(step+1,1)-nB;
                end
                if xy(step+1,1)<0
                    xy(step+1,1)=xy(step+1,1)+nB;
                end
	            % Walk in the y direction.
	            xy(step+1,2)=xy(step,2)+deltay(step);
                if xy(step+1,2)>=0 && xy(step+1,2)<=nB
                    xy(step+1,2)=xy(step+1,2);
                end
                if xy(step+1,2)>nB
                    xy(step+1,2)=xy(step+1,2)-nB;
                end
                if xy(step+1,2)<0
                    xy(step+1,2)=xy(step+1,2)+nB;
                end
            end
            x(l)=xy(nt,1); y(l)=xy(nt,2);
            x_left(l)=x(l)-1;x_right(l)=x(l)+1;
            y_down(l)=y(l)-1;
            if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
                count12=count12+1;
            end
            if x_left(l)>20 && x_left(l)<80
                count2=count2+1;
            end
            if x_right(l)>20 && x_right(l)<80
                count2=count2+1;
            end
            if y_down(l)>20 && y_down(l)<80
                count2=count2+1;
            end
        end
        acc2(i)=count2*100/3600;
        ac2(i)=count12*100/nl;
            
        for l=1:nl
            % set the initial point
            xy0=randi([0,100],1,2);
            if 20<xy0(1,1) && xy0(1,1)<80
                num=randi([0,100],1);
                xy0(1,1)=num;
            end
            if 20<xy0(1,2) && xy0(1,2)<80
                num=randi([0,100],1);
                xy0(1,2)=num;
            end
            xy=zeros(nt,2);
            xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
            % generate the random size of steps
            deltax=randi([-ns,ns],1,nt);
            deltay=randi([-ns,ns],1,nt);
            for step=1:nt-1
                % Walk in the x direction.
                if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
                    p=rand(1);
                    P(step)=p;
                    if np3>p
                        xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                        break;
                    end
                end
                xy(step+1,1)=xy(step,1)+deltax(step);
                % apply periodic boundary condition
                if xy(step+1,1)>=0 && xy(step+1,1)<=nB
                    xy(step+1,1)=xy(step+1,1);
                end
                if xy(step+1,1)>nB
                    xy(step+1,1)=xy(step+1,1)-nB;
                end
                if xy(step+1,1)<0
                    xy(step+1,1)=xy(step+1,1)+nB;
                end
	            % Walk in the y direction.
	            xy(step+1,2)=xy(step,2)+deltay(step);
                if xy(step+1,2)>=0 && xy(step+1,2)<=nB
                    xy(step+1,2)=xy(step+1,2);
                end
                if xy(step+1,2)>nB
                    xy(step+1,2)=xy(step+1,2)-nB;
                end
                if xy(step+1,2)<0
                    xy(step+1,2)=xy(step+1,2)+nB;
                end
            end
            x(l)=xy(nt,1); y(l)=xy(nt,2);
            x_left(l)=x(l)-1;x_right(l)=x(l)+1;
            y_down(l)=y(l)-1;y_top(l)=y(l)+1;
            if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
                count13=count13+1;
            end
            if x_left(l)>20 && x_left(l)<80
                count3=count3+1;
            end
            if x_right(l)>20 && x_right(l)<80
                count3=count3+1;
            end
            if y_down(l)>20 && y_down(l)<80
                count3=count3+1;
            end
            if y_top(l)>20 && y_top(l)<80
                count3=count3+1;
            end
        end
        acc3(i)=count3*100/3600;
        ac3(i)=count13*100/nl;
    end
    accuracy11=accuracy11+acc1;
    accuracy21=accuracy21+ac1;
    accuracy12=accuracy12+acc2;
    accuracy22=accuracy22+ac2;
    accuracy13=accuracy13+acc3;
    accuracy23=accuracy23+ac3;
end
accuracy11=accuracy11/5;
accuracy21=accuracy21/5;
accuracy12=accuracy12/5;
accuracy22=accuracy22/5;
accuracy13=accuracy13/5;
accuracy23=accuracy23/5;
figure;
hold on
plot(size,accuracy11,'bo-');
plot(size,accuracy12,'Color','#D95319','marker','x');
plot(size,accuracy13,'Color','#77AC30','marker','^');
xlabel('Number of steps')
ylabel('% Organoclay site occupied')
legend('NPGE(2-chain)','TMPE(3-chain)','PEE(4-chain)')
figure;
hold on
plot(size,accuracy21,'bo-');
plot(size,accuracy22,'Color','#D95319','marker','x');
plot(size,accuracy23,'Color','#77AC30','marker','^');
xlabel('Number of steps')
ylabel('% Surfactant adsorbed')
legend('NPGE(2-chain)','TMPE(3-chain)','PEE(4-chain)')


% percent depend on the number of surfactant molecules
accuracy11=zeros(1,5);accuracy12=zeros(1,5);accuracy13=zeros(1,5);
accuracy21=zeros(1,5);accuracy22=zeros(1,5);accuracy23=zeros(1,5);
for iter=1:3
    size=zeros(1,5);
    ns=5;np=0.4;nt=100;np2=0.7*np; np3=0.5*np;
    acc1=zeros(1,5);acc2=zeros(1,5);acc3=zeros(1,5);
    ac1=zeros(1,5);ac2=zeros(1,5);ac3=zeros(1,5);
    for i=1:5
        count1=0; count2=0; count3=0;
        count11=0; count12=0; count13=0;
        nl=(i-1)*50;
        if i==1
            nl=nl+25;
        end
        for l=1:nl
            % set the initial point
            xy0=randi([0,100],1,2);
            if 20<xy0(1,1) && xy0(1,1)<80
                num=randi([0,100],1);
                xy0(1,1)=num;
            end
            if 20<xy0(1,2) && xy0(1,2)<80
                num=randi([0,100],1);
                xy0(1,2)=num;
            end
            xy=zeros(nt,2);
            xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
            % generate the random size of steps
            deltax=randi([-ns,ns],1,nt);
            deltay=randi([-ns,ns],1,nt);
            for step=1:nt-1
                % Walk in the x direction.
                if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
                    p=rand(1);
                    P(step)=p;
                    if np>p
                        xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                        break;
                    end
                end
                xy(step+1,1)=xy(step,1)+deltax(step);
                % apply periodic boundary condition
                if xy(step+1,1)>=0 && xy(step+1,1)<=nB
                    xy(step+1,1)=xy(step+1,1);
                end
                if xy(step+1,1)>nB
                    xy(step+1,1)=xy(step+1,1)-nB;
                end
                if xy(step+1,1)<0
                    xy(step+1,1)=xy(step+1,1)+nB;
                end
	            % Walk in the y direction.
	            xy(step+1,2)=xy(step,2)+deltay(step);
                if xy(step+1,2)>=0 && xy(step+1,2)<=nB
                    xy(step+1,2)=xy(step+1,2);
                end
                if xy(step+1,2)>nB
                    xy(step+1,2)=xy(step+1,2)-nB;
                end
                if xy(step+1,2)<0
                    xy(step+1,2)=xy(step+1,2)+nB;
                end
            end
            x(l)=xy(nt,1); y(l)=xy(nt,2);
            x_left(l)=x(l)-1;x_right(l)=x(l)+1;
            if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
                count11=count11+1;
            end
            if x_left(l)>20 && x_left(l)<80
                count1=count1+1;
            end
            if x_right(l)>20 && x_right(l)<80
                count1=count1+1;
            end
        end
        acc1(i)=count1*100/3600;
        ac1(i)=count11*100/nl;
        size(i)=nl;
            
        for l=1:nl
            % set the initial point
            xy0=randi([0,100],1,2);
            if 20<xy0(1,1) && xy0(1,1)<80
                num=randi([0,100],1);
                xy0(1,1)=num;
            end
            if 20<xy0(1,2) && xy0(1,2)<80
                num=randi([0,100],1);
                xy0(1,2)=num;
            end
            xy=zeros(nt,2);
            xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
            % generate the random size of steps
            deltax=randi([-ns,ns],1,nt);
            deltay=randi([-ns,ns],1,nt);
            for step=1:nt-1
                % Walk in the x direction.
                if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
                    p=rand(1);
                    P(step)=p;
                    if np2>p
                        xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                        break;
                    end
                end
                xy(step+1,1)=xy(step,1)+deltax(step);
                % apply periodic boundary condition
                if xy(step+1,1)>=0 && xy(step+1,1)<=nB
                    xy(step+1,1)=xy(step+1,1);
                end
                if xy(step+1,1)>nB
                    xy(step+1,1)=xy(step+1,1)-nB;
                end
                if xy(step+1,1)<0
                    xy(step+1,1)=xy(step+1,1)+nB;
                end
	            % Walk in the y direction.
	            xy(step+1,2)=xy(step,2)+deltay(step);
                if xy(step+1,2)>=0 && xy(step+1,2)<=nB
                    xy(step+1,2)=xy(step+1,2);
                end
                if xy(step+1,2)>nB
                    xy(step+1,2)=xy(step+1,2)-nB;
                end
                if xy(step+1,2)<0
                    xy(step+1,2)=xy(step+1,2)+nB;
                end
            end
            x(l)=xy(nt,1); y(l)=xy(nt,2);
            x_left(l)=x(l)-1;x_right(l)=x(l)+1;
            y_down(l)=y(l)-1;
            if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
                count12=count12+1;
            end
            if x_left(l)>20 && x_left(l)<80
                count2=count2+1;
            end
            if x_right(l)>20 && x_right(l)<80
                count2=count2+1;
            end
            if y_down(l)>20 && y_down(l)<80
                count2=count2+1;
            end
        end
        acc2(i)=count2*100/3600;
        ac2(i)=count12*100/nl;
            
        for l=1:nl
            % set the initial point
            xy0=randi([0,100],1,2);
            if 20<xy0(1,1) && xy0(1,1)<80
                num=randi([0,100],1);
                xy0(1,1)=num;
            end
            if 20<xy0(1,2) && xy0(1,2)<80
                num=randi([0,100],1);
                xy0(1,2)=num;
            end
            xy=zeros(nt,2);
            xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
            % generate the random size of steps
            deltax=randi([-ns,ns],1,nt);
            deltay=randi([-ns,ns],1,nt);
            for step=1:nt-1
                % Walk in the x direction.
                if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
                    p=rand(1);
                    P(step)=p;
                    if np3>p
                        xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                        break;
                    end
                end
                xy(step+1,1)=xy(step,1)+deltax(step);
                % apply periodic boundary condition
                if xy(step+1,1)>=0 && xy(step+1,1)<=nB
                    xy(step+1,1)=xy(step+1,1);
                end
                if xy(step+1,1)>nB
                    xy(step+1,1)=xy(step+1,1)-nB;
                end
                if xy(step+1,1)<0
                    xy(step+1,1)=xy(step+1,1)+nB;
                end
	            % Walk in the y direction.
	            xy(step+1,2)=xy(step,2)+deltay(step);
                if xy(step+1,2)>=0 && xy(step+1,2)<=nB
                    xy(step+1,2)=xy(step+1,2);
                end
                if xy(step+1,2)>nB
                    xy(step+1,2)=xy(step+1,2)-nB;
                end
                if xy(step+1,2)<0
                    xy(step+1,2)=xy(step+1,2)+nB;
                end
            end
            x(l)=xy(nt,1); y(l)=xy(nt,2);
            x_left(l)=x(l)-1;x_right(l)=x(l)+1;
            y_down(l)=y(l)-1;y_top(l)=y(l)+1;
            if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
                count13=count13+1;
            end
            if x_left(l)>20 && x_left(l)<80
                count3=count3+1;
            end
            if x_right(l)>20 && x_right(l)<80
                count3=count3+1;
            end
            if y_down(l)>20 && y_down(l)<80
                count3=count3+1;
            end
            if y_top(l)>20 && y_top(l)<80
                count3=count3+1;
            end
        end
        acc3(i)=count3*100/3600;
        ac3(i)=count13*100/nl;
    end
    accuracy11=accuracy11+acc1;
    accuracy21=accuracy21+ac1;
    accuracy12=accuracy12+acc2;
    accuracy22=accuracy22+ac2;
    accuracy13=accuracy13+acc3;
    accuracy23=accuracy23+ac3;
end
accuracy11=accuracy11/3;
accuracy21=accuracy21/3;
accuracy12=accuracy12/3;
accuracy22=accuracy22/3;
accuracy13=accuracy13/3;
accuracy23=accuracy23/3;
figure;
hold on
plot(size,accuracy11,'bo-');
plot(size,accuracy12,'Color','#D95319','marker','x');
plot(size,accuracy13,'Color','#77AC30','marker','^');
xlabel('Surfactant molecules')
ylabel('% Organoclay site occupied')
legend('NPGE(2-chain)','TMPE(3-chain)','PEE(4-chain)')
figure;
hold on
plot(size,accuracy21,'bo-');
plot(size,accuracy22,'Color','#D95319','marker','x');
plot(size,accuracy23,'Color','#77AC30','marker','^');
xlabel('Surfactant molecules')
ylabel('% Surfactant adsorbed')
legend('NPGE(2-chain)','TMPE(3-chain)','PEE(4-chain)')


% percent depend on various length of random walks(step size)
accuracy11=zeros(1,5);accuracy12=zeros(1,5);accuracy13=zeros(1,5);
accuracy21=zeros(1,5);accuracy22=zeros(1,5);accuracy23=zeros(1,5);
for iter=1:3
    size=zeros(1,5);
    nl=100;np=0.4;nt=100;np2=0.7*np; np3=0.5*np;
    acc1=zeros(1,5);acc2=zeros(1,5);acc3=zeros(1,5);
    ac1=zeros(1,5);ac2=zeros(1,5);ac3=zeros(1,5);
    for i=1:5
        count1=0; count2=0; count3=0;
        count11=0; count12=0; count13=0;
        ns=i;
        for l=1:nl
            % set the initial point
            xy0=randi([0,100],1,2);
            if 20<xy0(1,1) && xy0(1,1)<80
                num=randi([0,100],1);
                xy0(1,1)=num;
            end
            if 20<xy0(1,2) && xy0(1,2)<80
                num=randi([0,100],1);
                xy0(1,2)=num;
            end
            xy=zeros(nt,2);
            xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
            % generate the random size of steps
            deltax=randi([-ns,ns],1,nt);
            deltay=randi([-ns,ns],1,nt);
            for step=1:nt-1
                % Walk in the x direction.
                if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
                    p=rand(1);
                    P(step)=p;
                    if np>p
                        xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                        break;
                    end
                end
                xy(step+1,1)=xy(step,1)+deltax(step);
                % apply periodic boundary condition
                if xy(step+1,1)>=0 && xy(step+1,1)<=nB
                    xy(step+1,1)=xy(step+1,1);
                end
                if xy(step+1,1)>nB
                    xy(step+1,1)=xy(step+1,1)-nB;
                end
                if xy(step+1,1)<0
                    xy(step+1,1)=xy(step+1,1)+nB;
                end
	            % Walk in the y direction.
	            xy(step+1,2)=xy(step,2)+deltay(step);
                if xy(step+1,2)>=0 && xy(step+1,2)<=nB
                    xy(step+1,2)=xy(step+1,2);
                end
                if xy(step+1,2)>nB
                    xy(step+1,2)=xy(step+1,2)-nB;
                end
                if xy(step+1,2)<0
                    xy(step+1,2)=xy(step+1,2)+nB;
                end
            end
            x(l)=xy(nt,1); y(l)=xy(nt,2);
            x_left(l)=x(l)-1;x_right(l)=x(l)+1;
            if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
                count11=count11+1;
            end
            if x_left(l)>20 && x_left(l)<80
                count1=count1+1;
            end
            if x_right(l)>20 && x_right(l)<80
                count1=count1+1;
            end
        end
        acc1(i)=count1*100/3600;
        ac1(i)=count11*100/nl;
        size(i)=ns;
            
        for l=1:nl
            % set the initial point
            xy0=randi([0,100],1,2);
            if 20<xy0(1,1) && xy0(1,1)<80
                num=randi([0,100],1);
                xy0(1,1)=num;
            end
            if 20<xy0(1,2) && xy0(1,2)<80
                num=randi([0,100],1);
                xy0(1,2)=num;
            end
            xy=zeros(nt,2);
            xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
            % generate the random size of steps
            deltax=randi([-ns,ns],1,nt);
            deltay=randi([-ns,ns],1,nt);
            for step=1:nt-1
                % Walk in the x direction.
                if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
                    p=rand(1);
                    P(step)=p;
                    if np2>p
                        xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                        break;
                    end
                end
                xy(step+1,1)=xy(step,1)+deltax(step);
                % apply periodic boundary condition
                if xy(step+1,1)>=0 && xy(step+1,1)<=nB
                    xy(step+1,1)=xy(step+1,1);
                end
                if xy(step+1,1)>nB
                    xy(step+1,1)=xy(step+1,1)-nB;
                end
                if xy(step+1,1)<0
                    xy(step+1,1)=xy(step+1,1)+nB;
                end
	            % Walk in the y direction.
	            xy(step+1,2)=xy(step,2)+deltay(step);
                if xy(step+1,2)>=0 && xy(step+1,2)<=nB
                    xy(step+1,2)=xy(step+1,2);
                end
                if xy(step+1,2)>nB
                    xy(step+1,2)=xy(step+1,2)-nB;
                end
                if xy(step+1,2)<0
                    xy(step+1,2)=xy(step+1,2)+nB;
                end
            end
            x(l)=xy(nt,1); y(l)=xy(nt,2);
            x_left(l)=x(l)-1;x_right(l)=x(l)+1;
            y_down(l)=y(l)-1;
            if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
                count12=count12+1;
            end
            if x_left(l)>20 && x_left(l)<80
                count2=count2+1;
            end
            if x_right(l)>20 && x_right(l)<80
                count2=count2+1;
            end
            if y_down(l)>20 && y_down(l)<80
                count2=count2+1;
            end
        end
        acc2(i)=count2*100/3600;
        ac2(i)=count12*100/nl;
            
        for l=1:nl
            % set the initial point
            xy0=randi([0,100],1,2);
            if 20<xy0(1,1) && xy0(1,1)<80
                num=randi([0,100],1);
                xy0(1,1)=num;
            end
            if 20<xy0(1,2) && xy0(1,2)<80
                num=randi([0,100],1);
                xy0(1,2)=num;
            end
            xy=zeros(nt,2);
            xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
            % generate the random size of steps
            deltax=randi([-ns,ns],1,nt);
            deltay=randi([-ns,ns],1,nt);
            for step=1:nt-1
                % Walk in the x direction.
                if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
                    p=rand(1);
                    P(step)=p;
                    if np3>p
                        xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                        break;
                    end
                end
                xy(step+1,1)=xy(step,1)+deltax(step);
                % apply periodic boundary condition
                if xy(step+1,1)>=0 && xy(step+1,1)<=nB
                    xy(step+1,1)=xy(step+1,1);
                end
                if xy(step+1,1)>nB
                    xy(step+1,1)=xy(step+1,1)-nB;
                end
                if xy(step+1,1)<0
                    xy(step+1,1)=xy(step+1,1)+nB;
                end
	            % Walk in the y direction.
	            xy(step+1,2)=xy(step,2)+deltay(step);
                if xy(step+1,2)>=0 && xy(step+1,2)<=nB
                    xy(step+1,2)=xy(step+1,2);
                end
                if xy(step+1,2)>nB
                    xy(step+1,2)=xy(step+1,2)-nB;
                end
                if xy(step+1,2)<0
                    xy(step+1,2)=xy(step+1,2)+nB;
                end
            end
            x(l)=xy(nt,1); y(l)=xy(nt,2);
            x_left(l)=x(l)-1;x_right(l)=x(l)+1;
            y_down(l)=y(l)-1;y_top(l)=y(l)+1;
            if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
                count13=count13+1;
            end
            if x_left(l)>20 && x_left(l)<80
                count3=count3+1;
            end
            if x_right(l)>20 && x_right(l)<80
                count3=count3+1;
            end
            if y_down(l)>20 && y_down(l)<80
                count3=count3+1;
            end
            if y_top(l)>20 && y_top(l)<80
                count3=count3+1;
            end
        end
        acc3(i)=count3*100/3600;
        ac3(i)=count13*100/nl;
    end
    accuracy11=accuracy11+acc1;
    accuracy21=accuracy21+ac1;
    accuracy12=accuracy12+acc2;
    accuracy22=accuracy22+ac2;
    accuracy13=accuracy13+acc3;
    accuracy23=accuracy23+ac3;
end
accuracy11=accuracy11/3;
accuracy21=accuracy21/3;
accuracy12=accuracy12/3;
accuracy22=accuracy22/3;
accuracy13=accuracy13/3;
accuracy23=accuracy23/3;
figure;
hold on
plot(size,accuracy11,'bo-');
plot(size,accuracy12,'Color','#D95319','marker','x');
plot(size,accuracy13,'Color','#77AC30','marker','^');
xlabel('Step size(Lattice units)')
ylabel('% Organoclay site occupied')
legend('NPGE(2-chain)','TMPE(3-chain)','PEE(4-chain)')
figure;
hold on
plot(size,accuracy21,'bo-');
plot(size,accuracy22,'Color','#D95319','marker','x');
plot(size,accuracy23,'Color','#77AC30','marker','^');
xlabel('Step size(Lattice units)')
ylabel('% Surfactant adsorbed')
legend('NPGE(2-chain)','TMPE(3-chain)','PEE(4-chain)')


% percent depend on the adsorption probabilty
accuracy11=zeros(1,5);accuracy12=zeros(1,5);accuracy13=zeros(1,5);
accuracy21=zeros(1,5);accuracy22=zeros(1,5);accuracy23=zeros(1,5);
for iter=1:10
    size=zeros(1,5);
    ns=5;nl=100;nt=100;
    acc1=zeros(1,5);acc2=zeros(1,5);acc3=zeros(1,5);
    ac1=zeros(1,5);ac2=zeros(1,5);ac3=zeros(1,5);
    for i=1:5
        count1=0; count2=0; count3=0;
        count11=0; count12=0; count13=0;
        np=(i+3)*0.1;np2=0.5*np; np3=0.3*np;
        for l=1:nl
            % set the initial point
            xy0=randi([0,100],1,2);
            if 20<xy0(1,1) && xy0(1,1)<80
                num=randi([0,100],1);
                xy0(1,1)=num;
            end
            if 20<xy0(1,2) && xy0(1,2)<80
                num=randi([0,100],1);
                xy0(1,2)=num;
            end
            xy=zeros(nt,2);
            xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
            % generate the random size of steps
            deltax=randi([-ns,ns],1,nt);
            deltay=randi([-ns,ns],1,nt);
            for step=1:nt-1
                % Walk in the x direction.
                if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
                    p=rand(1);
                    P(step)=p;
                    if np>p
                        xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                        break;
                    end
                end
                xy(step+1,1)=xy(step,1)+deltax(step);
                % apply periodic boundary condition
                if xy(step+1,1)>=0 && xy(step+1,1)<=nB
                    xy(step+1,1)=xy(step+1,1);
                end
                if xy(step+1,1)>nB
                    xy(step+1,1)=xy(step+1,1)-nB;
                end
                if xy(step+1,1)<0
                    xy(step+1,1)=xy(step+1,1)+nB;
                end
	            % Walk in the y direction.
	            xy(step+1,2)=xy(step,2)+deltay(step);
                if xy(step+1,2)>=0 && xy(step+1,2)<=nB
                    xy(step+1,2)=xy(step+1,2);
                end
                if xy(step+1,2)>nB
                    xy(step+1,2)=xy(step+1,2)-nB;
                end
                if xy(step+1,2)<0
                    xy(step+1,2)=xy(step+1,2)+nB;
                end
            end
            x(l)=xy(nt,1); y(l)=xy(nt,2);
            x_left(l)=x(l)-1;x_right(l)=x(l)+1;
            if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
                count11=count11+1;
            end
            if x_left(l)>20 && x_left(l)<80
                count1=count1+1;
            end
            if x_right(l)>20 && x_right(l)<80
                count1=count1+1;
            end
        end
        ac1(i)=count11*100/nl;
        acc1(i)=count11*200/3600;
        size(i)=np-0.2;
            
        for l=1:nl
            % set the initial point
            xy0=randi([0,100],1,2);
            if 20<xy0(1,1) && xy0(1,1)<80
                num=randi([0,100],1);
                xy0(1,1)=num;
            end
            if 20<xy0(1,2) && xy0(1,2)<80
                num=randi([0,100],1);
                xy0(1,2)=num;
            end
            xy=zeros(nt,2);
            xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
            % generate the random size of steps
            deltax=randi([-ns,ns],1,nt);
            deltay=randi([-ns,ns],1,nt);
            for step=1:nt-1
                % Walk in the x direction.
                if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
                    p=rand(1);
                    P(step)=p;
                    if np2>p
                        xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                        break;
                    end
                end
                xy(step+1,1)=xy(step,1)+deltax(step);
                % apply periodic boundary condition
                if xy(step+1,1)>=0 && xy(step+1,1)<=nB
                    xy(step+1,1)=xy(step+1,1);
                end
                if xy(step+1,1)>nB
                    xy(step+1,1)=xy(step+1,1)-nB;
                end
                if xy(step+1,1)<0
                    xy(step+1,1)=xy(step+1,1)+nB;
                end
	            % Walk in the y direction.
	            xy(step+1,2)=xy(step,2)+deltay(step);
                if xy(step+1,2)>=0 && xy(step+1,2)<=nB
                    xy(step+1,2)=xy(step+1,2);
                end
                if xy(step+1,2)>nB
                    xy(step+1,2)=xy(step+1,2)-nB;
                end
                if xy(step+1,2)<0
                    xy(step+1,2)=xy(step+1,2)+nB;
                end
            end
            x(l)=xy(nt,1); y(l)=xy(nt,2);
            x_left(l)=x(l)-1;x_right(l)=x(l)+1;
            y_down(l)=y(l)-1;
            if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
                count12=count12+1;
            end
            if x_left(l)>20 && x_left(l)<80
                count2=count2+1;
            end
            if x_right(l)>20 && x_right(l)<80
                count2=count2+1;
            end
            if y_down(l)>20 && y_down(l)<80
                count2=count2+1;
            end
        end
        ac2(i)=count12*100/nl;
        acc2(i)=count12*300/3600;
            
        for l=1:nl
            % set the initial point
            xy0=randi([0,100],1,2);
            if 20<xy0(1,1) && xy0(1,1)<80
                num=randi([0,100],1);
                xy0(1,1)=num;
            end
            if 20<xy0(1,2) && xy0(1,2)<80
                num=randi([0,100],1);
                xy0(1,2)=num;
            end
            xy=zeros(nt,2);
            xy(1,1)=xy0(1,1); xy(1,2)=xy0(1,2);
            % generate the random size of steps
            deltax=randi([-ns,ns],1,nt);
            deltay=randi([-ns,ns],1,nt);
            for step=1:nt-1
                % Walk in the x direction.
                if xy(step,1)>20 && xy(step,1)<80 && xy(step,2)>20 && xy(step,2)<80
                    p=rand(1);
                    P(step)=p;
                    if np3>p
                        xy(nt,1)=xy(step,1); xy(nt,2)=xy(step,2);
                        break;
                    end
                end
                xy(step+1,1)=xy(step,1)+deltax(step);
                % apply periodic boundary condition
                if xy(step+1,1)>=0 && xy(step+1,1)<=nB
                    xy(step+1,1)=xy(step+1,1);
                end
                if xy(step+1,1)>nB
                    xy(step+1,1)=xy(step+1,1)-nB;
                end
                if xy(step+1,1)<0
                    xy(step+1,1)=xy(step+1,1)+nB;
                end
	            % Walk in the y direction.
	            xy(step+1,2)=xy(step,2)+deltay(step);
                if xy(step+1,2)>=0 && xy(step+1,2)<=nB
                    xy(step+1,2)=xy(step+1,2);
                end
                if xy(step+1,2)>nB
                    xy(step+1,2)=xy(step+1,2)-nB;
                end
                if xy(step+1,2)<0
                    xy(step+1,2)=xy(step+1,2)+nB;
                end
            end
            x(l)=xy(nt,1); y(l)=xy(nt,2);
            x_left(l)=x(l)-1;x_right(l)=x(l)+1;
            y_down(l)=y(l)-1;y_top(l)=y(l)+1;
            if x(l)>20 && x(l)<80 && y(l)>20 && y(l)<80
                count13=count13+1;
            end
            if x_left(l)>20 && x_left(l)<80
                count3=count3+1;
            end
            if x_right(l)>20 && x_right(l)<80
                count3=count3+1;
            end
            if y_down(l)>20 && y_down(l)<80
                count3=count3+1;
            end
            if y_top(l)>20 && y_top(l)<80
                count3=count3+1;
            end
        end
        ac3(i)=count13*100/nl;
        acc3(i)=count13*400/3600;
    end
    accuracy11=accuracy11+acc1;
    accuracy21=accuracy21+ac1;
    accuracy12=accuracy12+acc2;
    accuracy22=accuracy22+ac2;
    accuracy13=accuracy13+acc3;
    accuracy23=accuracy23+ac3;
end
accuracy11=accuracy11/10;
accuracy21=accuracy21/10;%accuracy11=accuracy21/18;
accuracy12=accuracy12/10;
accuracy22=accuracy22/10;%accuracy12=accuracy22/12;
accuracy13=accuracy13/10;
accuracy23=accuracy23/10;%accuracy13=accuracy23/9;
figure;
hold on
plot(size,accuracy11,'bo-');
plot(size,accuracy12,'Color','#D95319','marker','x');
plot(size,accuracy13,'Color','#77AC30','marker','^');
xlabel('Adsorption probability')
ylabel('% Organoclay site occupied')
legend('NPGE(2-chain)','TMPE(3-chain)','PEE(4-chain)')
figure;
hold on
plot(size,accuracy21,'bo-');
plot(size,accuracy22,'Color','#D95319','marker','x');
plot(size,accuracy23,'Color','#77AC30','marker','^');
xlabel('Adsorption probability')
ylabel('% Surfactant adsorbed')
legend('NPGE(2-chain)','TMPE(3-chain)','PEE(4-chain)')