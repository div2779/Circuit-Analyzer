syms s t f(t);
%Reading Input from File
line_spec = "%d %d %d %d %d %s %s";
sizeA = [Inf 7];
circuit_file=fopen('Input File.txt', 'r');
A=textscan(circuit_file, line_spec, sizeA);
%A=A';
fclose(circuit_file);

vg=[A{6}];  %A{6} represents cell column of input voltages
ig=[A{7}];  %A{7} represents cell column of input currents
num_branches=size(vg,1);

numOfNodes = findNumOfNodes(A);
[parent,rank] = makeSet(numOfNodes);
[parent, rank, tree] = makeUnions(A, parent, rank, num_branches);
% tree = [1 3 6 10];

%we have been given tree branches. links are given by:
x=1:1:num_branches;
links = setdiff(x, tree); 

%The input read from the file was cell column of strings
%Converting strings to symbolic expression with symbolic
%variable s
for i = 1:num_branches
    vg{i}=eval(str2sym(vg{i}));
    ig{i}=eval(str2sym(ig{i}));
end

%Let x denote the sequence of twig branches followed
%by link branches
x=[tree links];

%The sybolic expression input is in the form of cells.
%Converting cell array to normal symbolic array
Vg=[];
for i=1:num_branches
    Vg=[Vg;vg{i}];
end
Vgt=Vg(tree);
Vg=Vg(x); %Rearranging elements of Vg as per sequence in x

Ig=[];
for i=1:num_branches
    Ig=[Ig;ig{i}/1000];  %Input will be in mA
end
Igl=Ig(links);
Ig=Ig(x); %Rearranging elements of Ig as per sequence in x


keySet = 1:1:numOfNodes;
valueSet = zeros(1,numOfNodes);
adjList = containers.Map(keySet, valueSet, 'UniformValues', false);
% fill map like ; src -> [dst direction branchNum]'

for i = tree
    if adjList(A{1}(i)) == 0
        adjList(A{1}(i)) = [A{2}(i) 1 i]';      %[from dirn to]'
    else
        adjList(A{1}(i)) = [adjList(A{1}(i)) [A{2}(i) 1 i]'];
    end
    if adjList(A{2}(i)) == 0
        adjList(A{2}(i)) = [A{1}(i) -1 i]';
    else
        adjList(A{2}(i)) = [adjList(A{2}(i)) [A{1}(i) -1 i]'];
    end
end

% Fundamental Loop Matrix
Bf = zeros(num_branches, num_branches);
for link = links
    Bf(link, link)=1;
    Bf = fillBft(Bf, A{1}(link), link, A{2}(link), link, adjList);
end
Bf = Bf(links,[tree links])
Bft = Bf(:,1:1:numel(tree));

% Fundamental Cutset Matrix
Qfl=-Bft';
Qft=eye(numel(tree));
Qf=[Qft Qfl];

%Extracting & expressing the inductive reactance in array 'l'
%where each input inductance is assumed to be in mH
l=[];
for i=1:num_branches
    l=[l A{3}(i)];
end
l=double(l);
l=l.*0.001.*s;

%Extracting & expressing the resistance in array 'r'
%where each input resistance is assumed to be in kOhm
r=[];
for i=1:num_branches
    r=[r A{4}(i)];
end
r=double(r);
r=r.*1000;

%Extracting & expressing the capacitive reactance in array 'c'
%where each input capacitance is assumed to be in uF
c=[];
for i=1:num_branches
    temp=A{5}(i);
    if temp
        c=[c 1/temp];
    else
        c=[c 0];
    end
end
c=double(c);
c=(c.*0.000001)./s;

%Net impedance in Laplace domain
z=l+c+r;
z=z(x);

%For increasing efficiency, we need as small matrix as possible for 
%inverting.Thus we choose Z path or Y path depending on which is lower,
% # of links or # of twigs.
if(numel(links)<=numel(tree))
    Z=diag(z);
    Zf=Bf*Z*(Bf');
    Zf_inv=inv(Zf);
    Il=Zf_inv*Bf*(Vg+(Z*(Bf')*Igl)-Z*Ig);
    I=(Bf')*Il+Ig-(Bf')*Igl;
    I=I.*1000; %output current is in mA
    V=Z*I;
else
    y=z.^(-1);
    Y=diag(y);
    Yf=Qf*Y*(Qf');
    Yf_inv=inv(Yf);
    Vt=Yf_inv*Qf*(Ig+(Y*(Qf')*Vgt)-Y*Vg);
    V=(Qf')*Vt+Vg-(Qf')*Vgt;
    I=Y*V;
    I=I.*1000; %output current is in mA
end
 
% for i=1:num_branches
%     V_new(i)=partfrac(V(i));
%     V_new(i)=vpa(V_new(i));
%     V_t(i)=ilaplace(V_new(i));
%     V_t(i)=vpa(V_t(i));
%     I_new(i)=partfrac(I(i));
%     I_t(i)=ilaplace(I_new(i));
%     I_t(i)=vpa(I_t(i));
% end


outputVoltageBranch = 10;   %HardCoded
outputCurrentBranch = 7;   %HardCoded
posV = find(x==outputVoltageBranch);
posI = find(x==outputCurrentBranch);



vPos(t) = ilaplace(V(posV));
vPos(t) = simplify(vPos);
vPos(t) = vpa(vPos);
figure(1)
fplot(vPos);
xlim([0 0.006]);
hold off;

iPos(t) = ilaplace(I(posI));
iPos(t) = simplify(iPos);
iPos(t) = vpa(iPos);
figure(2)
fplot(iPos);
xlim([0 0.006]);

function numOfNodes = findNumOfNodes(A)
    numOfNodes=max([A{1}' A{2}']);
end


function [parent, rank] = makeSet(numOfNodes)
    parent = 1:1:numOfNodes;
    rank = ones(1,numOfNodes);
end


function [parentOfA, parent] = findParent(parent, a)
    if parent(a)==a
        parentOfA = a;
        return
    end
    
    parent(a) = findParent(parent,parent(a));
    parentOfA = parent(a);
end


function [parent,rank] = union(a, b, parent, rank)
    a = parent(a);
    b = parent(b);
    
    if a~=b
        if rank(a)>rank(b)
            parent(b)=a;
            rank(a) = rank(a) + rank(b);
        else
            parent(a)=b;
            rank(b) = rank(b) + rank(a);
        end
    end 
end

function [parent, rank, tree] = makeUnions(A, parent, rank, num_branches)
    tree=[];
    for i=1:1:num_branches
        src=A{1}(i);
        dst=A{2}(i);
        [parentOfA, parent] = findParent(parent, src);
        [parentOfB, parent] = findParent(parent, dst);
        if(parentOfA~=parentOfB)
            tree = [tree i];
            [parent, rank] = union(src, dst, parent, rank);
        end
    end    
end

%Recursive Function for finding the Bf matrix.
function [Bf, flag] = fillBft(Bf, linkSrc, prevBranch, currNode, link, adjList)
    flag=0;
    if(currNode==linkSrc)
        flag=1;
        return;
    end
    
    for neighbourDetails = adjList(currNode)
        neighbour = neighbourDetails(1,1);
        direction = neighbourDetails(2,1);
        branchNum = neighbourDetails(3,1);
%      if(neighbour~=prevNode)   -> will give wrong result in parallel loop
        if(prevBranch~=branchNum)
            [Bf, flagNext] = fillBft(Bf, linkSrc, branchNum, neighbour, link, adjList);
            if(flagNext==1)
                Bf(link, branchNum) = direction;
                flag=1;
            end
        end
    end
end
