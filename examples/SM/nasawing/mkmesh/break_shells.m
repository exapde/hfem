function [ q2sh,sh2th,e2b,ib ] = break_shells(e2q,q2e,q2sh,sh2th)
%BREAK_SHELLS Summary of this function goes here
%   Detailed explanation goes here

% Number of elements (quads), edges, and shells
ne  = size(q2sh,1);
ned = size(e2q,1);
nsh = size(sh2th,1);

% Init vector of boolean (true if edge on line)
eOnLine = false([ned,1]);
% Init vector of boolean (true if edge on Cutting line)
eOnCutL = false([ned,1]);

% Mark edges on intersections between different shells
for i=1:ned
    nearSh = q2sh(e2q{i});
    if any(diff(nearSh))
        eOnLine(i) = true;
        if ~all(diff(sort(nearSh)))
            eOnCutL(i) = true;
        end
    end
end

% Initialized reorganized quads vector
reorgQ = cell(1);
doneQuad = [];
% New shell index
shInd = 0;

for i =1:ned
    if eOnCutL(i)
        % If Edge on cutting line, get separated quads
        [q12,shBroken] = qOnSameShell(e2q{i},q2sh);
        for j=1:2
            qj = q12(j);
            if ~any(doneQuad==qj)
                % Get all the quads in the same subshell
                qSet = findQinSameShell(e2q,q2e,eOnLine,qj);
                % Store all the quads in current subshell
                shInd = shInd+1;
                reorgQ{shInd}.Quads = qSet;
                % Save Id of the broken Parent Shell
                reorgQ{shInd}.ParSh = shBroken;
                % Store all treated quads
                doneQuad = [doneQuad ; qSet];
            end
        end
    end
end

% Assign the new Shells Id to quads previously in broken shell
% And update vector sh2th for new shells
for i=1:length(reorgQ)
    shInd = nsh + i;
    q2sh(reorgQ{i}.Quads) = shInd;
    % Thickness of the parent shell
    newTh = sh2th(reorgQ{i}.ParSh);
    sh2th = [sh2th; newTh];
end

% Building list of shell intersections (i.e "beams") ib,
% and the edge to beam (e2b) connectivity
e2b = zeros(ned,1);
ib  = {};
for i =1:ned
    if eOnLine(i)
        % Get shells connected to current edge
        nearSh = sort(q2sh(e2q{i}));
        if (~all(diff(nearSh))) error('Some Shells have not been broken'), end
        % Check if connection set hasnt been saved yet
        isNewInter = true;
        for j=1:length(ib)
            if(isequal(ib{j},nearSh))
                isNewInter = false;
                e2b(i) = j;
                break;
            end
        end
        % Enrich ib with a new connexion set
        if isNewInter
            j = length(ib)+1;
            ib{j} = nearSh;
            e2b(i) = j;
        end
    end
end

end


function [ q12,shell ] = qOnSameShell(e2q,q2sh)
% Gets the quads q1 and q2 that are on the same shell and that needs to be
% separated. Gets also the shell ID where these quads are.

% Shells connected to current edge
shells = q2sh(e2q);

% Get the first repeated shell and corresp quads
for i=1:length(e2q)
    f = find(shells==shells(i));
    if length(f)==2
        shell = shells(i);
        q12 = e2q(f);
    end
end
end


function [ qSet ] = findQinSameShell(e2q,q2e,eOnLine,q)
% Gets the set of quads qSet that are on the same shell than q

qSet = [];
qToDo = [q];
NotOnLine = ~eOnLine;

while (~isempty(qToDo))
    qi = qToDo(end);
    % Get neighbouring edges
    edges = q2e(qi,:);
    % remove edges on an intersection line
    edges = edges' .* NotOnLine(edges);
    edges(edges==0) = [];
    % Get neighbouring quads
    nquads = [];
    for i=1:length(edges)
        nquads = [nquads; e2q{edges(i)}];
    end
    in_qSet  = ~ismember(nquads,qSet);
    in_qToDo = ~ismember(nquads,qToDo);
    % Remove qi from qToDo
    qToDo(end) = [];
    % Update qToDo
    qnew  = in_qSet.*in_qToDo.*nquads;
    qnew(qnew==0) = [];
    qToDo = [qToDo ; qnew];
    % Update qSet
    qSet  = [qSet ; qi];
end
end

