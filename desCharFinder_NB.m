% Design Characterization Script
% This script examines all fully-symmetric, flip-symmetric, and
% asymmetric designs within a design space for the presence of key
% characteristics associated with symmetry
function charCounts = desCharFinder_NB(CA,NC,sel,sidenum,desiredChars)
    % Initialize values
    charCounts = zeros(length(desiredChars),1);
    
    % Populate boolean values
    for q = 1:1:length(desiredChars)
        if desiredChars(q) == 1
            %%% Boolean Characteristic 1: Presence of Long Diagonals
            charCounts(q) = longDiagChecker(CA,NC,sel,sidenum);
        elseif desiredChars(q) == 2
            %%% Boolean Characteristic 2: Presence of Pivot Point(s)
            charCounts(q) = pivotPointChecker(NC,CA,sidenum);
        elseif desiredChars(q) == 4
            %%% Boolean Characteristic 4: Presence of Arrows
            charCounts(q) = arrowCheck(sidenum,NC,CA);
        elseif desiredChars(q) == 5
            %%% Boolean Characteristic 5: Presence of Spider Nodes
            charCounts(q) = spiderNodeChecker(sidenum,CA,NC,sel);
        elseif desiredChars(q) == 6
            %%% Boolean Characteristic 6: Presence of Stacked Members
            charCounts(q) = stackedMembersChecker(NC,CA);
        end
    end
end

% FUNCTION FOR CHARACTERISTIC #1 (LONG DIAGONALS)
function ldCount = longDiagChecker(CA,NC,sel,sidenum)
    % Initialize boolean
    ldCount = 0;

    % Loop through each member in the current design
    for i = 1:size(CA,1)
        % Finding member length, angle from nodal coordinates
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        angle = abs(acosd((x2-x1)./L));
        
        shortest45Diag = sqrt(2)*(sel/(sidenum-1));
        if (angle ~= 0) && (angle ~= 90) && (angle ~= 180)
            if abs(L) > shortest45Diag
                ldCount = ldCount + 1;
            end
        end
    end
end

% FUNCTION FOR CHARACTERISTIC #2 (PIVOT POINTS)
function ppCount = pivotPointChecker(NC,CA,sidenum)
    % Initialize Boolean
    ppCount = 0;

    % Loop through each horizontal row of nodes
    [N,~] = histcounts(CA,size(NC,1));
    for j = 1:1:sidenum
        hrow = j:sidenum:((sidenum^2)-sidenum+j);
        connRow = N(hrow); numPivots = 0;
        % For each node in the row:
        for q = 1:1:length(connRow)
            % Isolate all elements starting/ending at each node
            id = connRow(q);
            indione = CA(:,1) == id; inditwo = CA(:,2) == id;
            mCAone = CA(indione,:); mCAtwo = CA(inditwo,:);
            mCA = ...
             [setdiff(mCAone,mCAtwo,'rows');setdiff(mCAtwo,mCAone,'rows')];
            % Identify angles of all elements, relative to the current node
            x1 = NC(mCA(:,1),1); x2 = NC(mCA(:,2),1);
            y1 = NC(mCA(:,1),2); y2 = NC(mCA(:,2),2);
            L = sqrt(((x2-x1).^2)+((y2-y1).^2));
            angles = acos((x2-x1)./L);
            % Determine how many members have angles [0,pi] and [-pi,0]
            posangs = 0; negangs = 0;
            for o = 1:1:length(angles)
                if y1(o) > y2(o)
                    angles(o) = angles(o) - pi;
                end
                if (angles(o) < 0) && (angles(o) > -pi)
                    negangs = negangs + 1;
                elseif (angles(o) > 0) && (angles(o) < pi)
                    posangs = posangs + 1;
                end
            end
            % Tally up potential pivot point nodes in row
            if (posangs ~= 0) && (negangs ~= 0)
                numPivots = numPivots + 1;
            end
        end
        % Add up all pivot points, if any, in current row
        ppCount = ppCount + numPivots;
    end
    % Loop through each vertical column of nodes
    for k = 1:sidenum:((sidenum^2)-sidenum+1)
        vcol = k:1:(k+sidenum-1);
        connRow = N(vcol); numPivots = 0;
        % For each node in the col:
        for q = 1:1:length(connRow)
            % Isolate all elements starting/ending at each node
            id = connRow(q);
            indione = CA(:,1) == id; inditwo = CA(:,2) == id;
            mCAone = CA(indione,:); mCAtwo = CA(inditwo,:);
            mCA = ...
             [setdiff(mCAone,mCAtwo,'rows');setdiff(mCAtwo,mCAone,'rows')];
            % Identify angles of all elements, relative to the current node
            x1 = NC(mCA(:,1),1); x2 = NC(mCA(:,2),1);
            y1 = NC(mCA(:,1),2); y2 = NC(mCA(:,2),2);
            L = sqrt(((x2-x1).^2)+((y2-y1).^2));
            angles = acos((x2-x1)./L);
            % Determine how many members have angles [0,pi] and [-pi,0]
            posangs = 0; negangs = 0;
            for o = 1:1:length(angles)
                if y1(o) > y2(o)
                    angles(o) = angles(o) - pi;
                end
                if (angles(o) < 0) && (angles(o) > -pi)
                    negangs = negangs + 1;
                elseif (angles(o) > 0) && (angles(o) < pi)
                    posangs = posangs + 1;
                end
            end
            % Tally up potential pivot point nodes in row
            if (posangs ~= 0) && (negangs ~= 0)
                numPivots = numPivots + 1;
            end
        end
        % Add up all pivot points, if any, in current row
        ppCount = ppCount + numPivots;
    end
end

% FUNCTION FOR CHARACTERISTIC #4 (ARROWS)
function arrowCount = arrowCheck(sidenum,NC,CA)
    % Initialize values
    arrowCount = 0;
    allnodes = 1:1:(sidenum^2);
    allquartets = nchoosek(allnodes,4);
    SortedCA = sortrows((sort(CA'))');
    
    % Loop through all possible combinations of 4 nodes
    for i = 1:1:size(allquartets,1)
        quartet = allquartets(i,:);
        allorders = perms(quartet);
        % Loop through each ordered series of 4 nodes for each combination
        for j = 1:1:size(allorders,1)
            % Determine if 4 line segments (i.e. members)exist between the 
            % four points, in order
            fn = allorders(j,:);
            mCA = [fn(1),fn(2);fn(2),fn(3);fn(3),fn(4);fn(4),fn(1)];
            smCA = sortrows((sort(mCA'))');
            if all(ismember(smCA,SortedCA,'rows'))
                % Determine interior angles between all 4 line segments
                m1 = (NC(mCA(1,2),2)-NC(mCA(1,1),2))/...
                     (NC(mCA(1,2),1)-NC(mCA(1,1),1));
                m1 = round(m1,3);
                m2 = (NC(mCA(2,2),2)-NC(mCA(2,1),2))/...
                     (NC(mCA(2,2),1)-NC(mCA(2,1),1));
                m2 = round(m2,3);
                m3 = (NC(mCA(3,2),2)-NC(mCA(3,1),2))/...
                     (NC(mCA(3,2),1)-NC(mCA(3,1),1));
                m3 = round(m3,3);
                m4 = (NC(mCA(4,2),2)-NC(mCA(4,1),2))/...
                     (NC(mCA(4,2),1)-NC(mCA(4,1),1));
                m4 = round(m4,4);
                slopes = [m1,m2,m3,m4,m1];
                angles = [];
                for k = 1:1:length(slopes)-1
                    ma = slopes(k); mb = slopes(k+1);
                    if ma*mb == -1
                        angles = [angles,(pi/2)];
                    elseif ma == mb
                        angles = [angles,pi];
                    elseif isinf(ma)
                        mc = 0;
                        temptheta = (atan(abs((mc-mb)/(1+(mc*mb)))));
                        angles = [angles,((pi/2)-temptheta)];
                    elseif isinf(mb)
                        mc = 0;
                        temptheta = (atan(abs((ma-mc)/(1+(ma*mc)))));
                        angles = [angles,((pi/2)-temptheta)];
                    else
                        angles = [angles,(atan(abs((ma-mb)/(1+(ma*mb)))))];
                    end
                end
                revangles = abs(angles-(2*pi));
                
                % Determine if there are 3 interior acute angles and one
                % interior obtuse angle
                for q = 1:1:length(angles)
                    angleset = angles; angles(q) = revangles(q);
                    acutes = angleset(angleset > 0 & angleset < (pi/2));
                    obtuse = angleset(angleset > pi);
                    if (length(acutes) == 3) && (length(obtuse) == 1)
                        arrowCount = arrowCount + 1;
                    end
                end
            end
        end
    end
end

% FUNCTION FOR CHARACTERISTIC #5 (SPIDER NODES)
function spiderCount = spiderNodeChecker(sidenum,CA,NC,sel)
    % Initialize values
    spiderCount = 0;
    
    % Use connectivity counter to find conenctivity
    [N,~] = connectivityCounter(sidenum,CA,NC,sel);
    
    % Determine if any nodes have greater than 5 connections
    for i = 1:1:length(N)
        if N(i) > 5
            % Find all elements connecting to/from the current node
            indone = CA(:,1) == i; indtwo = CA(:,2) == i;
            mCAone = CA(indone,:); mCAtwo = CA(indtwo,:);
            mCA = [setdiff(mCAone,mCAtwo,'rows');...
                    setdiff(mCAtwo,mCAone,'rows')];
            revidx = mCA(:,1) ~= i;
            mCA(revidx,:) = [mCA(revidx,2),mCA(revidx,1)];
            
            % Determine if any nodes have at least two angles between
            % members that are < 45 degrees
            slopes = [];
            for j = 1:1:(size(mCA,1)-1)
                m = (NC(mCA(j,2),2)-NC(mCA(j,1),2))/...
                     (NC(mCA(j,2),1)-NC(mCA(j,1),1));
                slopes = [slopes,round(m,3)];
            end
            if ~isempty(slopes)
                slopes = [slopes,slopes(1)];
                angles = [];
                for k = 1:1:length(slopes)-1
                    ma = slopes(k); mb = slopes(k+1);
                    if ma*mb == -1
                        angles = [angles,(pi/2)];
                    elseif ma == mb
                        angles = [angles,pi];
                    elseif isinf(ma)
                        mc = 0;
                        temptheta = (atan(abs((mc-mb)/(1+(mc*mb)))));
                        angles = [angles,((pi/2)-temptheta)];
                    elseif isinf(mb)
                        mc = 0;
                        temptheta = (atan(abs((ma-mc)/(1+(ma*mc)))));
                        angles = [angles,((pi/2)-temptheta)];
                    else
                        angles = [angles,(atan(abs((ma-mb)/(1+(ma*mb)))))];
                    end
                end
                if any(angles < (pi/4))
                    spiderCount = spiderCount + 1;
                end
            end
        end
    end
end

% FUNCTION FOR CHARACTERISTIC #6 (STACKED MEMBERS)
function stackedCount = stackedMembersChecker(NC,CA)
    % definition of stacked member: 3 or more members with orientations
    % deviating no more than X degrees from each other, which all
    % share no more start or end points than the number of members present
    
    % Initialize values
    stackedCount = 0;
    
    % Define deviation angle X;
    angX = deg2rad(5);
    
    % Sort CA
    SortedCA = sortrows((sort(CA'))');
    
    % Find angles of all members in CA
    x1 = NC(SortedCA(:,1),1); x2 = NC(SortedCA(:,2),1);
    y1 = NC(SortedCA(:,1),2); y2 = NC(SortedCA(:,2),2);
    L = sqrt(((x2-x1).^2)+((y2-y1).^2));
    %angles = atan((y2-y1)./(x2-x1));
    angles = acos((x2-x1)./L);
    
    % Order all angles (and save reordered indices)
    [sAngles,idxROs] = sort(angles);
    sortAngles = sAngles((sAngles ~= 0) & abs(sAngles) ~= deg2rad(90));
    idxRO = idxROs.*((sAngles ~= 0) & abs(sAngles) ~= deg2rad(90));
    idxRO = idxRO(idxRO ~= 0);
    
    % Collect all groups of members with angles that deviate less than X
    % degrees from each other
    groupsIdx = {}; numGroups = 1; 
    addingBool = false;
    for i = 1:1:length(sortAngles)
        for j = 1:1:length(sortAngles)
            if addingBool == false
                newGroupIdx = [];
            end
            if ((sortAngles(i) - sortAngles(j)) <= angX) && (i ~= j)
                newGroupIdx = [newGroupIdx,idxRO(j)];
                addingBool = true;
            else
                if addingBool == true
                    addingBool = false;
                    newGroupIdx = [newGroupIdx,idxRO(i)];
                    newGroupIdx = unique(newGroupIdx);
                    groupsIdx(numGroups) = {newGroupIdx};
                    numGroups = numGroups + 1;
                end
            end
        end
    end
    
    % Within each group, catalog all subgroups of members that share a set 
    % of startpoints/endpoints (by member index)
    for k = 1:1:length(groupsIdx)
        currentGroupIdx = cell2mat(groupsIdx(k));
        groupCA = SortedCA(currentGroupIdx,:);
        flowbool = feas_module3_binary(groupCA);
        if (flowbool == 1) 
            stackedCount = stackedCount + 1;
        end
    end

end
