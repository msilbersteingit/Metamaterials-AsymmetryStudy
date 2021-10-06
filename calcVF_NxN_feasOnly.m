% FUNCTION TO CALCULATE VOLUME FRACTION 
function volFrac = calcVF_NxN_feasOnly(NC,CA,r,sel,sidenum)
    % Initialize radii vector
    rvar = r.*ones(1,size(CA,1));
    
    % Calculate Avar & modify for edge members
    Avar = pi.*(rvar.^2); % Cross-sectional areas of truss members
    Avar = modifyAreas(Avar,CA,NC,sidenum);

    % Calculate cumulative volume of all individual members
    totalTrussVol = 0;
    totl = sel;
    for i = 1:size(CA,1)
        % Finding element length from nodal coordinates
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        % Adding current element to total volume of members
        totalTrussVol = totalTrussVol + (L*Avar(i));
    end
    
    % Modify total member volume based on overlaps at nodes
    totalTrussVol = subNodOLVol(NC,CA,totalTrussVol,rvar);
    
    % Finding average side "thickness" due to differing element radii
    thick = mean(rvar);
    
    % Calculating volume fraction (using a solid square with 2*(avg 
    %   thickness) as a baseline)
    volFrac = totalTrussVol/(2*thick*(totl^2)); 
end

% FUNCTION TO GENERATE NODAL COORDINATES BASED ON GRID SIZE
function NC = generateNC(sel,sidenum)
    notchvec = linspace(0,1,sidenum);
    NC = [];
    for i = 1:1:sidenum
        for j = 1:1:sidenum
            NC = [NC;notchvec(i),notchvec(j)];
        end
    end
    NC = sel.*NC;
end

% FUNCTION TO MODIFY AREAS FOR EDGE MEMBERS
function Avar = modifyAreas(Avar,CA,NC,sidenum)
    % Identify edge nodes
    edgenodes = [(1:1:sidenum),((2*sidenum):sidenum:...
                 ((sidenum^2)-sidenum)),((sidenum+1):sidenum:...
                 ((sidenum^2)-(2*sidenum)+1)),(((sidenum^2)-sidenum+1):...
                 1:(sidenum^2))];
             
    % Identify members connecting solely to edge nodes
    edgeconn1 = ismember(CA(:,1),edgenodes);
    edgeconn2 = ismember(CA(:,2),edgenodes);
    edgeconnectors = (edgeconn1 & edgeconn2);
    
    % Isolate edge members based on angle
    edgelogical = [edgeconnectors,edgeconnectors];
    CAedgenodes = CA.*(edgelogical);
    CAedgenodes = CAedgenodes(any(CAedgenodes,2),:);
    x1 = NC(CAedgenodes(:,1),1); x2 = NC(CAedgenodes(:,2),1);
    y1 = NC(CAedgenodes(:,1),2); y2 = NC(CAedgenodes(:,2),2);
    L = sqrt(((x2-x1).^2)+((y2-y1).^2));
    angles = rad2deg(abs(acos((x2-x1)./L)));
    CAedgy = [];
    for i = 1:1:size(CAedgenodes,1)
        if (angles(i) == 0) || (angles(i) == 90)
            CAedgy = [CAedgy;CAedgenodes(i,:)];
        end
    end
    
    % Find and modify areas belonging to edge members
    if isempty(CAedgy)
        % Do nothing
    else
        edgemembers = ismember(CA,CAedgy,'rows');
        selectAreas = Avar'.*edgemembers;
        k = find(selectAreas);
        Avar(k) = Avar(k)./2;
    end
end

% FUNCTION TO SUBTRACT VOLUME OVERLAP AT NODES
function tTV = subNodOLVol(NC,CA,tTV,rvar)
    for z = 1:1:size(NC,1)
        % Isolate members originating/ending at the current node
        indione = CA(:,1) == z; inditwo = CA(:,2) == z;
        mCAone = CA(indione,:); mCAtwo = CA(inditwo,:);
        mCA = ...
            [setdiff(mCAone,mCAtwo,'rows');setdiff(mCAtwo,mCAone,'rows')];
        
        % Loop through each pair of connecting members at the current node
        for i = 1:1:size(mCA,1)
            for j = 1:1:size(mCA,1)
                if i ~= j
                    % Find the sweep angle between the connecting members
                    p1 = [NC(mCA(i,1),1),NC(mCA(i,1),2)];
                    q1 = [NC(mCA(i,2),1),NC(mCA(i,2),2)];
                    p2 = [NC(mCA(j,1),1),NC(mCA(j,1),2)];
                    q2 = [NC(mCA(j,2),1),NC(mCA(j,2),2)];
                    v1 = [q1(1),q1(2),0] - [p1(1),p1(2),0];
                    v2 = [q2(1),q2(2),0] - [p2(1),p2(2),0];
                    alpha = ...
                       acos(max(min(dot(v1,v2)/(norm(v1)*norm(v2)),1),-1));
                    theta = pi - alpha;
                    
                    % Check if members overlap
                    if (abs(theta-pi) < (10^-5)) || (abs(theta) < (10^-5))
                        continue
                    else
                        % Find the volume of each overlap sphere sector,
                        % subtract from the total truss volume
                        frac = theta/(2*pi);
                        [~,idx1] = ismember(mCA(i,:),CA,'rows');
                        [~,idx2] = ismember(mCA(j,:),CA,'rows');
                        r1 = rvar(idx1); r2 = rvar(idx2);
                        avgrad = mean([r1,r2]);
                        VOL = ((4*pi)/3)*(avgrad^3)*frac;
                        tTV = tTV - VOL;
                        %disp(tTV);
                    end
                end
            end
        end
    end
end

