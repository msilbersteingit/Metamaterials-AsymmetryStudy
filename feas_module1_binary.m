% FIRST CONSTRAINT: members only intersect at nodes (no crossing)
% nocrossbool = true for no crossing, false for at least one crossing
% present
function nocrossbool = feas_module1_binary(CA,NC,sel)
    % Initialize values
    nocrossbool = true;
    SortedCA = sortrows(CA);
    ND = NC./sel;
    
    % Develop 4xM matrix of line segment endpoint coordinates, where M is 
    %   the number of truss members.  Each row of format (x1,y1,x2,y2),
    %   where point 1 is leftmost, point 2 is rightmost
    PosA = [ND(SortedCA(:,1),1),ND(SortedCA(:,1),2),...
            ND(SortedCA(:,2),1),ND(SortedCA(:,2),2)];
    
    % Loop through each pair of elements
    for i = 1:1:size(PosA,1)
        for j = 1:1:size(PosA,1)
            %printstring = strcat('i and j: ',num2str(i),', ',num2str(j));
            %disp(printstring);
            % Determine whether the given pair of elements intersects
            intersect = findLineSegIntersection([PosA(i,1),PosA(i,2)],...
                        [PosA(i,3),PosA(i,4)],[PosA(j,1),PosA(j,2)],...
                        [PosA(j,3),PosA(j,4)]);
            
            % Throw an error, given an intersection
            if intersect == true
                nocrossbool = false;
                return
            end
        end
    end
end

% FUNCTION TO DETERMINE PRESENCE OF INTERSECTION
function intersect = findLineSegIntersection(p1,q1,p2,q2)
    if (findOrientation(p1,q1,p2) ~= findOrientation(p1,q1,q2))&&...
            (findOrientation(p2,q2,p1) ~= findOrientation(p2,q2,q1)) 
        if isequal(p1,p2) || isequal(q1,q2) || ...
                isequal(p1,q2) || isequal(q1,p2)
            % Intersection due to shared start or end point
            intersect = false;
        else
            intersect = true;
        end
    else
        % Intersection not present
        intersect = false;
    end
end

% FUNCTION TO CALCULATE ORIENTATION FROM 3 POINTS
function orientation = findOrientation(p,q,r)
    val = ((q(2)-p(2))*(r(1)-q(1)))-((q(1)-p(1))*(r(2)-q(2)));
    %disp(strcat('value = ',num2str(val)));
    if val == 0
        orientation = 0;
    elseif val > 0
        orientation = 1;
    else
        orientation = 2;
    end
    %disp(strcat('orientation = ',num2str(orientation)));
end