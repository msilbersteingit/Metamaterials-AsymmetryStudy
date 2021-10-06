% SECOND CONSTRAINT: Elements (of either the same or different lengths)
% nooverlapbool = true for a design with no overlapping members and false
% for a design that contains at least one instance of overlapping members
function nooverlapbool = feas_module2_binary(CA,NC,sel)
    % Initialize values
    SortedCA = sortrows(CA);
    ND = NC./sel;
    nooverlapbool = 1;

    % Loop through each element
    for k = 1:1:size(SortedCA,1)
        % Loop through each element again, to consider each possible pair 
        %   of elements
        for q = 1:1:size(SortedCA,1)
            % Isolate startpoint/endpoint coordinates of both members,
            % calculate their slopes
            A = ND(SortedCA(k,1),:); B = ND(SortedCA(k,2),:);
            C = ND(SortedCA(q,1),:); D = ND(SortedCA(q,2),:);
            mk = (B(2)-A(2))/(B(1)-A(1));
            mq = (D(2)-C(2))/(D(1)-C(1));
            mk = round(mk,4); mq = round(mq,4);
            
            % Check if the same element is being compared twice
            if k == q
                continue
            % Check if the elements' slopes are the same (i.e. they are
            % parallel)
            elseif (mk == mq)
                % Check if the elements are horizontal and their
                % x-coordinates overlap
                if (mk == 0) && ((C(1) >= A(1)) && (C(1) < B(1)))
                    % Check if the elements' y-coordinates overlap
                    if C(2) == A(2)
                        nooverlapbool = 0;
                        return
                    end
                % Check if the elements are vertical and their coordinates
                % overlap
                elseif isinf(mk) && ((C(2) >= A(2)) && (C(2) < B(2)))
                    % Check if the elements' x-coordinates overlap
                    if C(1) == A(1)
                        nooverlapbool = 0;
                        return
                    end
                else
                    if k < q
                        t1 = (C(1)-A(1))/(B(1)-A(1));
                        t2 = (C(2)-A(2))/(B(2)-A(2));
                        t3 = (D(1)-A(1))/(B(1)-A(1));
                    elseif k > q
                        t1 = (A(1)-C(1))/(D(1)-C(1));
                        t2 = (A(2)-C(2))/(D(2)-C(2));
                        t3 = (B(1)-C(1))/(D(1)-C(1));
                    end
                    % Check if the diagonal elements overlap
                    if (t1 == t2)
                        if (t1 >= 0) && (t1 < 1)
                            nooverlapbool = 0;
                            return
                        elseif (t3 >= 0) && (t3 < 1)
                            nooverlapbool = 0;
                            return
                        end
                    end
                end
            end
        end
    end
end

