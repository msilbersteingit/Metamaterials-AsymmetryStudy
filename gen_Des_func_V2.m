%%% GENERATIVE DESIGN FUNCTION V2
function CA_des_all = gen_Des_func_V2(sidenum,sel,n_des,short_member_prob,short_members_only,numperdes)
    NC = generate_NC(sel,sidenum);
    n_total_members = nchoosek(sidenum*sidenum,2);
    CA_all = generate_full_CA(sidenum);

    % Create node labels
    labels = struct;
    for i = 1:sidenum*sidenum
        labels.(strcat('label',num2str(i))) = num2str(i);
    end

    % Generating designs
    CA_des_all = {};
    des_count = 0;

    all_satisfactory_designs_found = false;
    while ((des_count <= n_des) || ~all_satisfactory_designs_found)
        CA_des_current = zeros(n_total_members,2);
        member_count = 0;
        %{
        % Visualization of design generation
        figure
        % Plot node positions
        for k = 1:size(NC,1)
            plot(NC(k,1),NC(k,2),'*r')
            hold on
            text(NC(k,1),NC(k,2),labels.(strcat('label',num2str(k))),'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15)
            hold on
            drawnow limitrate
        end
        %}
        % Add random member
        %rand_member_index = randi(size(CA_all,1));
        %rand_member = CA_all(rand_member_index,:);

        % Choose random member with probability biasing towards shorter
        % members
        first_member = true;
        prob_array = get_probabilities(CA_all,NC,sidenum,sel,short_members_only,short_member_prob,first_member);
        feas_indices = 1:1:size(CA_all,1);
        rand_index = randsample(feas_indices,1,true,prob_array);
        rand_member = CA_all(rand_index,:);

        CA_des_current(member_count+1,:) = rand_member;
        member_count = member_count + 1;
        %{
        % Plotting member
        x1 = NC(CA_des_current(member_count,1),1);
        y1 = NC(CA_des_current(member_count,1),2);
        x2 = NC(CA_des_current(member_count,2),1);
        y2 = NC(CA_des_current(member_count,2),2);
        line1 = plot([x1,x2],[y1,y2],'-b','LineWidth',2);
        drawnow limitrate
        hold on
        %}
        repeat_member = get_corresponding_repeated_member(rand_member,sidenum);
        if (~isempty(repeat_member))
            CA_des_current(member_count+1,:) = repeat_member;
            member_count = member_count + 1;
            %{
            % Plotting repeat member
            x1 = NC(CA_des_current(member_count,1),1);
            y1 = NC(CA_des_current(member_count,1),2);
            x2 = NC(CA_des_current(member_count,2),1);
            y2 = NC(CA_des_current(member_count,2),2);
            line2 = plot([x1,x2],[y1,y2],'-b','LineWidth',2);
            drawnow limitrate
            hold on
            %}
        end

        % Start generative design process
        random_member_old = rand_member;
        satisfactory_design = false;
        first_member = false;

        while ~satisfactory_design

            % Add feasible connected member at random
            conn_members_rand = get_feasible_connected_members(random_member_old,CA_des_current(1:member_count,:),CA_all,NC,sidenum,sel,short_members_only);

            % If no feasible connected member is found, find feasible connected
            % members at node with lowest connections
            if (isempty(conn_members_rand))

                [n_conn_nodes_current, ~] = connectivityCounter(sidenum,CA_des_current(1:member_count,:),NC,sel);

                % Add up counters based on nodal connectivities (sans repeatability)
                %binedges = (1:1:(size(NC,1)+1)) - 0.5;
                %[n_conn_nodes_current,~] = histcounts(CA_des_current(1:member_count,:),binedges);

                unused_nodes = n_conn_nodes_current == 0;
                n_conn_nodes_current(unused_nodes) = 50; 
                [~,min_index] = min(n_conn_nodes_current);
                conn_members_rand = get_feasible_connected_members(min_index,CA_des_current(1:member_count,:),CA_all,NC,sidenum,sel,short_members_only);
            end

            %rand_member_index = randi(size(conn_members_rand,1));
            %random_member_new = conn_members_rand(rand_member_index,:);

            % Choose random member with probability biasing towards shorter
            % members
            prob_array = get_probabilities(conn_members_rand,NC,sidenum,sel,short_members_only,short_member_prob,first_member);
            feas_indices = 1:1:size(conn_members_rand,1);
            rand_index = randsample(feas_indices,1,true,prob_array);
            random_member_new = conn_members_rand(rand_index,:);

            CA_des_current(member_count+1,:) = random_member_new;
            member_count = member_count + 1;
            %{
            % Plotting member
            x1 = NC(CA_des_current(member_count,1),1);
            y1 = NC(CA_des_current(member_count,1),2);
            x2 = NC(CA_des_current(member_count,2),1);
            y2 = NC(CA_des_current(member_count,2),2);
            line2 = plot([x1,x2],[y1,y2],'-b','LineWidth',2);
            drawnow limitrate
            hold on
            %}
            repeat_member_new = get_corresponding_repeated_member(random_member_new,sidenum);
            if (~isempty(repeat_member_new))

                CA_des_current(member_count+1,:) = repeat_member_new;
                member_count = member_count + 1;

                % Check for intersection and overlap satisfaction
                no_inters = feas_module1_binary(CA_des_current(1:member_count,:),NC,sel);
                no_overlaps = feas_module2_binary(CA_des_current(1:member_count,:),NC,sel);

                if (no_inters && no_overlaps)
                    %{
                    % Plotting repeated member
                    x1 = NC(CA_des_current(member_count,1),1);
                    y1 = NC(CA_des_current(member_count,1),2);
                    x2 = NC(CA_des_current(member_count,2),1);
                    y2 = NC(CA_des_current(member_count,2),2);
                    line2_rep = plot([x1,x2],[y1,y2],'-b','LineWidth',2);
                    drawnow limitrate
                    hold on
                    %}
                else
                    % Remove both previous members since repeated member is not
                    % feasible
                    CA_des_current(member_count-1:member_count,:) = zeros(2,2);
                    member_count = member_count - 2;

                    %delete(line2)
                    %drawnow limitrate
                    %hold on
                end
            end

            disp(strcat('Number of members added: ',num2str(member_count)))

            % Check for intercell and intracell connectivity and nodal
            % connection rules satisfaction
            intracell_conn = feas_module3_binary(CA_des_current(1:member_count,:));
            intercell_conn = feas_module4_binary(CA_des_current(1:member_count,:),sidenum);
            %{
            nodes_rule_sat = true;
            [n_conn_nodes, n_holes] = connectivityCounter(sidenum,CA_des_current(1:member_count,:),NC,sel);

    %         for j = 1:length(n_conn_nodes)
    %             if (n_conn_nodes(j) <= 2 && n_conn_nodes(j) > 0)
    %                 nodes_rule_sat = false;
    %                 break
    %             end
    %         end

            if(any(n_conn_nodes == 1))
                nodes_rule_sat = false;
            end
            %}
            % If the three rules are satisfied, deem the design as satisfactory
            %if(intercell_conn && intracell_conn && nodes_rule_sat)
            if(intercell_conn && intracell_conn)
                satisfactory_design = true;
                disp('Satisfactory design found')
                hold off
            end

            random_member_old = random_member_new;
        end

        % Add the generated design to the design structure
        CA_des_sorted = sortrows(CA_des_current(1:member_count,:));
        des_count = des_count + 1;
        CA_des_all(des_count) = {CA_des_sorted};

        % Keep adding feasible members till list of feasible members is
        % exhausted
        additional_feas_members = get_all_feasible_connections(CA_des_sorted,CA_all,NC,sidenum,sel,short_members_only);
        all_feas_exhausted = false;
        CA_des_add = {};
        des_created = 0;
        while ~all_feas_exhausted
            % Select first member from array
            %current_feas_member = additional_feas_members(1,:);

            % Select member at random from array
            prob_array = get_probabilities(additional_feas_members,NC,sidenum,sel,short_members_only,short_member_prob,first_member);
            feas_indices = 1:1:size(additional_feas_members,1);
            rand_index = randsample(feas_indices,1,true,prob_array);
            current_feas_member = additional_feas_members(rand_index,:);

            % Add member to last design from series
            CA_des_new = vertcat(CA_des_sorted,current_feas_member);
            repeat_member = get_corresponding_repeated_member(current_feas_member,sidenum);
            if (~isempty(repeat_member))
                CA_des_new = vertcat(CA_des_new,repeat_member);
            end
            %visualize_CA(CA_des_new,NC,labels)

            % Save design
            des_created = des_created + 1;
            CA_des_sorted = sortrows(CA_des_new);
            CA_des_add(des_created) = {CA_des_sorted};

            disp('Additional candidate design created')
            additional_feas_members = get_all_feasible_connections(CA_des_sorted,CA_all,NC,sidenum,sel,short_members_only);

            % Complete loop if all feasible members are used up, and select
            % some of the new designs to save permanently
            if isempty(additional_feas_members)
                all_feas_exhausted = true;
                if des_created <= numperdes
                    for j = 1:1:des_created
                        des_count = des_count + 1;
                        CA_des_all(des_count) = CA_des_add(des_created);
                    end
                else
                    idxbin = randperm(des_created,numperdes);
                    for k = 1:1:length(idxbin)
                        des_count = des_count + 1;
                        CA_des_all(des_count) = CA_des_add(idxbin(k));
                    end
                end
            end
        end
        all_satisfactory_designs_found = true;
    end
end

%%% REQUISITE SUBFUNCTIONS
function NC = generate_NC(sel,sidenum)
    notchvec = linspace(0,1,sidenum);
    NC = [];
    for i = 1:1:sidenum
        for j = 1:1:sidenum
            NC = [NC;notchvec(i),notchvec(j)];
        end
    end
    NC = sel.*NC;
end

function CA_full = generate_full_CA(sidenum)
    member_count = 0;
    n_total_members = nchoosek(sidenum*sidenum,2);
    CA_full = zeros(n_total_members,2);
    for i = 1:((sidenum*sidenum)-1)
        for j = i+1:(sidenum*sidenum)
            CA_full(member_count+1,:) = [i,j];
            member_count = member_count + 1;
        end
    end
end

function repeat_member = get_corresponding_repeated_member(member,sidenum)
    repeat_member = [];
    left_edge_nodes = 1:1:sidenum;
    bottom_edge_nodes = 1:sidenum:(sidenum*sidenum - sidenum)+1;
    top_edge_nodes = sidenum:sidenum:sidenum*sidenum+1;
    right_edge_nodes = (sidenum*sidenum - sidenum)+1:1:sidenum*sidenum;
    
    if (all(ismember(member,left_edge_nodes))) % member is a left edge member
        repeat_member = member + (sidenum*sidenum - sidenum);
    elseif (all(ismember(member,bottom_edge_nodes))) % member is a bottom edge member
        repeat_member = member + sidenum - 1;
    elseif (all(ismember(member,right_edge_nodes))) % member is a right edge member
        repeat_member = member - (sidenum*sidenum - sidenum);
    elseif (all(ismember(member,top_edge_nodes))) % member is a top edge member
        repeat_member = member - (sidenum - 1);
    end
end

function prob_array = get_probabilities(feas_members,NC,sidenum,sel,only_short_members,prob_short_member,first_member_addition)
    prob_array = zeros(1,size(feas_members,1));
    
    if (only_short_members && ~first_member_addition)
        prob_array = 1/size(feas_members,1).*ones(1,size(feas_members,1));
    elseif (~only_short_members || first_member_addition)
        short_member_indices = zeros(1,size(feas_members,1));
        long_member_indices = zeros(1,size(feas_members,1));
        n_short_members = 0;
        n_long_members = 0;
        for i = 1:size(feas_members,1)
            current_member = feas_members(i,:);
            node1_pos = [NC(current_member(1),1), NC(current_member(1),2)];
            node2_pos = [NC(current_member(2),1), NC(current_member(2),2)];
            member_length = sqrt((node1_pos(1) - node2_pos(1))^2 + (node1_pos(2) - node2_pos(2))^2);
            if member_length <= (sqrt(2)+1e-5)*(sel/(sidenum - 1)) % additional 1e-5 to resolve precision issue
                n_short_members = n_short_members + 1;
                short_member_indices(n_short_members) = i;
            else
                n_long_members = n_long_members + 1;
                long_member_indices(n_long_members) = i;
            end
        end
        short_member_indices_red = short_member_indices(1:n_short_members);
        long_member_indices_red = long_member_indices(1:n_long_members);
        
        if only_short_members
            prob_array(short_member_indices_red) = 1/n_short_members;
            prob_array(long_member_indices_red) = 0;
        else
            prob_array(short_member_indices_red) = prob_short_member/n_short_members;
            prob_array(long_member_indices_red) = (1 - prob_short_member)/n_long_members;
        end
    end
    
end

function conn_members_red = get_feasible_connected_members(member,CA_des,CA_full,nodal_conn,sidenum,sel,only_short_members)
    % CA_des is the current CA without the connections
    feas_conn_members = zeros(size(CA_full,1),2);
    n_feas_conn_members = 0;
    for i = 1:size(CA_full,1)
        if any(ismember(member,CA_full(i,:)))
            current_member = CA_full(i,:);
            CA_des_current = vertcat(CA_des,current_member);
            
            % Check for intersection and overlap satisfaction
            no_inters = feas_module1_binary(CA_des_current,nodal_conn,sel);
            no_overlaps = feas_module2_binary(CA_des_current,nodal_conn,sel);
        
            if (no_inters && no_overlaps)
                if only_short_members
                    node1_pos = [nodal_conn(current_member(1),1), nodal_conn(current_member(1),2)];
                    node2_pos = [nodal_conn(current_member(2),1), nodal_conn(current_member(2),2)];
                    member_length = sqrt((node1_pos(1) - node2_pos(1))^2 + (node1_pos(2) - node2_pos(2))^2);
                    if member_length <= (sqrt(2)+1e-5)*(sel/(sidenum - 1)) % additional 1e-5 to resolve precision issue
                        feas_conn_members(n_feas_conn_members+1,:) = current_member;
                        n_feas_conn_members = n_feas_conn_members + 1;
                    end
                else
                    feas_conn_members(n_feas_conn_members+1,:) = current_member;
                    n_feas_conn_members = n_feas_conn_members + 1;
                end
            end
        end
    end
    conn_members_red = feas_conn_members(1:n_feas_conn_members,:);
end

function all_feasible_members_un = get_all_feasible_connections(CA_sat,CA_complete,nodal_conn,sidenum,sel,only_short_members)
    all_feasible_members = zeros(size(CA_complete,1),2);
    n_feas_conns = 0;
    for i = 1:size(CA_sat,1)
        current_member = CA_sat(i,:);
        current_conn_members = get_feasible_connected_members(current_member,CA_sat,CA_complete,nodal_conn,sidenum,sel,only_short_members);
        all_feasible_members(n_feas_conns+1:n_feas_conns+size(current_conn_members,1),:) = current_conn_members; 
        n_feas_conns = n_feas_conns + size(current_conn_members,1);
    end
    all_feasible_members_red = all_feasible_members(1:n_feas_conns,:);
    all_feasible_members_un = unique(all_feasible_members_red,'rows');
end

% Function creates a figure visualizing the corresponding CA
function [] = visualize_CA(CA_current,NC,labels_struct)
    figure
    % Plot node positions
    for k = 1:size(NC,1)
        plot(NC(k,1),NC(k,2),'*r')
        hold on
        text(NC(k,1),NC(k,2),labels_struct.(strcat('label',num2str(k))),'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15)
        hold on
        drawnow limitrate
    end
    % Plot the members of CA_current
    for m = 1:size(CA_current,1)
        x1 = NC(CA_current(m,1),1);
        y1 = NC(CA_current(m,1),2);
        x2 = NC(CA_current(m,2),1);
        y2 = NC(CA_current(m,2),2);
        plot([x1,x2],[y1,y2],'-b','LineWidth',2);
        drawnow limitrate
        hold on
    end
    hold off
end