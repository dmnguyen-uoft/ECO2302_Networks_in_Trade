function out = network(N,p)

    %% generate network

    % initialize adjacency matrix
    M = zeros(N);

    % form edges
    for i = 1:N-1
        for j = i+1:N

            % form edge if uniform RV on [0,1] is no greater than formation probability p
            M(i,j) = rand <= p; 
            
            % form reverse edge
            M(j,i) = M(i,j);

        end
    end

    
    %% compute components
    
    % initialize set of visited nodes
    visited = zeros(N,1);
    
    % initialize components
    C = {};
    
    % cycle through nodes
    for n = 1:N
        
        % if node has not been visited
        if ~visited(n)
           
            % flag node as visited
            visited(n) = 1;
            
            % start new component and add current node as member
            members = n;
            
            % initialize index
            index = 1;
            
            while index <= length(members)
            
                % find neighbors of current node
                nbrs = find(M(:,members(index)));
                
                % find unvisited neighbors
                new_nbrs = nbrs(visited(nbrs)==0);
                
                % add unvisited neighbors to member list
                members(end+1:end+length(new_nbrs)) = new_nbrs;
                
                % flag visited neighbors as visited
                visited(new_nbrs) = 1;
                
                % increment index
                index = index + 1;
                
            end
            
            % store members of current component
            C{end+1} = members; %#ok<AGROW>
           
        end
        
    end
    
     %% compute distances
     
     % initialize distance matrix
     D = zeros(N);
     
     % compute distances for each component
     Ncomps = length(C);
     for nc = 1:Ncomps
        
         % get component
         comp = C{nc};
         Nc = length(comp);
         
         % if component has more than one node
         
         if Nc > 1
             
            % get adjacency matrix for current component
            Mc = M(comp,comp);
            
            % compute eigenvector decomposition of adjacency matrix
            [V,Lam] = eig(Mc);
            Vinv = inv(V);
            
            % initialize distance matrix
            Dc = Inf*ones(Nc);
            
            % set self distance to zero
            for i = 1:Nc
                Dc(i,i) = 0;
            end
            
            % initialize path length
            k = 1;
            Lam_k= eye(Nc);

            while max(Dc(:)) == Inf
                %k

                % compute M^k
                Lam_k = Lam_k*Lam;   
                Mck = V*Lam_k*Vinv;

                % store path length for connected nodes
                for i = 1:Nc-1
                    for j = i+1:Nc

                        if Mck(i,j) > 1E-3 && Dc(i,j) == Inf
                            Dc(i,j) = k;
                            Dc(j,i) = k;
                        end

                    end
                end

                % increment path length
                k = k+1;

            end
            
            % store distances for current component
            D(comp,comp) = Dc;
            
         end
         
     end
    
    %% store results
    
    out.M = M;
    out.D = D;
    out.C = C;


end