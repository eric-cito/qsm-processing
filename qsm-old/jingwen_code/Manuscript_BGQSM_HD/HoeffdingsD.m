function [HoeffdingsD] = HoeffdingsD(X,Y)
%calculate the Hoeffdings' D measurement between two series.
%   X and Y is a N*J matrix, using this function will return HoeffdingsD
%   with length N. For each rows, we will calculate one HoeffdingsD.
[N, J] = size(X);
%initialize;
R = zeros(N,J);
S = zeros(N,J);
Q = zeros(N,J);
for j =1:J
    Xuse=X(:,j);
    Yuse=Y(:,j);
    [ Xsort , pos ]=sort(Xuse);
    Ysort=Yuse(pos);
    R(:,j) = tiedrank(Xsort);    %Fill the j^th R
    S(:,j) = tiedrank(Ysort);    %Fill the j^th S
    
    xLen = numel(Xuse) ;                     %get the length
    x_ties = Xsort(1:xLen-1) >= Xsort(2:xLen);     %get all the point tie to x
    x_tieloc = [find(x_ties); xLen+2];       %reset the site   
    maxTies = numel(x_tieloc);  
    
    [ y_gt , y_et ] = tiedrank_sum(Ysort); % calculate the sum of tied rank
    Q(:,j) = 1 + y_gt + ( 0.5 * y_et ) ;    %ideal statution
    tiecount = 1;                        %remedy for the statution occuring the tie of x
    while (tiecount < maxTies)      
        tiestart = x_tieloc(tiecount);   %tiestart: the stie of tie 
        ntied = 1;                                             %the number of tie
        while(x_tieloc(tiecount+1) == x_tieloc(tiecount)+1)    %jump to next tie
            tiecount = tiecount+1;
            ntied = ntied+1;
        end
        
        Y_temp = Q( (tiestart):(tiestart + ntied) ,j);
        [ y_gt_temp_a , y_et_temp_a ] = tiedrank_sum(Y_temp) ;
        Y_inv_temp = flipud(Y_temp);
        [ y_gt_temp_b , y_et_temp_b ] = tiedrank_sum(Y_inv_temp) ;
        y_gt_temp_b = flipud(y_gt_temp_b);
        y_et_temp_b = flipud(y_et_temp_b);
        Q( (tiestart):(tiestart + ntied) ,j) = Y_temp - 0.5*y_gt_temp_a - 0.25*y_et_temp_a;
        Q( (tiestart):(tiestart + ntied) ,j) = Y_temp + 0.5*y_gt_temp_b + 0.25*y_et_temp_b; 
        tiecount = tiecount + 1;
    end
 
    
    Q(pos,j)=Q(:,j);
    R(pos,j)=R(:,j);
    S(pos,j)=S(:,j);
end
D1 = sum( (Q-1).*(Q-2) );
D2 = sum( (R-1).*(R-2).*(S-1).*(S-2) );
D3 = sum( (R-2).*(S-2).*(Q-1) );
D = 30*((N-2)*(N-3)*D1 + D2 - 2*(N-2)*D3) / (N*(N-1)*(N-2)*(N-3)*(N-4));
HoeffdingsD = D';
  
end
% --------------------------------
function [rsum,t] = tiedrank_sum(x)
% rsum: sum of rank
% t:    number of tie
    [r,t] = Mytr(x);
    rsum=zeros( size(r) );
    n_x = length(x);
    div_a = ceil( sqrt( n_x/(2*log2(n_x)) ) );
    
    %calculate the trQ
    sort_size = floor(n_x/div_a);
    sort_rem = rem( n_x,div_a );
    sort_sum = zeros(sort_size,div_a); %initialize;
    for i = 1:div_a
        temp_x = Mytr( x(1:i*sort_size) );
        sort_sum(:,i) = temp_x( ((i-1)*sort_size + 1) : i*sort_size );
    end
    if sort_rem > 0 
        temp_x = Mytr(x);
        sort_sum_rem = temp_x( (div_a*sort_size+1):n_x );
        n_check = 1;
    else
        n_check = 0;
    end
    
    %calucate each division
    for i = 1:div_a
        i_clock = sort_size;
        rtemp = sort_sum(:,i);
        i_cut = 1;
        start_temp = (i-1)*sort_size;
        while ( i_clock>0 )
            rsum( start_temp + i_clock ) = rtemp(i_clock) - 1 ;
            r_index = rtemp > rtemp(i_clock);
            rtemp = rtemp - r_index;
            %boost inside
            if i_clock == ceil( sort_size - i_cut*sqrt(sort_size) )
                rtemp = rtemp(1:i_clock);
                i_cut = i_cut + 1;
                if i_cut > sqrt( sort_size )
                    i_cut = 0; 
                end
            end
            i_clock = i_clock -1;
        end
    end
    if n_check
        i_clock = sort_rem;
        rtemp = sort_sum_rem;
        i_cut = 1;
        start_temp = div_a * sort_size;
        while ( i_clock >0 )
            rsum( start_temp + i_clock ) = rtemp(i_clock) - 1 ;
            r_index = rtemp > rtemp(i_clock);
            rtemp = rtemp - r_index;
            %boost inside
            if i_clock == ceil( sort_rem - i_cut*sqrt(sort_rem) )
                rtemp = rtemp(1:i_clock);
                i_cut = i_cut + 1;
                if i_cut > sqrt( sort_rem )
                    i_cut = 0; 
                end
            end
            i_clock = i_clock -1;
        end     
    end
end
% --------------------------------
function [r,tievec] = Mytr(x)
%TR Local tiedrank function to compute results for one column
% r:rank
% tievector£ºthe number of tie
% Sort, then leave the NaNs (which are sorted to the end) alone
[sx, rowidx] = sort(x(:));
numNaNs = sum(isnan(x));
xLen = numel(x) - numNaNs;
% Use ranks counting from low end
ranks = [1:xLen NaN(1,numNaNs)]';
tie_temp = zeros( xLen + numNaNs ,1);
% Adjust for ties.  Avoid using diff(sx) here in case there are infs.
ties = sx(1:xLen-1) >= sx(2:xLen);   
tieloc = [find(ties); xLen+2];      
maxTies = numel(tieloc);          
tiecount = 1;
while (tiecount < maxTies)
    tiestart = tieloc(tiecount);   
    ntied = 2;
    while(tieloc(tiecount+1) == tieloc(tiecount)+1)  
        tiecount = tiecount+1;
        ntied = ntied+1;
    end
  
    % Compute mean of tied ranks
    ranks(tiestart:tiestart+ntied-1) = ...
                  sum(ranks(tiestart:tiestart+ntied-1)) / ntied;
    ranks(tiestart:tiestart+ntied-1) = ranks(tiestart:tiestart+ntied-1) - ( ntied -1 )/2;
    tie_temp(tiestart:tiestart+ntied-1) = (1 : ntied)-1  ;
    tiecount = tiecount + 1;
end
% Broadcast the ranks back out, including NaN where required.
r(rowidx) = ranks;
r = reshape(r,size(x));
tievec(rowidx) = tie_temp;
tievec = reshape(tievec,size(x));
end
