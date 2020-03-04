function possible_sequence = nucleotide_seq(desired_sequence)
    % Generate nt to aa map
    map = revgeneticcode(1, 'Alphabet', 'DNA');
    % Desired AA should be in three code format
    possible_sequence = {};
    for (i=1:strlength(desired_sequence))
        desired_aa = desired_sequence(i);
        switch desired_aa
            case 'a'
                possible_nt_seq = map.A;
            case 'r'
                possible_nt_seq = map.R;
            case 'n'
                possible_nt_seq = map.N;
            case 'd'
                possible_nt_seq = map.D;
            case 'c'
                possible_nt_seq = map.C;
            case 'q'
                possible_nt_seq = map.Q;
            case 'e'
                possible_nt_seq = map.E;
            case 'g'
                possible_nt_seq = map.G;
            case 'h'
                possible_nt_seq = map.H;
            case 'i'
                possible_nt_seq = map.I;
            case 'l'
                possible_nt_seq = map.L;
            case 'k'
                possible_nt_seq = map.K;
            case 'm'
                possible_nt_seq = map.M;
            case 'f'
                possible_nt_seq = map.F;
            case 'p'
                possible_nt_seq = map.P;
            case 's'
                possible_nt_seq = map.S;
            case 't'
                possible_nt_seq = map.T;
            case 'w'
                possible_nt_seq = map.W;
            case 'y'
                possible_nt_seq = map.Y;
            case 'v'
                possible_nt_seq = map.V;
            case 'Starts'
                possible_nt_seq = map.Starts;
            case 'Stops'
                possible_nt_seq = map.Stops;
            otherwise
                possible_nt_seq = "Please input valid sequence";
            
        end
        possible_sequence = vertcat(possible_sequence, possible_nt_seq');
    end
        
    

end