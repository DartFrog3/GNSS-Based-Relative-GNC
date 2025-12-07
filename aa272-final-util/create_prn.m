function prn_code = create_prn(idx1, idx2, len)
    % generate prn code as in homework

    % G1 and G2 registers initialized to all ones
    G1 = ones(1, 10, 'uint8');
    G2 = ones(1, 10, 'uint8');
    
    prn_code = zeros(1, len, 'uint8');

    for n = 1:len
        G1_out = G1(end);
        G2_phase = bitxor(G2(idx1), G2(idx2));

        prn_code(n) = bitxor(G1_out, G2_phase);

        G1_feedback = bitxor(G1(10), G1(3));
        G1(2:end) = G1(1:end-1);
        G1(1)  = G1_feedback;

        G2_feedback = mod( ...
            G2(10) + G2(9) + G2(8) + G2(6) + G2(3) + G2(2), 2);
        G2(2:end) = G2(1:end-1);
        G2(1)  = G2_feedback;
    end
end
