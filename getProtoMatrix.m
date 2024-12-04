function [cfgLDPCEnc, cfgLDPCDec] = getProtoMatrix(block_length, num_info)
    X = -1;
    switch block_length
        case 648
            Z = 27;
            switch num_info
                case 324
                    proto = [...
                             0 X X X 0 0 X X 0 X X 0 1 0 X X X X X X X X X X; ...
                             22 0 X X 17 X 0 0 12 X X X X 0 0 X X X X X X X X X ; ...
                             6 X 0 X 10 X X X 24 X 0 X X X 0 0 X X X X X X X X ; ...
                             2 X X 0 20 X X X 25 0 X X X X X 0 0 X X X X X X X ; ...
                             23 X X X 3 X X X 0 X 9 11 X X X X 0 0 X X X X X X ; ...
                             24 X 23 1 17 X 3 X 10 X X X X X X X X 0 0 X X X X X; ...
                             25 X X X 8 X X X 7 18 X X 0 X X X X X 0 0 X X X X ; ...
                             13 24 X X 0 X 8 X 6 X X X X X X X X X X 0 0 X X X ; ...
                             7 20 X 16 22 10 X X 23 X X X X X X X X X X X 0 0 X X ; ...
                             11 X X X 19 X X X 13 X 3 17 X X X X X X X X X 0 0 X; ...
                             25 X 8 X 23 18 X 14 9 X X X X X X X X X X X X X 0 0; ...
                             3 X X X 16 X X 2 25 5 X X 1 X X X X X X X X X X 0];
                case 432
                    proto = [...
                             25 26 14 X 20 X 2 X 4 X X 8 X 16 X 18 1 0 X X X X X X; ...
                             10 9 15 11 X 0 X 1 X X 18 X 8 X 10 X X 0 0 X X X X X ; ...
                             16 2 20 26 21 X 6 X 1 26 X 7 X X X X X X 0 0 X X X X ; ...
                             10 13 5 0 X 3 X 7 X X 26 X X 13 X 16 X X X 0 0 X X X; ...
                             23 14 24 X 12 X 19 X 17 X X X 20 X 21 X 0 X X X 0 0 X X; ...
                             6 22 9 20 X 25 X 17 X 8 X 14 X 18 X X X X X X X 0 0 X ; ...
                             14 23 21 11 20 X 24 X 18 X 19 X X X X 22 X X X X X X 0 0; ...
                             17 11 11 20 X 21 X 26 X 3 X X 18 X 26 X 1 X X X X X X 0];
                case 486
                    proto = [...
                             16 17 22 24 9 3 14 X 4 2 7 X 26 X 2 X 21 X 1 0 X X X X;...
                             25 12 12 3 3 26 6 21 X 15 22 X 15 X 4 X X 16 X 0 0 X X X  ;...
                             25 18 26 16 22 23 9 X 0 X 4 X 4 X 8 23 11 X X X 0 0 X X;...
                             9 7 0 1 17 X X 7 3 X 3 23 X 16 X X 21 X 0 X X 0  0 X ;...
                             24 5 26 7 1 X X 15 24 15 X 8 X 13 X 13 X 11 X X X X 0 0;...
                             2 2 19 14 24 1 15 19 X 21 X 2 X 24 X 3 X 2 1 X X X X 0];
                case 540
                    proto = [...
                             17 13  8 21  9  3 18 12 10  0  4 15 19  2  5 10 26 19 13 13 1 0 X X ;...
                             03 12 11 14 11 25  5 18  0  9  2 26 26 10 24  7 14 20  4  2 X 0 0 X ;...
                             22 16  4  3 10 21 12  5 21 14 19  5  X  8  5 18 11  5  5 15 0 X 0 0 ;...
                             07  7 14 14  4 16 16 24 24 10  1  7 15  6 10 26  8 18 21 14 1 X X 0 ];
            end

        case 1296
            Z = 54;
            switch num_info
                case 648
                    proto = [ ...
                             40 X X X 22 X 49 23 43 X X X 1 0 X X X X X X X X X X;...
                             50 1 X X 48 35 X X 13 X 30 X X 0 0 X X X X X X X X X;...
                             39 50 X X 4 X 2 X X X X 49 X X 0 0 X X X X X X X X ;...
                             33 X X 38 37 X X 4 1 X X X X X X 0 0 X X X X X X X;...
                             45 X X X 0 22 X X 20 42 X X X X X X 0 0 X X X X X X;...
                             51 X X 48 35 X X X 44 X 18 X X X X X X 0 0 X X X X X;...
                             47 11 X X X 17 X X 51 X X X 0 X X X X X 0 0 X X X X ;...
                             5 X 25 X 6 X 45 X 13 40 X X X X X X X X X 0 0 X X X;...
                             33 X X 34 24 X X X 23 X X 46 X X X X X X X X 0 0 X X;...
                             1 X 27 X 1 X X X 38 X 44 X X X X X X X X X X 0 0 X ;...
                             X 18 X X 23 X X 8 0 35 X X X X X X X X X X X X 0 0;...
                             49 X 17 X 30 X X X 34 X X 19 1 X X X X X X X X X X 0 ];
                case 864
                    proto = [ ...
                             39 31 22 43 X 40 4 X 11 X X 50 X X X 6 1 0 X X X X X X;...
                             25 52 41 2 6 X 14 X 34 X X X 24 X 37 X X 0 0 X X X X X;...
                             43 31 29 0 21 X 28 X X 2 X X 7 X 17 X X X 0 0 X X X X;...
                             20 33 48 X 4 13 X 26 X X 22 X X 46 42 X X X X 0 0 X X X;...
                             45 7 18 51 12 25 X X X 50 X X 5 X X X 0 X X X 0 0 X X;...
                             35 40 32 16 5 X X 18 X X 43 51 X 32 X X X X X X X 0 0 X;...
                             9  24 13 22 28 X X 37 X X 25 X X 52 X 13 X X X X X X 0 0;...
                             32 22 4 21 16 X X X 27 28 X 38 X X X 8 1 X X X X X X 0];
                case 972
                    proto = [...
                             39 40 51 41 3 29 8 36 X 14 X 6 X 33 X 11 X 4 1 0 X X X X;...
                             48 21 47 9 48 35 51 X 38 X 28 X 34 X 50 X 50 X X 0 0 X X X;...
                             30 39 28 42 50 39 5 17 X 6 X 18 X 20 X 15 X 40 X X 0 0 X X;...
                             29 0 1 43 36 30 47 X 49 X 47 X 3 X 35 X 34 X 0 X X 0 0 X;...
                             1 32 11 23 10 44 12 7 X 48 X 4 X 9 X 17 X 16 X X X X 0 0;...
                             13 7 15 47 23 16 47 X 43 X 29 X 52 X 2 X 53 X 1 X X X X 0];
                case 1080
                    proto = [...
                             48 29 37 52 2 16 6 14 53 31 34 5 18 42 53 31 45 X 46 52 1 0 X X;...
                             17 4 30 7 43 11 24 6 14 21 6 39 17 40 47 7 15 41 19 X X 0 0 X;...
                             7 2 51 31 46 23 16 11 53 40 10 7 46 53 33 35 X 25 35 38 0 X 0 0;...
                             19 48 41 1 10 7 36 47 5 29 52 52 31 10 26 6 3 2 X 51 1 X X 0];
            end
        case 1944
            Z = 81;
            switch num_info
                case 972
                    proto = [...
                             57 X X X 50 X 11 X 50 X 79 X 1 0 X X X X X X X X X X;...
                             3 X 28 X 0 X X X 55 7 X X X 0 0 X X X X X X X X X;...
                             30 X X X 24 37 X X 56 14 X X X X 0 0 X X X X X X X X;...
                             62 53 X X 53 X X 3 35 X X X X X X 0 0 X X X X X X X;...
                             40 X X 20 66 X X 22 28 X X X X X X X 0 0 X X X X X X;...
                             0 X X X 8 X 42 X 50 X X 8 X X X X X 0 0 X X X X X;...
                             69 79 79 X X X 56 X 52 X X X 0 X X X X X 0 0 X X X X;...
                             65 X X X 38 57 X X 72 X 27 X X X X X X X X 0 0 X X X;...
                             64 X X X 14 52 X X 30 X X 32 X X X X X X X X 0 0 X X;...
                             X 45 X 70 0 X X X 77 9 X X X X X X X X X X X 0 0 X;...
                             2 56 X 57 35 X X X X X 12 X X X X X X X X X X X 0 0;...
                             24 X 61 X 60 X X 27 51 X X 16 1 X X X X X X X X X X 0];
                case 1296
                    proto = [...
                             61 75 4 63 56 X X X X X X 8 X 2 17 25 1 0 X X X X X X;...
                             56 74 77 20 X X X 64 24 4 67 X 7 X X X X 0 0 X X X X X;...
                             28 21 68 10 7 14 65 X X X 23 X X X 75 X X X 0 0 X X X X;...
                             48 38 43 78 76 X X X X 5 36 X 15 72 X X X X X 0 0 X X X;...
                             40 2 53 25 X 52 62 X 20 X X 44 X X X X 0 X X X 0 0 X X;...
                             69 23 64 10 22 X 21 X X X X X 68 23 29 X X X X X X 0 0 X;...
                             12 0 68 20 55 61 X 40 X X X 52 X X X 44 X X X X X X 0 0;...
                             58 8 34 64 78 X X 11 78 24 X X X X X 58 1 X X X X X X 0];
                case 1458
                    proto = [...
                             48 29 28 39 9 61 X X X 63 45 80 X X X 37 32 22 1 0 X X X X;...
                             4 49 42 48 11 30 X X X 49 17 41 37 15 X 54 X X X 0 0 X X X;...
                             35 76 78 51 37 35 21 X 17 64 X X X 59 7 X X 32 X X 0 0 X X;...
                             9 65 44 9 54 56 73 34 42 X X X 35 X X X 46 39 0 X X 0 0 X;...
                             3 62 7 80 68 26 X 80 55 X 36 X 26 X 9 X X 72 X X X X 0 0;...
                             26 75 33 21 69 59 3 38 X X X 35 X 62 36 26 X X 1 X X X X 0];
                case 1620
                    proto = [...
                             13 48 80 66 4 74 7 30 76 52 37 60 X 49 73 31 74 73 23 X 1 0 X X;...
                             69 63 74 56 64 77 57 65 6 16 51 X 64 X 68 9 48 62 54 27 X 0 0 X;...
                             51 15 0 80 24 25 42 54 44 71 71 9 67 35 X 58 X 29 X 53 0 X 0 0;...
                             16 29 36 41 44 56 59 37 50 24 X 65 4 65 52 X 4 X 73 52 1 X X 0];
            end
    end
    pcmatrix = ldpcQuasiCyclicMatrix(Z,proto);
    cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);
    cfgLDPCDec = ldpcDecoderConfig(pcmatrix,"bp");
end
