%
% prime_generator function
%
function p = prime_generator(L)
    %
    % prime generation by Miller-Rabin primality test
    %
    % prime: p, 2^(L-1) < p < 2^L
    % L = 27; % for constant time L = 28 would take too long to compute
    %
    % on the basis of Miller-Rabin test, we have 25% chance that p is not a
    % prime in a single test. if we constantly repeat the test for like 40
    % times to a numer p, and if every test result all states "inconclusive," 
    % then the probability of p is not a prime number is close to (0.25)^40 ~
    % 2^(-80) ~ 8.271806125530277e-25, so we can firmly believe that this
    % number p is a prime number!
    %
    starting_time = cputime;
    %
    iteration = 4000; % max iteration == 40
    %
    % twoL1 = 2^(L-1), the lowerbound of p
    %
    L1 = L - 1;
    twoL1 = 1;
    twoL1q = floor(L1 / 32);
    twoL1r = mod(L1, 32);
    for i = 1 : twoL1q
        twoL1 = twoL1 * 2^32;
    end
    twoL1 = twoL1 * 2^twoL1r;
    %
    % twoL = 2^L, the upperbound of p
    %
    twoL = twoL1 * 2;
    %
    % randomly select a number p, 2^(L-1) < p < 2^L
    %
    dd = twoL - twoL1 - 1;
    %
    flag_prime = 0; % set the flag to 0, means False
    while flag_prime == 0
        p = floor(dd * rand(1)) + 1;
        p = p + twoL1;
        if mod(p, 2) == 0
            p = p - 1; % p must be odd, for obvious reason
        end
        %
        it = 0; % iteration index
        result = 'inconclusive';  % default setting 'inconclusive'
        while strcmp(result, 'inconclusive') == 1 && it <= iteration
            it = it + 1;
            result = Miller_Rabin_test(p);
        end
        if strcmp(result, 'inconclusive') == 1
            % after we find the possible valid prime number
            flag_prime = 1; % set the flag to 1, means True
        end
    end
    %
    ending_time = cputime;
    comp_time = ending_time - starting_time;
    %
    fprintf('the current prime is: %d\n', p);
    fprintf('the computation time is: %f\n', comp_time);
    %
    % testing results:
    %     p = 0917505
    %         0786433
    %         1048561
    %         1015809
    %         1048545
    %         
return
%
% L = 28
% the current prime is: 139081193
% the computation time is: 120.015625
%

