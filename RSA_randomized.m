%
% RSA simple example with randomized p, q, e
%
clear;
clc;
%
% purpose: Alice wanna sends a message to Bob
%
% step 1: Bob choose two very large distinct primes p and q, and computes 
%         n = p * q and £p(n) = (p-1) * (q-1)
%
% step 2: Bob choose a random integer e, 1 < e < phi_n, such that gcd(e, £p(n)) = 1 
%
% step 3: Bib finds d, 1 < d < £p(n), such that e * d = 1 mod £p(n)
%
% step 4: Bob makes (n, e) public, and keeps (p, q, d) private
%
% step 5: Alice computes c = m^e mod n, 1 < m < n, and sends c to Bob
%
% step 6: Bob decrypts c by computing c^d mod n
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% step 1: Bob choose two very large distinct primes p and q, and computes 
%         n = p * q and phi_n = (p-1) * (q-1)
%
% generate p, 2^(L-1) < p < 2^L
L = 13; % 13
p = prime_generator(L);
% generate q, 2^(N-1) < q < 2^N
N = L + 1; % 12, 14
q = prime_generator(N);
%
n = p * q;
phi_n = (p-1) * (q-1);
%
% step 2: Bob choose a random integer e, 1 < e < £p(n), such that gcd(e, £p(n)) = 1
%
flag = 0;
while flag == 0
    %
    e = floor((phi_n - 2) * rand(1)) + 2;
    %
    % compute gcd(e, £p(n))
    %
    % Algorithm 2.19 Extended Euclidean algorithm for integers
    %
    ue = e;     % u
    ve = phi_n; % v
    x1e = 1;
    y1e = 0;
    x2e = 0;
    y2e = 1;
    while ue ~= 0
        % 3.1
        qe = floor(ve / ue); % qq
        re = ve - qe * ue;
        xe = x2e - qe * x1e;
        ye = y2e - qe * y1e;
        % 3.2
        ve = ue;
        ue = re;
        x2e = x1e;
        x1e = xe;
        y2e = y1e;
        y1e = ye;
    end
    % gcd, Greatest Common Divisor
    gcd = ve;
    xe = x2e; % x
    % ye = y2e;
    if gcd == 1
        flag = 1;
    end
end
%
% step 3: Bob finds d, 1 < d < £p(n), such that e * d = 1 mod £p(n)
%         because m^£p(n) = 1 mod n
% check if xe is positive, cause we take xe positive only
if xe < 0
    d = phi_n + xe;
else
    d = xe;
end
%
% step 4: Bob makes (n, e) public, and keeps (p, q, d) private
%
% step 5: Alice computes c = m^e mod n, 1 < m < n, and sends c to Bob
%
% input massage, valid ASCII code for texts or symbols range from 32 to 126
%
k = floor(log10(n) / log10(95)); % 126 - 31 = 95, (126 - 31) ^ k <= n
% plaintext = 'Do not judge a book by its cover. Do not judge me from my outside.';
plaintext = '11223344556677889900aabbccddeeff';
plaintext_dec = double(plaintext) - 31;
plaintext_len = length(plaintext_dec);
m = zeros(1, ceil(plaintext_len / k));
for i = 1 : floor(plaintext_len / k)
    for j = 1 : k
        m(i) = m(i) + plaintext_dec((i-1)*k + j) * (95)^(j-1); % 126 - 31 = 95
    end
end
%
if mod(plaintext_len, k) ~= 0
    i = i + 1;
    for j = 1 : mod(plaintext_len, k)
        m(i) = m(i) + plaintext_dec((i-1)*k + j) * (95)^(j-1); % 126 - 31 = 95
    end
end
%
% compute c = m^e mod n
%
c = ones(ceil(plaintext_len / k), 1);
for i = 1 : ceil(plaintext_len / k)
    ee = e;    % keep e    unchanged, use ee for following computation
    mm = m(i); % keep m(i) unchanged, use mm for following computation
    while ee ~= 0
        if mod(ee, 2) == 1
            c(i) = mod(c(i) * mm, n);
        end
        ee = floor(ee / 2);
        mm = mod(mm * mm, n);
    end
end
%
% step 6: Bob decrypts c by computing c^d mod n
%
m_size = size(m, 2); % size of m
rm = ones(1, m_size);
for i = 1 : m_size
    dd = d;    % keep d    unchanged, use dd for following computation
    cc = c(i); % keep c(i) unchanged, use cc for following computation
    while dd ~= 0
        if mod(dd, 2) == 1
            rm(i) = mod(rm(i) * cc, n);
        end
        dd = floor(dd / 2);
        cc = mod(cc * cc, n);
    end
end
%
R_plaintext = zeros(1, 1);
j = 1;
for i = 1 : m_size
    rmm = rm(i); % keep rm(i) unchanged, use rmm for following computation
    while rmm ~= 0
        R_plaintext(1, j) = mod(rmm, 126 - 31);
        rmm = floor(rmm / (126-31));
        j = j + 1;
    end
end
R_plaintext = R_plaintext + 31; 
R_plaintext = char(R_plaintext);
%
% print out the resulting data
%
fprintf('the original  message is: %s \n', plaintext);
fprintf('the recovered message is: %s \n', R_plaintext);
%
fprintf('the public keys [n, e] are: [%d, %d] \n', n, e);
fprintf('the private keys [p, q, d] are: [%d, %d, %d] \n', p, q, d);

