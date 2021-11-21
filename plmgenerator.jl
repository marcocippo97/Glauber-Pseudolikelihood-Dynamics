function rand_generator(marg) 
    u = rand()
    F= cumsum(marg)
    @assert F[end] ≈ 1
    flag = false
    k = 0
    while flag == false
        k += 1
        if  u <= F[k]
            flag = true
        end
    end 
    return k
end

function P_marg(i, J, h, seq)  
    N = length(seq)
    q = length(h[1])
    marg = zeros(q)
    somma = 0.0
    for a in 1:q
        ene = h[i][a]
        for j = 1:N
            ene += J[i,j][a, seq[j]]
        end
        marg[a] = exp(ene)
        somma += marg[a]
    end    
    return marg./somma
end

function log_pseudo(seq, J, h) 
    N = length(seq)
    q = length(h[1])
    tot = 0.0 
    for i in 1:N
        tot += h[i][seq[i]]
        for j in 1:N
            tot += J[i,j][seq[i], seq[j]]
        end
        Z = 0.0
        for a in 1:q
            ene = h[i][a]
            for j in 1:N
                ene += J[i,j][a, seq[j]]
            end
            Z += exp(ene)
        end
        tot -= log(Z)
    end
    return tot 
end

function mean_pseudo(natural, J, h) 
    n = length(natural)
    all_pseudo = [log_pseudo(natural[i], J, h) for i in 1:n]
    μ = mean(all_pseudo)
    σ = std(all_pseudo)
    return μ, σ
end

function dist(vec_1, vec_2)
    n = length(vec_1)
    out = 0
    for i in 1:n
        vec_1[i] != vec_2[i] && (out += 1)
    end
    out
end

function mean_dist(vecdata)
    n = length(vecdata)
    media = 0.0
    for i in 1:n
        for j in 1:n
            media += dist(vecdata[i], vecdata[j])
        end
    end
    media /= n*(n-1)
end

#################################################################################

function dynamic(J, h, natural; 
                 r = true, 
                 nterm = 1e4,
                 nsweep = 1e3,
                 nsamp = 1e4)
    N = length(natural[1])
    q = length(h[1])
    all_seq = Array{Int64, 1}[]
    r == true ? seq = rand(1:q, N) : seq = rand(natural)
    for _ in 1:nterm
        i = rand(1:N)
        marg = P_marg(i, J, h, seq)
        seq[i] = rand_generator(marg)
    end
    for _ in 1:nsamp
        for _ in 1:nsweep
            i = rand(1:N)
            marg = P_marg(i, J, h, seq)
            seq[i] = rand_generator(marg)
        end
        push!(all_seq, copy(seq))
    end
    return all_seq
end

function blockmean(data)
    top = length(data)
    bottom = top >> 1
    μ = Float64[]
    σ = Float64[]
    t = Int64[]
    while bottom > 1 
        indx = bottom+1:top
        push!(t, (bottom + top) >> 1)
        val = data[indx] 
        push!(μ, mean(val))
        push!(σ, std(val)/sqrt(length(val)))
        top = bottom
        bottom = bottom >> 1
    end
    return μ, σ, t
end 

function dynamic_term(J, h, natural; r = true, nterm = 1e5)
    N = length(natural[1])
    q = length(h[1])
    all_log = Float64[]
    r == true ? seq = rand(1:q, N) : seq = rand(natural)
    for _ in 1:nterm
        i = rand(1:N)
        marg = P_marg(i, J, h, seq)
        seq[i] = rand_generator(marg)
        push!(all_log, log_pseudo(seq, J, h))
    end
    return all_log
end

#################################################################################

function freq_res(vecdata; W = fill(1/length(vecdata), length(vecdata)), q = 21)
    @assert sum(W) ≈ 1
    N = length(vecdata[1])
    M = length(vecdata)
    freq = zeros(Float64, (N, q))
    for i in 1:N
        for j in 1:M
            freq[i, vecdata[j][i]] += W[j]
        end
    end   
    return freq
end

function freq_res_2(vecdata; W = fill(1/length(vecdata), length(vecdata)), q = 21)
    @assert sum(W) ≈ 1
    N = length(vecdata[1])
    M = length(vecdata)
    freq = zeros(Float64, (N, N, q, q))
    for i in 1:N
        for j in 1:N
            for k in 1:M
                freq[i, j, vecdata[k][i], vecdata[k][j]] += W[k]
            end
        end
    end
    return freq
end

function connected_correlation(vecdata; W = fill(1/length(vecdata), length(vecdata)), q=21)
    @assert sum(W) ≈ 1
    N = length(vecdata[1])
    C = zeros(Float64, (N, N, q, q))
    freq_1 = freq_res(vecdata, W=W, q=q)
    freq_2 = freq_res_2(vecdata, W=W, q=q)
    for i in 1:N
        for j in 1:N
            for k in 1:q
                for l in 1:q
                    C[i, j, k, l] = freq_2[i, j, k, l] - freq_1[i, k]*freq_1[j, l]
                end
            end
        end
    end
    return C
end

function freq_sing_1(vecdata, indx::NTuple{2, Int64}; W = fill(1/length(vecdata), length(vecdata)))
    @assert sum(W) ≈ 1
    M = length(vecdata)
    i, a = indx[1], indx[2]
    freq = 0.0;
    for l in 1:M
        freq += (vecdata[l][i] == a)*W[l];
    end
    return freq
end

function freq_sing_2(vecdata, indx::NTuple{4, Int64}; W = fill(1/length(vecdata), length(vecdata)))
    @assert sum(W) ≈ 1
    M = length(vecdata)
    i, j, a, b = indx[1], indx[2], indx[3], indx[4]
    freq = 0.0;
    for l in 1:M
        freq += (vecdata[l][i] == a)*(vecdata[l][j] == b)*W[l];
    end
    return freq
end

function freq_sing_3(vecdata, indx::NTuple{6, Int64}; W = fill(1/length(vecdata), length(vecdata)))
    @assert sum(W) ≈ 1
    M = length(vecdata)
    i, j, k, a, b, c = indx[1], indx[2], indx[3], indx[4], indx[5], indx[6]
    freq = 0.0;
    for l in 1:M
        freq += (vecdata[l][i] == a)*(vecdata[l][j] == b)*(vecdata[l][k] == c)*W[l];
    end
    return freq
end

function connected_sing_3(vecdata, indx::NTuple{6, Int64}; W = fill(1/length(vecdata), length(vecdata)))
    @assert sum(W) ≈ 1
    i, j, k = indx[1], indx[2], indx[3]
    a, b, c = indx[4], indx[5], indx[6]
    C1 = freq_sing_3(vecdata, (i, j, k, a, b, c), W = W)
    C2 = -freq_sing_2(vecdata, (i, j, a, b), W = W)*freq_sing_1(vecdata, (k, c), W = W)
    C3 = -freq_sing_2(vecdata, (j, k, b, c), W = W)*freq_sing_1(vecdata, (i, a), W = W) 
    C4 = -freq_sing_2(vecdata, (k, i, c, a), W = W)*freq_sing_1(vecdata, (j, b), W = W)
    C5 = 2*freq_sing_1(vecdata, (i, a), W = W)*freq_sing_1(vecdata, (j, b), W = W)*freq_sing_1(vecdata, (k, c), W = W)
    return C1 + C2 + C3 + C4 + C5
end

function test_corr(corr, vecdata; toll = 0.001, W = fill(1/length(vecdata), length(vecdata)), ntest = 100)
    num_ok = 0
    N = size(corr, 1)
    for _ in 1:ntest
        i = rand(1:N)
        corr[i, 7] - toll < connected_sing_3(vecdata, (corr[i, 1], corr[i, 2], corr[i, 3], corr[i, 4], corr[i, 5], corr[i, 6]), W = W) < corr[i, 7] + toll && (num_ok += 1)
    end
    return num_ok/ntest
end

function connected_3(vecdata, M_indx; W = fill(1/length(vecdata), length(vecdata)))
    N = size(M_indx, 1)
    C = zeros(Float64, N)
    for i in 1:N
        C[i] = connected_sing_3(vecdata, (M_indx[i, 1], M_indx[i, 2], M_indx[i, 3], M_indx[i, 4], M_indx[i, 5], M_indx[i, 6]), W = W)
    end
    return C
end

#################################################################################

function autocorrelazione(data, Δt)
    Tmeas = length(data)
    Δt >= Tmeas && error("Δt = $Δt > Tmeas = $Tmeas")
    c0 =  _autocorrelazione(data,0)
    ac = zeros(Δt+1)
    ac[1]=c0
    for dt in 1:Δt
        ac[dt+1] = _autocorrelazione(data,dt)
    end
    return ac./c0
end

function _autocorrelazione(data, Δt) 
    Tmeas = length(data)
    c2  = 0.0
    c1a =  0.0
    c1b =  0.0
    ctr = 0
    for t in 1:Tmeas-Δt
        c2 += data[t] * data[t+Δt]
        c1a += data[t]
        c1b += data[t+Δt]
        ctr += 1
    end
    c2 /= ctr
    c1a /= ctr
    c1b /= ctr
    return (c2 - c1a*c1b)
end

function dist_auto(data, t, Δt)
    n = length(data[1])
    out = 0
    for i in 1:n
        data[t][i] != data[t+Δt][i] && (out += 1)
    end
    out
end

function autocorrelazione_dist(data, Δt)
    Tmeas = length(data)
    Δt >= Tmeas && error("Δt = $Δt > Tmeas = $Tmeas")
    ac = zeros(Δt)
    for dt in 1:Δt
        ac[dt] = _autocorrelazione_dist(data,dt)
    end
    return ac, ac[end]
end

function _autocorrelazione_dist(data, Δt) 
    Tmeas = length(data)
    out  = 0.0
    for t in 1:Tmeas-Δt
        out += dist_auto(data, t, Δt)
    end
    return out /= Tmeas-Δt 
end

#################################################################################

function contrastive(J, h, natural; power_max=4)
    M = length(natural)
    n = length(natural[1])
    out = Array{Array{Int64, 1}}(undef, M, power_max + 1)
    for i in 0:power_max
        for j in 1:M
            seq = copy(natural[j])
            for _ in 1:2^i
                for l in 1:n
                    marg = P_marg(l, J, h, seq)
                    seq[l] = rand_generator(marg)
                end
            end
            out[j, i+1] = seq
        end
    end
    return out
end

function bipartite(vecdata, training_set)
    n = length(vecdata)
    l = length(vecdata[1])
    out = zeros(Int64, n)
    for i in 1:n 
        out[i] = minimum([dist(vecdata[i], x) for x in training_set])
    end
    out/l
end

function get_first(vec, m; toll = 10^-2, point = 0.8)
    n = length(vec)
    for i in 1:n
       point - toll < vec[i] < point + toll && return i/m
    end
    return println("non trovato")
end

function mean_change(vecdata_1, vecdata_2)
    @assert length(vecdata_1) == length(vecdata_1)
    n = length(vecdata_1)
    l = length(vecdata_1[1])
    out = 0.0
    for i in 1:n
        out += dist(vecdata_1[i], vecdata_2[i])
    end
    out/(n*l)    
end

function num_contact(file::String)
        M = readdlm(file)
        n = 0
        for i in 1:size(M, 1)
            if abs(M[i, 1] - M[i, 2]) > 4 && M[i, 4] < 7
                n += 1
            end
        end
        n
    end