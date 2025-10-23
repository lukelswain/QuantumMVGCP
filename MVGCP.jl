using LinearAlgebra, SparseArrays, ITensors

c6 = [
    10 5;
    5 10
]
rydberg_states = size(c6)[1]
node_locations = [0 0 ; 1 0]
detunings = [-2 -4]
rabis = [2 4]
initial_state = [x > 1 ? 0 : x for x in 1:((rydberg_states+1)^size(node_locations)[1])]

id = 1 * Matrix{Float64}(I, rydberg_states+1, rydberg_states+1)

function gen_drive(detunings, rabis)
    n = length(detunings)
    single = zeros(n+1, n+1)
    for i in 1:n
        single[i+1, i+1] = detunings[i]
        single[1, i+1] = rabis[i]
        single[i+1, 1] = rabis[i]
    end

    n_qudits = size(node_locations)[1]
    drive = zeros((rydberg_states+1)^n_qudits, (rydberg_states+1)^n_qudits)

    for i in 1:n_qudits
        left_product = id
        drive_temp = id
        i == 1 ? drive_temp = single : drive_temp = id
        for m in 2:n_qudits
            m == i ? left_product = single : left_product = id
            drive_temp = kron(left_product, drive_temp)
        end
        drive += drive_temp
    end
    return drive
end

function gen_interaction(c6, node_locations)
    n = size(node_locations)[1]
    dim = size(id)[1]^n
    h_int_total = zeros(dim, dim)

    for p in 1:rydberg_states
        for q in 1:rydberg_states
            for i in 1:n 
                for j in i+1:n
                    h_int_temp = id
                    left_product = id
                    int_i = zeros(rydberg_states+1, rydberg_states+1)
                    int_j = zeros(rydberg_states+1, rydberg_states+1)
                    distance = 1
                    int_i[p+1, p+1] = c6[p, q]/distance
                    int_j[q+1, q+1] = c6[p, q]/distance
                    i == 1 ? h_int_temp = int_i : h_int_temp = id
                    for m in 2:n
                        m == i ? left_product = int_i : m == j ? left_product = int_j : left_product = id
                        h_int_temp = kron(left_product, h_int_temp)
                    end
                    h_int_total += h_int_temp
                end
            end
        end
    end
    return h_int_total
end

gen_drive(detunings, rabis)

gen_interaction(c6, node_locations)
