module binaryQuestionaireFunctions


export binaryPowerSet
function binaryPowerSet(set)
    cardinality = size(set)[1]
    powerSet = Array{Array,1}(undef,2^cardinality)

    #load the power set with empty vectors
    for subset in 1:2^cardinality
        powerSet[subset] = Array{Int64,1}(undef,cardinality)
    end

    #fill those in by binary counting rule
    for subset in 1:2^cardinality
        number = subset -1

        for b in cardinality:-1:1
            (powerSet[subset])[b] = number ÷ 2^(b-1)
            number = number%2^(b-1)
        end
    end


    return powerSet
end




export b_i_j
function b_i_j(frequencySet)
    value = Array{Int64,2}(undef,size(frequencySet)[2]-1,size(frequencySet)[2]-1)
    fill!(value, 0)

    for k in 1:size(frequencySet)[1]
        for i in 1:size(frequencySet)[2]-1
            for j in 1:size(frequencySet)[2]-1
                if frequencySet[k,i] == 1 && frequencySet[k,j] == 0
                    value[i,j] += frequencySet[k,end]
                end
            end
        end

    end

    return value
end

export generateOrderHeirarchy
function generateOrderHeirarchy(counterexamples)
    #generate A_0 through A_m given bij matrix
    value = Array{Int64,3}(undef,size(counterexamples)[1],size(counterexamples)[2],1)
    fill!(value, 0)
    for i in 1:size(counterexamples)[1]
        for j in 1:size(counterexamples)[2]
            if counterexamples[i,j] == 0
                value[i,j,1] = 1
            end
        end
    end

    level = 1
    done = isTotal(value[:,:,level])
    while !done
        new = A_L_1(value[:,:,level], counterexamples, level)
        value = cat(dims = 3, value[:,:,:], new[:,:])
        level += 1
        done = isTotal(value[:,:,level])
    end


    return value
end

function isTotal(orderRelation)
    value = true
    for i in 1:size(orderRelation)[1], j in 1:size(orderRelation)[1]
        if orderRelation[i,j] == 0
            value = false
            break
        end
    end

    return value
end




function A_L_1(orderRelation, counterexamples, level)
    #generate A_L+1 given A_L
    value = orderRelation[:,:]

    #construct list of possible additions
    candidates = Array{Int64,2}(undef,0,2)
    for i in 1:size(orderRelation)[1]
        for j in 1:size(orderRelation)[2]
            if counterexamples[i,j] <= level && orderRelation[i,j] == 0
                candidates = cat(dims = 1, candidates, transpose([i,j]))
            end
        end
    end


    #remove items that cause intransitivity
    bigCheckAgainstRelation = orderRelation[:,:]
    for k in 1:size(candidates)[1]
        bigCheckAgainstRelation[candidates[k,1],candidates[k,2]] = 1
    end

    if level > 58
        println()
        print(level, " : ")
    end


    for k in 1:size(candidates)[1]

        #look for a->i && ~(a->j) case
        for a in 1:size(bigCheckAgainstRelation)[2]
            if bigCheckAgainstRelation[a,candidates[k,1]] == 1 && bigCheckAgainstRelation[a,candidates[k,2]] == 0
                if level > 58
                    print("a: ",a," ", candidates[k,:], " ")
                end
                candidates[k,:] = [1,1]
                break
            end
        end


        #look for j ->b and ~(i->b)
        for b in 1:size(bigCheckAgainstRelation)[2]
            if bigCheckAgainstRelation[candidates[k,2],b] == 1 && bigCheckAgainstRelation[candidates[k,1],b] == 0
                if level > 58
                    print("b: ",b," ", candidates[k,:], " ")
                end
                candidates[k,:] = [1,1]
                break
            end
        end
    end




    #add all the maintained additions to the order relation
    for k in 1:size(candidates)[1]
        value[candidates[k,1], candidates[k,2]] = 1
    end


    return value
end


export p_i
function p_i(frequencySet, i)
    N = sum(frequencySet[:,end])
    value = 0

    for k in 1:size(frequencySet)[1]
        if frequencySet[k,i] == 1
            value += frequencySet[k,end]
        end
    end

    value = value/N

    return value
end


export o_i
function o_i(frequencySet, i)
    N = sum(frequencySet[:,end])
    value = 0

    for k in 1:size(frequencySet)[1]
        if frequencySet[k,i] == 2
            value += frequencySet[k,end]
        end
    end

    value = value/N

    return value
end

export q_i_j
function q_i_j(frequencySet, i, j)
    N = sum(frequencySet[:,end])
    value = 0

    for k in 1:size(frequencySet)[1]
        if frequencySet[k,i] != 2 && frequencySet[k,j] == 1
            value += frequencySet[k,end]
        end
    end

    value = value/N
    return value
end


export gamma
function gamma(orderRelation, frequencySet, distances)
    value = 0


    for i in 1:size(orderRelation)[1]
        for j in 1:size(orderRelation)[1]
            if orderRelation[i,j] == 1 && i != j
                value += (distances[i,j]/q_i_j(frequencySet,i,j))
            end
        end
    end

    denominator = sum(orderRelation) - size(orderRelation)[1]

    value = value / denominator

    value = value/200

    return value
end

export t_i_j
function t_i_j(frequencySet, i, j, orderRelation, distances)
    value = 0

    if orderRelation[i,j] == 1
        value = q_i_j(frequencySet, i, j) * sum(frequencySet[:,end]) * gamma(orderRelation, frequencySet, distances)
    end
    if orderRelation[i,j] == 0
        value = (1 - (p_i(frequencySet,i) + o_i(frequencySet, i))) * p_i(frequencySet,j) * sum(frequencySet[:,end])
    end

    return value
end


export Diff
function Diff(orderRelation, frequencySet, counterexamples)
    value = 0
    length = size(orderRelation)[1]

    for i in 1:length
        for j in 1:length
            if i != j
                value += (counterexamples[i,j] - t_i_j(frequencySet, i, j, orderRelation, counterexamples)) ^2
            end
        end
    end

    value = value / (length ^2 - length)

    return value
end

end
