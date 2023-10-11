module correctedIITA


export binaryPowerSet
function binaryPowerSet(cardinality)
    powerSet = Array{Array,1}(undef,2^cardinality)

    #load the power set with empty vectors
    powerSet = Array{Int64,2}(undef,2^cardinality, cardinality)

    #fill those in by binary counting rule
    for subset in 1:2^cardinality
        number = subset -1

        for b in cardinality:-1:1
            powerSet[subset,b] = number ÷ 2^(b-1)
            number = number%2^(b-1)
        end
    end


    return powerSet
end

export transitiveClosure
function transitiveClosure(orderRelation)
    length = size(orderRelation)[1]
    transitive = false

    while !transitive
        transitive = true
        for i in 1:length
            for j in 1:length

                #see if there is an intransitive trio and add the needed elements
                if orderRelation[i,j] == 1
                    for a in 1:length
                        if orderRelation[j,a] == 1 && orderRelation[i,a] == 0
                            transitive = false
                            orderRelation[i,a] = 1
                        end
                    end
                end

                if orderRelation[i,j] == 1
                    for b in 1:length
                        if orderRelation[b,i] == 1 && orderRelation[b,j] == 0
                            transitive = false
                            orderRelation[b,j] = 1
                        end
                    end
                end

            end
        end
    end

    return orderRelation
end

export expectedCounterExamples
function expectedCounterExamples(initial, terminal, counterexamples, orderRelation, responseFrequencies)
    value = 0

    if orderRelation[terminal, initial] == 1
        value = errorRate(counterexamples, orderRelation,responseFrequencies) * solutionFrequency(terminal, responseFrequencies) * sum(responseFrequencies[:,end])
    end

    if orderRelation[terminal, initial] == 0
        if orderRelation[initial, terminal] == 0
            value = (1 - solutionFrequency(initial, responseFrequencies)) * solutionFrequency(terminal, responseFrequencies) * sum(responseFrequencies[:,end])
        end
        if orderRelation[initial, terminal] == 1
            value = (solutionFrequency(terminal, responseFrequencies) * sum(responseFrequencies[:,end])) - (solutionFrequency(initial, responseFrequencies) - solutionFrequency(initial, responseFrequencies) * errorRate(counterexamples, orderRelation,responseFrequencies))  * sum(responseFrequencies[:,end])
        end
    end


    return value
end

export solutionFrequency
function solutionFrequency(question, responseFrequencies)
    value = 0

    for response in 1:size(responseFrequencies)[1]
        if responseFrequencies[response,question] == 1
            value += responseFrequencies[response,end]
        end
    end

    return value
end


export errorRate
function errorRate(counterexamples, orderRelation,responseFrequencies)
    value = 0
    m = sum(responseFrequencies[:,end])

    for i in 1:size(orderRelation)[1]
        for j in 1:size(orderRelation)[1]
            if i != j && orderRelation[j,i] == 1
                value += (counterexamples[j,i]/ (solutionFrequency(j, responseFrequencies)*m))
            end
        end
    end


    denominator = sum(orderRelation) - size(orderRelation)[1]

    if denominator != 0
        value /= denominator
    end

    return value
end


export computeDiff
function computeDiff(orderRelation, counterExamples, responseFrequencies)
    value = 0

    n = size(orderRelation)[1]

    for i in 1:n
        for j in 1:n
            if i != j
                value += (counterExamples[j,i] - expectedCounterExamples(i, j, counterExamples, orderRelation, responseFrequencies)) ^ 2 / (n * (n-1))
            end
        end
    end



    return value
end

export diffGammeCoeffThingy
function diffGammeCoeffThingy(orderRelation, counterExamples, responseFrequencies)
    value1 = 0
    value2 = 0
    value3 = 0
    value4 = 0
    m = sum(responseFrequencies[:,end])

    for i in 1:size(orderRelation)[1]
        for j in 1:size(orderRelation)[2]

            #x1 case
            if orderRelation[i,j] == 0 && orderRelation[j,i] == 1
                value1 += (-2)*(counterExamples[j,i]) * solutionFrequency(i, responseFrequencies)*(m) + 2*solutionFrequency(i, responseFrequencies)*solutionFrequency(j, responseFrequencies)*m*m - 2*(solutionFrequency(i, responseFrequencies)^2) * m*m
                value3 += errorRate(counterExamples, orderRelation) * 2 * (solutionFrequency(i, responseFrequencies)^2) * m * m
            end
            if orderRelation[i,j] == 1
                value2 += (-2) * (counterExamples[j,i]) * solutionFrequency(j, responseFrequencies) * m
                value4 += errorRate(counterExamples, orderRelation) * 2 * (solutionFrequency(i, responseFrequencies)^2) * m * m
            end
        end
    end


    value = -(value1 + value2)/(value3 + value4)

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
        level = level + 1
        done = isTotal(value[:,:,level])

    end


    return value
end

function isTotal(orderRelation)
    value = true
    for i in 1:size(orderRelation)[1], j in 1:size(orderRelation)[2]
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
    for i in 1:size(counterexamples)[1]
        for j in 1:size(counterexamples)[2]
            if counterexamples[i,j] <= level && orderRelation[i,j] == 0
                candidates = cat(dims = 1, candidates, transpose([i,j]))
            end
        end
    end

    stop = false

    #load A_l+1
    bigCheck = orderRelation[:,:]
    for k in 1:size(candidates)[1]
        bigCheck[candidates[k,1], candidates[k,2]] = 1
    end


while !stop
    upForRemoval = Array{Int64,2}(undef,0,2)

    for k in 1:size(candidates)[1]

        #look for a->i && ~(a->j) case
        for a in 1:size(bigCheck)[2]
            if bigCheck[candidates[k,2],a] == 1 && bigCheck[candidates[k,1],a] == 0

                upForRemoval = cat(dims = 1, upForRemoval, transpose(candidates[k,:]))

            end
        end


        #look for j ->b and ~(i->b)
        for b in 1:size(bigCheck)[2]
            if bigCheck[b,candidates[k,1]] == 1 && bigCheck[b,candidates[k,2]] == 0

                upForRemoval = cat(dims = 1, upForRemoval, transpose(candidates[k,:]))

            end
        end
    end


    #remove what's up for removal
    for c in 1:size(upForRemoval)[1]
        bigCheck[upForRemoval[c,1],upForRemoval[c,2]] = 0
    end

    stop = isTransitive(bigCheck)
end



    return bigCheck
end


export isTransitive
function isTransitive(orderRelation)
    value = true

    for a in 1:size(orderRelation)[1]
        for b in 1:size(orderRelation)[1]
            #see if (a,b) is in the relation, if it is see if it is involved in an intransitive triple
            if orderRelation[a,b] == 1
                for i in 1:size(orderRelation)[1]
                    if orderRelation[i,a] == 1 && orderRelation[i,b] == 0
                        value = false
                    end
                end

                for j in 1:size(orderRelation)[1]
                    if orderRelation[b,j] == 1 && orderRelation[a,j] == 0
                        value = false
                    end
                end
            end
        end
    end

    return value
end

export computeSymmetricDifference
function computeSymmetricDifference(order1, order2)
    value = 0
    num = size(order1)[1]

    for i in 1:num
        for j in 1:num
            if (order1[i,j] == 1 && order2[i,j] == 0) || (order2[i,j] == 1 && order1[i,j] == 0)
                value +=1
            end
        end
    end

    return value
end




end
