module fromDAKS

export oIITA
function oIITA(data)

    value = Array{Any, 1}(undef, 2)
    counterexamples = determineCounterexamples(data)
    relations = constructOrderRelations(data, counterexamples)


    value[1] = 0
    value[2] = Diff(counterexamples, expectedCounterexamplesO(data, relations[:,:,1], counterexamples))

    for k in 2:size(relations)[3]

        if relations[:,:,k-1] != relations[:,:,k]
            newValue = Diff(counterexamples, expectedCounterexamplesO(data, relations[:,:,k], counterexamples))

            if  newValue < value[2]
                value[1] = k-1
                value[2] = newValue
            end
        end

    end


    return value
end


export cIITA
function cIITA(data)

    value = Array{Any, 1}(undef, 2)
    counterexamples = determineCounterexamples(data)
    relations = constructOrderRelations(data, counterexamples)


    value[1] = 0
    value[2] = Diff(counterexamples, expectedCounterexamplesC(data, relations[:,:,1], counterexamples))

    for k in 2:size(relations)[3]

        if relations[:,:,k-1] != relations[:,:,k]
            newValue = Diff(counterexamples, expectedCounterexamplesC(data, relations[:,:,k], counterexamples))

            if  newValue < value[2]
                value[1] = k-1
                value[2] = newValue
            end
        end

    end


    return value
end



export mcIITA
function mcIITA(data)

    value = Array{Any, 1}(undef, 2)
    counterexamples = determineCounterexamples(data)
    relations = constructOrderRelations(data, counterexamples)


    value[1] = 0
    value[2] = Diff(counterexamples, expectedCounterexamplesMC(data, relations[:,:,1], counterexamples))

    for k in 2:size(relations)[3]

        if relations[:,:,k-1] != relations[:,:,k]
            newValue = Diff(counterexamples, expectedCounterexamplesMC(data, relations[:,:,k], counterexamples))

            if  newValue < value[2]
                value[1] = k-1
                value[2] = newValue
            end
        end

    end


    return value
end






export Diff
function Diff(counters, expectedCounters)
    numberQuestions = size(counters)[1]

    value = 0

    for i in 1:numberQuestions
        for j in 1:numberQuestions
            if i != j
                value += (counters[i,j] - expectedCounters[i,j]) ^ 2
            end
        end
    end

    value /= (numberQuestions * (numberQuestions - 1))

    return value

end

#checked
export constructOrderRelations
function constructOrderRelations(dataSet, counters)
    numberQuestions = size(dataSet)[2]
    counterexamples = counters

    orderRelations = Array{Bool, 3}(undef, numberQuestions, numberQuestions, 1)
    levelPlusOne  = 1
    for j in 1:numberQuestions
        for i in 1:numberQuestions
            if counterexamples[i,j] == 0
                orderRelations[i,j, levelPlusOne] = true
            end
            if counterexamples[i,j] != 0
                orderRelations[i,j, levelPlusOne] = false
            end
        end
    end


    while !isTotal(orderRelations[:,:,levelPlusOne])

        candidates = Array{Int64, 2}(undef, 0, 2)
        for i in 1:numberQuestions
            for j in 1:numberQuestions
                if (counterexamples[i,j] <= levelPlusOne) &&  (!orderRelations[i,j,levelPlusOne])
                    candidates = cat(dims = 1,  candidates[:,:], transpose([i,j]))
                end
            end
        end

        proposalRelation = orderRelations[:,:,levelPlusOne]
        for k in 1:size(candidates)[1]
            proposalRelation[candidates[k,1], candidates[k,2]] = true
        end

        while !isTransitive(proposalRelation)
            #remove all problematic candidates
            for k in 1:size(candidates)[1]
                for b in 1:numberQuestions
                    if proposalRelation[candidates[k,2], b] && !proposalRelation[candidates[k,1], b]
                        candidates[k,:] = [1,1]
                    end
                end

                for a in 1:numberQuestions
                    if proposalRelation[a, candidates[k,1]] && !proposalRelation[a,candidates[k,2]]
                        candidates[k,:] = [1,1]
                    end
                end
            end

            proposalRelation = orderRelations[:,:,levelPlusOne]
            for k in 1:size(candidates)[1]
                proposalRelation[candidates[k,1], candidates[k,2]] = true
            end
        end

        orderRelations = cat(dims = 3,  orderRelations[:,:,:], proposalRelation[:,:])

        levelPlusOne += 1
    end


    return orderRelations
end


#checked
export isTotal
function isTotal(orderRelation)

    numEntries = size(orderRelation)[1] * size(orderRelation)[2]

    if sum(orderRelation) < numEntries
        return false
    else
        return true
    end
end

#checked
export isTransitive
function isTransitive(orderRelation)

    value = true

    for i in 1:size(orderRelation)[1]
        for j in 1:size(orderRelation)[2]

            if orderRelation[i,j]

                for k in 1:size(orderRelation)[1]
                    if orderRelation[k,i] && !orderRelation[k,j]
                        value = false
                    end
                end
            end
        end
    end

    return value
end

#checked
export determineCounterexamples
function determineCounterexamples(dataSet)
    numQuestions = size(dataSet)[2]
    value = Array{Int64, 2}(undef, numQuestions, numQuestions)
    fill!(value, 0)

    for response in 1:size(dataSet)[1]
        for j in 1:numQuestions
            for i in 1:numQuestions
                if dataSet[response,i] == 0 && dataSet[response,j] == 1
                    value[i,j] += 1
                end
            end
        end
    end

    return value
end



#checked
export binaryPowerSet
function binaryPowerSet(cardinality)
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

#checked
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


#
#functions of oIITA
#
export expectedCounterexamplesO
function expectedCounterexamplesO(data, orderRelation, counterexamples)
    numQs = size(data)[2]
    numRs = size(data)[1]
    value = Array{Float64, 2}(undef, numQs, numQs)

    #determine error rate
    gamma = 0
    for i in 1:numQs
        for j in 1:numQs
            if orderRelation[i,j] && i != j
                gamma += counterexamples[i,j] / ((sum(data[:,j])/numRs) * numRs)
            end
        end
    end
    if gamma != 0
        gamma /= (sum(orderRelation) - numQs)
    end


    for i in 1:numQs
        for j in 1:numQs

            if orderRelation[i,j]
                value[i,j] = numRs * (sum(data[:,j])/numRs) * gamma
            end

            if !orderRelation[i,j]
                value[i,j] = (1-(sum(data[:,i])/numRs)) * (sum(data[:,j])/numRs) * numRs * (1 - gamma)
            end

        end
    end

    return value
end

#
#functions of cIITA
#
export expectedCounterexamplesC
function expectedCounterexamplesC(data, orderRelation, counterexamples)
    numQs = size(data)[2]
    numRs = size(data)[1]
    value = Array{Float64, 2}(undef, numQs, numQs)

    #determine error rate
    gamma = 0
    for i in 1:numQs
        for j in 1:numQs
            if orderRelation[i,j] && i != j
                gamma += counterexamples[i,j] / ((sum(data[:,j])/numRs) * numRs)
            end
        end
    end
    if gamma != 0
        gamma /= (sum(orderRelation) - numQs)
    end


    for i in 1:numQs
        for j in 1:numQs

            if orderRelation[i,j]
                value[i,j] = (sum(data[:,j])/numRs) * gamma * numRs
            end


            if !orderRelation[i,j]
                if !orderRelation[j,i]
                    value[i,j] = (1-(sum(data[:,i])/numRs)) * (sum(data[:,j])/numRs) * numRs
                end
                if orderRelation[j,i]
                    value[i,j] = ((sum(data[:,j])/numRs) - (sum(data[:,i])/numRs) + (gamma * (sum(data[:,i])/numRs))) * numRs
                end

            end

        end
    end

    return value
end


#
#functions of mcIITA
#
export expectedCounterexamplesMC
function expectedCounterexamplesMC(data, orderRelation, counterexamples)
    numQs = size(data)[2]
    numRs = size(data)[1]
    value = Array{Float64, 2}(undef, numQs, numQs)

    #determine error rate
    gamma = gammaMC(counterexamples, data, orderRelation)


    for i in 1:numQs
        for j in 1:numQs

            if orderRelation[i,j]
                value[i,j] = (sum(data[:,j])/numRs) * gamma * numRs
            end


            if !orderRelation[i,j]
                if !orderRelation[j,i]
                    value[i,j] = (1-(sum(data[:,i])/numRs)) * (sum(data[:,j])/numRs) * numRs
                end
                if orderRelation[j,i]
                    value[i,j] = ((sum(data[:,j])/numRs) - (sum(data[:,i])/numRs) + (gamma * (sum(data[:,i])/numRs))) * numRs
                end

            end

        end
    end

    return value
end

export gammaMC
function gammaMC(counters, responses, orderRelation)
    numberQuestions = size(counters)[1]
    numRs = size(responses)[1]

    x = zeros(4)

    for i in 1: numberQuestions
        for j in 1:numberQuestions
            if !orderRelation[i,j] && orderRelation[j,i]
                x[1] += (-2 * counters[i,j] * (sum(responses[:,i])/numRs) * numRs + 2 * (sum(responses[:,i])/numRs) * (sum(responses[:,j])/numRs) * numRs * numRs - 2 * (sum(responses[:,i])/numRs) * (sum(responses[:,i])/numRs) * numRs * numRs)
                x[3] += 2 * (sum(responses[:,i])/numRs) * (sum(responses[:,i])/numRs) * numRs * numRs
            end

            if orderRelation[i,j]
                x[2] += -2 * counters[i,j] * (sum(responses[:,j])/numRs) * numRs
                x[4] += 2 * (sum(responses[:,j])/numRs) * (sum(responses[:,j])/numRs) * numRs * numRs
            end
        end
    end

    ErrorProbability = - ((x[1] + x[2]) / (x[3] + x[4]))

    return ErrorProbability

end

end
