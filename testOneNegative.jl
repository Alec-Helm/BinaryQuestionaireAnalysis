
#the version herein determines a Diff which is the average of the diffs for each individual relation

module testOneNegative

export nmIITA
function nmIITA(data)

    value = Array{Any, 1}(undef, 2)
    counterexamplesSurmise = determineSurmiseCounterexamples(data)
    counterexamplesExcludes = determineExcludesCounterexamples(data)

    relations = constructOrderRelations(data, counterexamplesSurmise, counterexamplesExcludes)


    value[1] = 0
    value[2] = modifiedDiff(counterexamplesSurmise, expectedCounterexamplesSurmise(data, relations[:,:,1], counterexamplesSurmise), counterexamplesExcludes, expectedCounterexamplesExcludes(data, relations[:,:,1], counterexamplesExcludes))

    for k in 2:size(relations)[3]

        if relations[:,:,k-1] != relations[:,:,k]
            newValue = modifiedDiff(counterexamplesSurmise, expectedCounterexamplesSurmise(data, relations[:,:,k], counterexamplesSurmise), counterexamplesExcludes, expectedCounterexamplesExcludes(data, relations[:,:,k], counterexamplesExcludes))

            if  newValue < value[2]
                value[1] = k-1
                value[2] = newValue
            end
        end

    end


    return value
end


export modifiedDiff
function modifiedDiff(positiveCounters, expectedPositiveCounters, negativeCounters, expectedNegativeCounters)
    value = (Diff(positiveCounters, expectedPositiveCounters) + Diff(negativeCounters, expectedNegativeCounters))/2

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
function constructOrderRelations(dataSet, positiveCounters, negativeCounters)
    numberQuestions = size(dataSet)[2]
    surmiseCounterexamples = positiveCounters
    excludeCounterexamples = negativeCounters

    orderRelations = Array{Int64, 3}(undef, numberQuestions, numberQuestions, 1)

    levelPlusOne  = 1
    for j in 1:numberQuestions
        for i in 1:numberQuestions
            if surmiseCounterexamples[i,j] == 0
                orderRelations[i,j, levelPlusOne] = 1
            elseif excludeCounterexamples[i,j] == 0
                orderRelations[i,j, levelPlusOne] = -1
            else
                orderRelations[i,j, levelPlusOne] = 0
            end
        end
    end

    while !isDone(surmiseCounterexamples, excludeCounterexamples, levelPlusOne)
        candidates = Array{Int64, 2}(undef, 0, 3)
        for i in 1:numberQuestions
            for j in 1:numberQuestions
                if (surmiseCounterexamples[i,j] <= levelPlusOne) &&  (orderRelations[i,j,levelPlusOne] == 0)
                    candidates = cat(dims = 1,  candidates[:,:], transpose([i,j,1]))
                end
                if (excludeCounterexamples[i,j] <= levelPlusOne) &&  (orderRelations[i,j,levelPlusOne] == 0)
                    candidates = cat(dims = 1,  candidates[:,:], transpose([i,j,-1]))
                end
            end
        end

        proposalRelation = orderRelations[:,:,levelPlusOne]
        for k in 1:size(candidates)[1]
            #handle case where in one step we try to both add (i,j) to Surmise and to Excludes, we include it in neither
            if proposalRelation[candidates[k,1], candidates[k,2]] != 0
                proposalRelation[candidates[k,1], candidates[k,2]] = 0
            else
                proposalRelation[candidates[k,1], candidates[k,2]] = candidates[k,3]
            end
        end

        while !isCompatable(proposalRelation)
            #remove all problematic candidates

            for k in 1:size(candidates)[1]

                if candidates[k,3] == 1

                    #handles intransitive triples
                    for a in 1:numberQuestions
                        if proposalRelation[a, candidates[k,1]] == 1 && proposalRelation[a,candidates[k,2]] != 1
                            candidates[k,:] = [1,1,1]
                        end
                        if proposalRelation[candidates[k,2], a] == 1 && proposalRelation[candidates[k,1], a] != 1
                            candidates[k,:] = [1,1,1]
                        end
                    end


                    #handles half of incompatability
                    for b in 1:numberQuestions
                        if proposalRelation[candidates[k,1], b] == -1 && proposalRelation[candidates[k,2], b] != -1
                            candidates[k,:] = [1,1,1]
                        end
                    end

                end


                if candidates[k,3] == -1

                    #handles failure of symmetry
                    if proposalRelation[candidates[k,2], candidates[k,1]] != -1
                        candidates[k,:] = [1,1,1]
                    end


                    #handles other half of incompatability
                    for c in 1:numberQuestions
                        if proposalRelation[candidates[k,1], c] == 1 && proposalRelation[candidates[k,2], c] != -1
                            candidates[k,:] = [1,1,1]
                        end
                        if proposalRelation[candidates[k,2], c] == 1 && proposalRelation[candidates[k,1], c] != -1
                            candidates[k,:] = [1,1,1]
                        end
                    end
                end
            end




            #reload our proposed relation using shortened candidate list
            proposalRelation = orderRelations[:,:,levelPlusOne]
            for k in 1:size(candidates)[1]
                proposalRelation[candidates[k,1], candidates[k,2]] = candidates[k,3]
            end
        end

        orderRelations = cat(dims = 3,  orderRelations[:,:,:], proposalRelation[:,:])

        levelPlusOne += 1
    end


    return orderRelations
end


#checked
export isDone
function isDone(surmiseCounterexamples, excludeCounterexamples, levelPlusOne)

    value = true

    if maximum(surmiseCounterexamples) > levelPlusOne + 1
        value = false
    elseif maximum(excludeCounterexamples) > levelPlusOne + 1
        value = false
    end

    return value

end

#checked
export isCompatable
function isCompatable(orderRelation)

    value = true

    for i in 1:size(orderRelation)[1]
        for j in 1:size(orderRelation)[2]

            if orderRelation[i,j] == 1

                for k in 1:size(orderRelation)[1]
                    if orderRelation[k,i] == 1 && orderRelation[k,j] != 1
                        value = false
                    end

                    if orderRelation[k,i] == -1 && orderRelation[k,j] != -1
                        value = false
                    end
                end
            end


            if orderRelation[i,j] == -1

                if orderRelation[j,i] != -1
                    value = false
                end
            end
        end
    end

    return value
end


export determineSurmiseCounterexamples
function determineSurmiseCounterexamples(dataSet)
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



export determineExcludesCounterexamples
function determineExcludesCounterexamples(dataSet)
    numQuestions = size(dataSet)[2]
    value = Array{Int64, 2}(undef, numQuestions, numQuestions)
    fill!(value, 0)

    for response in 1:size(dataSet)[1]
        for j in 1:numQuestions
            for i in 1:numQuestions
                if dataSet[response,i] == 1 && dataSet[response,j] == 1
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





export expectedCounterexamplesSurmise
function expectedCounterexamplesSurmise(data, orderRelation, counterexamples)
    numQs = size(data)[2]
    numRs = size(data)[1]
    value = Array{Float64, 2}(undef, numQs, numQs)

    #determine error rate
    gamma = 0
    denom = 0
    for i in 1:numQs
        for j in 1:numQs
            if orderRelation[i,j] == 1 && i != j
                gamma += counterexamples[i,j] / (sum(data[:,j]))
                denom += 1
            end
        end
    end
    if gamma != 0
        gamma /= denom
    end


    for i in 1:numQs
        for j in 1:numQs

            if orderRelation[i,j] == 1
                value[i,j] = (sum(data[:,j])/numRs) * gamma * numRs
            end


            if orderRelation[i,j] != 1
                if orderRelation[j,i] != 1
                    value[i,j] = (1-(sum(data[:,i])/numRs)) * (sum(data[:,j])/numRs) * numRs
                end
                if orderRelation[j,i] == 1
                    value[i,j] = ((sum(data[:,j])/numRs) - (sum(data[:,i])/numRs) + (gamma * (sum(data[:,i])/numRs))) * numRs
                end

            end

        end
    end

    return value
end




export expectedCounterexamplesExcludes
function expectedCounterexamplesExcludes(data, orderRelation, counterexamples)
    numQs = size(data)[2]
    numRs = size(data)[1]
    value = Array{Float64, 2}(undef, numQs, numQs)

    #determine error rate
    gammaNumerator = 0
    gammaDenominator = 0
    for i in 1:numQs
        for j in 1:numQs
            if orderRelation[i,j] == -1

                #find number of people that knew either
                denom = 0
                for k in 1:numRs
                    if data[k,i] == 1 || data[k,j] == 1
                        denom += 1
                    end
                end
                gammaNumerator += counterexamples[i,j] / (denom)
                gammaDenominator += 1
            end
        end
    end
    if gammaDenominator != 0
        gamma = gammaNumerator / gammaDenominator
    else
        gamma = 0
    end


    for i in 1:numQs
        for j in 1:numQs

            if orderRelation[i,j] == -1
                value[i,j] =  sum(data[:,i])* gamma + sum(data[:,j]) * gamma
            end


            if orderRelation[i,j] != -1
                value[i,j] = (sum(data[:,i])/numRs) * (sum(data[:,j])/numRs) * numRs * (1 - gamma)
            end

        end
    end

    return value
end






end
