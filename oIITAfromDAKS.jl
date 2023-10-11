module oIITAfromDAKS

export performIITA
function performIITA(dataSet)
    println("getting started...")

    #construct initial variables
    numberQuestions = size(dataSet)[2]
    numberResponses = size(dataSet)[1]

    println("determining counterexamples...")


    #construct b_i_j matrix
    counterexamples = determineCounterexamples(dataSet)


    println("constructing order heirarchy....")
    #construct order heirarchy
    constructOrderRelations(dataSet)

    println("done constructing order heirarchy")


    numberRelations = size(orderRelations)[3]

    diffValues = zeros(numberRelations)
    println("Calculating diff...")

    for l in 1:numberRelations

        expectedCounterExamples = expectedCounters(dataSet, orderRelations[:,:,l], counterexamples)

        for i in 1:numberQuestions
            for j in 1:numberQuestions
                diffValues[l] += (counterexamples[i,j] - expectedCounterExamples[i,j])^2
            end
        end

        diffValues[l] /= (numberQuestions * (numberQuestions -1))

    end


    return diffValues
end


export constructOrderRelations
function constructOrderRelations(dataSet)
    numberQuestions = size(dataSet)[2]
    counterexamples = determineCounterexamples(dataSet)
    println("counterexamples constructed")
    orderRelations = Array{Bool, 3}(undef, numberQuestions, numberQuestions, 1)
    level  = 1
    for j in 1:numberQuestions
        for i in 1:numberQuestions
            if counterexamples[i,j] == 0
                orderRelations[i,j, level] = true
            end
        end
    end

    println("level one complete")

    while !isTotal(orderRelations[:,:,level])
        println("currently at level ", level)

        level += 1

        candidates = Array{Int64, 2}(undef, 0, 2)
        for i in 1:numberQuestions
            for j in 1:numberQuestions
                if (counterexamples[i,j] < level ) &&  (!orderRelations[i,j,level-1])
                    candidates = cat(dims = 1,  candidates, transpose([i,j]))
                end
            end
        end

        proposalRelation = orderRelations[:,:,level-1]
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

            proposalRelation = orderRelations[:,:,level-1]
            for k in 1:size(candidates)[1]
                proposalRelation[candidates[k,1], candidates[k,2]] = true
            end
        end

        orderRelations = cat(dims = 3,  orderRelations[:,:,:], proposalRelation[:,:])
    end


    return orderRelations
end



function isTotal(orderRelation)

    numEntries = size(orderRelation)[1] * size(orderRelation)[2]

    if sum(orderRelation) < numEntries
        return false
    else
        return true
    end
end


function isTransitive(orderRelation)

    value = true

    for i in 1:size(orderRelation)[1]
        for j in 1:size(orderRelation)[2]

            if orderRelation[i,j]

                for k in 1:size(orderRelation)[1]
                    if orderRelation[j,k] && !orderRelation[i,k]
                        value = false
                    end
                end
            end
        end
    end

    return value
end

function determineCounterexamples(dataSet)
    num = size(dataSet)[2]
    value = Array{Int64, 2}(undef, num, num)
    fill!(value, 0)

    for response in 1:size(dataSet)[1]
        for j in 1:num
            for i in 1:num
                if dataSet[response,i] == 0 && dataSet[response,j] == 1
                    value[i,j] += 1
                end
            end
        end
    end

    return value
end


function expectedCounters(data, orderRelation, counterexamples)
    numQs = size(data)[2]
    numRs = size(data)[1]
    value = Array{Int64, 2}(undef, numQs, numQs)

    gamma = 0
    for i in 1:numQs
        for j in 1:numQs
            if orderRelation[i,j] && i != j
                gamma += counterexamples[i,j] / (sum(data[:,j]) * numRs)
            end
        end
    end

    if gamma != 0
        gamma /= (sum(orderRelation) - numQs)
    end

    for i in 1:numQs
        for j in 1:numQs

            if orderRelation[i,j]

                value[i,j] = numRs * sum(data[:,j]) * gamma
            end

            if !orderRelation[i,j]
                value[i,j] = (1-sum(data[:,i])) * sum(data[:,j]) * numRs * (1- gamma)
            end

        end
    end

    return value
end




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

end
