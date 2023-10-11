module ideaForNegativeImplication



#checked
export constructOrderRelations
function constructOrderRelations(dataSet, positiveCounters, negativeCounters)
    numberQuestions = size(dataSet)[2]

    orderRelations = Array{Bool, 4}(undef, numberQuestions, numberQuestions, 1, 2)

    levelPlusOne  = 1
    for j in 1:numberQuestions
        for i in 1:numberQuestions
            if positiveCounters[i,j] == 0
                orderRelations[i,j, levelPlusOne,1] = true
            end
            if positiveCounters[i,j] != 0
                orderRelations[i,j, levelPlusOne,1] = false
            end

            if negativeCounters[i,j] == 0
                orderRelations[i,j, levelPlusOne,2] = true
            end
            if negativeCounters[i,j] != 0
                orderRelations[i,j, levelPlusOne,2] = false
            end
        end
    end


    while !isTotal(orderRelations[:,:,levelPlusOne,:])

        candidates = Array{Int64, 2}(undef, 0, 3)
        for i in 1:numberQuestions
            for j in 1:numberQuestions
                if (positiveCounters[i,j] <= levelPlusOne) &&  (!orderRelations[i,j,levelPlusOne,1])
                    candidates = cat(dims = 1,  candidates[:,:], transpose([i,j,1]))
                end
                if (negativeCounters[i,j] <= levelPlusOne) &&  (!orderRelations[i,j,levelPlusOne,2])
                    candidates = cat(dims = 1,  candidates[:,:], transpose([i,j,2]))
                end
            end
        end


        proposalRelation = orderRelations[:,:,levelPlusOne,:]
        for k in 1:size(candidates)[1]
            proposalRelation[candidates[k,1], candidates[k,2],candidates[k,3]] = true
        end

        while !satisfiesProperties(proposalRelation)
            #remove all problematic candidates

            for k in 1:size(candidates)[1]
                #remove intransitive triples

                if candidates[k,3] == 1
                    for b in 1:numberQuestions
                        if proposalRelation[candidates[k,2], b,1] && !proposalRelation[candidates[k,1], b,1]
                            candidates[k,:] = [1,1,1]
                        end
                    end
                    for a in 1:numberQuestions
                        if proposalRelation[a, candidates[k,1],1] && !proposalRelation[a,candidates[k,2],1]
                            candidates[k,:] = [1,1,1]
                        end
                    end
                end

                #remove incompatabilities
                for b in 1:numberQuestions
                    if proposalRelation[candidates[k,2], b,2] && !proposalRelation[candidates[k,1], b,2]
                        candidates[k,:] = [1,1,1]
                    end
                end





                if candidates[k,3] == 2
                    #remove non-symmetric pairs
                    for i in 1:size(candidates)[1]
                        if !proposalRelation[candidates[k,2], candidates[k,1], 2]
                            candidates[k,:] = [1,1,1]
                        end
                    end


                    #remove imcompatabilities
                    for a in 1:numberQuestions
                        if proposalRelation[candidates[k,1],a,1] && !proposalRelation[a,candidates[k,2],2]
                            candidates[k,:] = [1,1,1]
                        end
                        if proposalRelation[candidates[k,2],a,1] && !proposalRelation[a,candidates[k,1],2]
                            candidates[k,:] = [1,1,1]
                        end
                    end

                end
            end



            proposalRelation = orderRelations[:,:,levelPlusOne]
            for k in 1:size(candidates)[1]
                proposalRelation[candidates[k,1], candidates[k,2]] = true
            end
        end

        orderRelations = cat(dims = 3,  orderRelations[:,:,:], proposalRelation[:,:], :)

        levelPlusOne += 1
    end


    return orderRelations
end


export satisfiesProperties
function satisfiesProperties(orderRelation)

    value = true

    #make sure < is transitive
    for i in 1:size(orderRelation)[1]
        for j in 1:size(orderRelation)[2]

            if orderRelation[i,j,1]

                for k in 1:size(orderRelation)[1]
                    if orderRelation[k,i,1] && !orderRelation[k,j,1]
                        value = false
                    end
                    if orderRelation[j,k,1] && !orderRelation[i,k,1]
                        value = false
                    end
                end
            end
        end
    end


    #make sure X is symmetric
    for i in 1:size(orderRelation)[1]
        for j in 1:size(orderRelation)[2]

            if orderRelation[i,j,2]

                if !orderRelation[j,i,2]
                    value = false
                end

            end
        end
    end



    #check compatability
    for i in 1:size(orderRelation)[1]
        for j in 1:size(orderRelation)[1]

            if orderRelation[i,j,1]

                for k in 1:size(orderRelation)[1]
                    if orderRelation[k,i,2] && !orderRelation[k,j,2]
                        value = false
                    end
                end
            end
        end
    end


    return value
end



export isTotal
function isTotal(orderRelation)
    value = true
    number = size(orderRelation)[1]

    for i in 1:number
        for j in 1:number
            if !orderRelation[i,j,1] && !orderRelation[i,j,2]
                value = false
            end
        end
    end


    return value
end


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


export determinePositiveCounterexamples
function determinePositiveCounterexamples(dataSet)
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
export determineNegativeCounterexamples
function determineNegativeCounterexamples(dataSet)
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



export relationsClosure
function relationsClosure(orderRelation)
    length = size(orderRelation)[1]
    closed = false
    newRelation = orderRelation[:,:,:]

    while !closed
        closed = true
        #make sure < is transitive
        for i in 1:size(newRelation)[1]
            for j in 1:size(newRelation)[2]

                if newRelation[i,j,1]

                    for k in 1:size(newRelation)[1]
                        if newRelation[k,i,1] && !newRelation[k,j,1]
                            newRelation[k,j,1] = true
                            closed = false
                        end
                        if newRelation[j,k,1] && !newRelation[i,k,1]
                            newRelation[i,k,1] = true
                            closed = false
                        end
                    end
                end
            end
        end


        #make sure X is symmetric
        for i in 1:size(newRelation)[1]
            for j in 1:size(newRelation)[2]

                if newRelation[i,j,2]

                    if !newRelation[j,i,2]
                        newRelation[j,i,2] = true
                        closed = false
                    end

                end
            end
        end



        #check compatability
        for i in 1:size(newRelation)[1]
            for j in 1:size(newRelation)[1]

                if newRelation[i,j,1]

                    for k in 1:size(newRelation)[1]
                        if newRelation[k,i,2] && !newRelation[k,j,2]
                            newRelation[k,j,2] = true
                            closed = false
                        end
                    end
                end
            end
        end
    end


    return newRelation
end


end
