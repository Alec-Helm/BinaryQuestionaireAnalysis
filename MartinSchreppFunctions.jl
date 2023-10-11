module martinSchreppFunctions

export checkClosure
function checkClosure(knowledgeSpace) #knowledge space is a 1xn matrix of px1 vectors
    closed = true
    cardinality = size(knowledgeSpace)[1]
    numberOfQuestions = size(knowledgeSpace[1])[1]

    #check for closure under union
    for i in 2:cardinality
        for j in 1:i-1

            #construct the union of the two
            union = zeros(numberOfQuestions)
            for p in 1:numberOfQuestions
                union[p] = trunc(Int,maximum([(knowledgeSpace[i])[p], (knowledgeSpace[j])[p]]))
            end

            #see if the union belongs to the set
            found = false

            for k in 1:cardinality
                if knowledgeSpace[k] == union
                    found = true
                end
            end

            #if it did not, update accordingly
            if !found
                closed = false
            end
        end
    end

    #check for closure under intersection
    for i in 2:cardinality
        for j in 1:i-1

            #construct the intersection of the two
            intersection = zeros(numberOfQuestions)
            for p in 1:numberOfQuestions
                intersection[p] = trunc(Int,minimum([(knowledgeSpace[i])[p], (knowledgeSpace[j])[p]]))
            end

            #see if the union belongs to the set
            found = false

            for k in 1:cardinality
                if knowledgeSpace[k] == intersection
                    found = true
                end
            end

            #if it did not, update accordingly
            if !found
                closed = false
            end
        end
    end

    return closed
end


export closeSet
function closeSet(knowledgeSpace) #knowledge space is a 1xn matrix of px1 vectors
done = false
while !done

    startOver = false
    cardinality = size(knowledgeSpace)[1]
    numberOfQuestions = size(knowledgeSpace[1])[1]

    #check for closure under union
    for i in 2:cardinality
        for j in 1:i-1

            #construct the union of the two
            union = Array{Int64,1}(undef, numberOfQuestions)
            for p in 1:numberOfQuestions
                union[p] = maximum([(knowledgeSpace[i])[p], (knowledgeSpace[j])[p]])
            end

            #see if the union belongs to the set
            found = false

            for k in 1:cardinality
                if knowledgeSpace[k] == union
                    found = true
                end
            end

            #if it did not, add it to the set and make it run through again
            if !found && !startOver
                knowledgeSpace = [knowledgeSpace;[union]]
                startOver = true
            end
        end
    end

    #check for closure under intersection
    for i in 2:cardinality
        for j in 1:i-1

            #construct the intersection of the two
            intersection = Array{Int64,1}(undef, numberOfQuestions)
            for p in 1:numberOfQuestions
                intersection[p] = minimum([(knowledgeSpace[i])[p], (knowledgeSpace[j])[p]])
            end

            #see if the union belongs to the set
            found = false

            for k in 1:cardinality
                if knowledgeSpace[k] == intersection
                    found = true
                end
            end

            #if it did not, then add it and start over
            if !found && !startOver
                knowledgeSpace = [knowledgeSpace;[intersection]]
                startOver = true
            end
        end
    end
    if !startOver
        done = true
    end
end

    return knowledgeSpace
end


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


export binaryToLetters
function binaryToLetters(array)
    alpha = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"]

    hold = ""

    for index in 1:size(array)[1]
        if array[index] == 1
            hold = string(hold, alpha[index])
        end
    end

    letter = "{"

    if length(hold) > 0
        letter = string(letter, hold[1])
    end


    for index in 2:length(hold)
        letter = string(letter, ",",hold[index])
    end

    letter = string(letter,"}")

    return letter
end


export K_L
function K_L(frequencySet, threshold) #frequencySet is nx2, a column of all subsets and a column of frequencies observed, threshold is frequency threshold
    leveled = zeros(0,size(frequencySet)[2])

    for k in 1:size(frequencySet)[1]
        if frequencySet[k,end] >= threshold
            leveled = [leveled;transpose(frequencySet[k,:])]
        end
    end

    return leveled
end


export d_alpha
function d_alpha(response,knowledge)
    distance = 0
    numQuestions = maximum(size(response))
    for i in 1:numQuestions
        if knowledge[i] == 1 && response[i] == 0
            distance += 1
        end
    end

    return distance
end


export d_beta
function d_beta(response,knowledge)
    distance = 0
    numQuestions = maximum(size(response))

    for i in 1:numQuestions
        if response[i] == 1 && knowledge[i] == 0
            distance += 1
        end
    end

    return distance
end


export K_L_D
function K_L_D(frequencySet, response)
    numLevels = size(frequencySet)[1]
    currentBestDistance = 2*size(response)[1] + 1
    value = zeros(0,size(frequencySet)[2])

    for level in 1:numLevels
        if d_alpha(response,frequencySet[level,1:end-1]) + d_beta(response,frequencySet[level,1:end-1]) == currentBestDistance
            newRow = frequencySet[level,:]
            value = [value;newRow]
        end
        if d_alpha(response,frequencySet[level,1:end-1]) + d_beta(response,frequencySet[level,1:end-1]) < currentBestDistance
            currentBestDistance = d_alpha(response,frequencySet[level,1:end-1]) + d_beta(response,frequencySet[level,1:end-1])
            value = frequencySet[level,:]
        end
    end

    return value
end


export alpha_L_D
function alpha_L_D(levelSet, response)
    leveledSet = K_L_D(levelSet, response)

    levels = size(leveledSet)[1]
    errorFactor = 0

    if ndims(leveledSet) != 1
        for k in 1:levels
            errorFactor += d_alpha(response, leveledSet[k,1:end-1])
        end
        errorFactor /= levels
    end
    if ndims(leveledSet) == 1
        errorFactor += d_alpha(response, leveledSet)
    end



    return errorFactor
end



export beta_L_D
function beta_L_D(levelSet, response)
    leveledSet = K_L_D(levelSet, response)



    levels = size(leveledSet)[1]
    errorFactor = 0



    if ndims(leveledSet) != 1
        for k in 1:levels
            errorFactor += d_beta(response, leveledSet[k,1:end-1])
        end
        errorFactor /= levels
    end
    if ndims(leveledSet) == 1
        errorFactor += d_beta(response, leveledSet)
    end


    return errorFactor
end



export alpha_L_min
function alpha_L_min(levelSet, responseSet)
    numResponses = size(responseSet)[1]
    errorFactor = 0

    for k in 1:numResponses
        errorFactor +=  (alpha_L_D(levelSet, responseSet[k,1:end-1]) * responseSet[k,end])
    end

    errorFactor /= (maximum(size(responseSet[1,1:end-1])) * sum(responseSet[:,end]))

    return errorFactor
end



export beta_L_min
function beta_L_min(levelSet, responseSet)
    numResponses = size(responseSet)[1]
    errorFactor = 0

    for k in 1:numResponses
        errorFactor +=  (beta_L_D(levelSet, responseSet[k,1:end-1]) * responseSet[k,end])
    end

    errorFactor /= (maximum(size(responseSet[1,1:end-1])) * sum(responseSet[:,end]))

    return errorFactor
end



export P_k_d
function P_k_d(knowledgeSet ,response, levelSet, responseSet)
    KintD = 0
    QcomKD = 0

    for k in 1:size(response)[1]
        if knowledgeSet[k] == response[k] == 1
            KintD += 1
        end
        if knowledgeSet[k] == response[k] == 0
            QcomKD += 1
        end
    end

    value = ( (alpha_L_min(levelSet, responseSet)) ^ d_alpha(response,knowledgeSet)) * ( (beta_L_min(levelSet, responseSet)) ^ d_beta(response, knowledgeSet)) * ( (1 - alpha_L_min(levelSet, responseSet)) ^ KintD)  * ( (1 - beta_L_min(levelSet, responseSet)) ^ QcomKD)

    return value
end



export F_L_D
function F_L_D(response, levelSet, responseSet)
    levels = size(levelSet)[1]
    value = 0

    for k in 1:levels
        value += P_k_d(levelSet[k,1:end-1] ,response, levelSet, responseSet)
    end

    value /= levels

    return value
end



export app
function app(levelSet, frequencySet)
    index = size(frequencySet)[1]
    value = 0
    numResponses = sum(frequencySet[:,end])

    for k in 1:index
        value += ( frequencySet[k,end]/numResponses - F_L_D(frequencySet[k,1:end-1], levelSet, frequencySet)) ^ 2
    end

    value /= (index)

    return value
end




export calculateLopt
function calculateLopt(frequencySet)
    bestTolerance = 2
    Lopt = [-1,2]  #L x bestTolerance

    #find set of actually present tolerance levels
    maxFrequency = trunc(Int,maximum(frequencySet[:,end]))
    frequencyLevels = zeros(maxFrequency + 1)
    for level in 1:size(frequencySet)[1]
        index = trunc(Int,frequencySet[level,end]+1)
        frequencyLevels[index] = 1
    end

    #store array for the table
    number = trunc(Int,sum(frequencyLevels))
    toleranceSet = zeros(number)
    tableCounter = 0
    for level in 1:size(frequencyLevels)[1]
        if frequencyLevels[level] == 1
            tableCounter += 1
            toleranceSet[tableCounter] = level-1
        end
    end

    for level in toleranceSet
        levelSet = K_L(frequencySet,level)
        value = app(levelSet, frequencySet)

        if value == bestTolerance
            Lopt = [Lopt;[level,bestTolerance]]
        end
        if value < bestTolerance
            bestTolerance = value
            Lopt = [level,bestTolerance]
        end
    end

    return Lopt
end



export exampleFourSimulation
function exampleFourSimulation(numberQuestions, gamma, numTrials, alpha, beta) #numQuestions is an array of number of questions 1s, gamma is inclusion probability

    #randomly select a knowledge structure (it sounds like they ignored closure?)
    numberOfQuestions = maximum(size(numberQuestions))
    powerSet = binaryPowerSet(numberQuestions)
    knowledgeStructure = fill(0,1,numberOfQuestions)

    for i in 1:size(powerSet)[1]
        if rand() < gamma
            knowledgeStructure = [knowledgeStructure;transpose(powerSet[i][:])]
        end
    end

    #generate numTrials simulated responses by randomly choosing response patterns from K
    responseSet = zeros(numTrials,numberOfQuestions)
    upperBound = size(knowledgeStructure)[1]
    for subject in 1:numTrials
        knowledge = knowledgeStructure[rand(1:upperBound),:]
        for q in 1:numberOfQuestions
            if knowledge[q] == 1
                if rand()<alpha
                    knowledge[q] = 0
                end
            else
                if rand()< beta
                    knowledge[q] = 1
                end
            end
        end
        responseSet[subject,:] = knowledge[:]
    end

    #analyse D with algorithm
    frequencies = Array{Any, 2}(UndefInitializer(),2^numberOfQuestions,numberOfQuestions+1)
    for i in 1:2^numberOfQuestions
        frequencies[i,1:end-1] = powerSet[i]
    end
    frequencies[:,end] = zeros(2^numberOfQuestions)

    for subject in 1:size(responseSet)[1]
        for subset in 1:2^numberOfQuestions
            if responseSet[subject,:] == frequencies[subset,1:end-1]
                frequencies[subset,end] += 1
            end
        end
    end

    sortedFrequencies = frequencies[sortperm(frequencies[:, end]), :]



    Lopt = transpose(calculateLopt(sortedFrequencies))

    #determine diff(K,K_L)
    return diff(knowledgeStructure,K_L(sortedFrequencies, Lopt[1])[:,1:end-1])


end





export diff
function diff(trueStructure, optimumStructure)
    trueSize = size(trueStructure)[1]
    optimumSize = size(optimumStructure)[1]

    value = 0


    for i in 1:trueSize
        found = false

        for j in 1:optimumSize
            if trueStructure[i,:] == optimumStructure[j,:]
                found = true
            end
        end

        if found == false
            value += 1
        end
    end



    for i in 1:optimumSize
        found = false

        for j in 1:trueSize
            if optimumStructure[i,:] == trueStructure[j,:]
                found = true
            end
        end

        if found == false
            value += 1
        end
    end


    return value
end


export exampleSevenSimulation
function exampleSevenSimulation(numberQuestions, gamma, numTrials, alpha, beta) #numQuestions is an array of number of questions 1s, gamma is inclusion probability

    #randomly select a knowledge structure (it sounds like they ignored closure?)
    numberOfQuestions = maximum(size(numberQuestions))
    powerSet = binaryPowerSet(numberOfQuestions)
    knowledgeStructure = zeros(0,numberOfQuestions)

    for i in 1:size(powerSet)[1]
        if rand() < gamma
            knowledgeStructure = [knowledgeStructure;transpose(powerSet[i][:])]
        end
    end

    #find r(K) for each knowledgeState
    probabilities = zeros(size(knowledgeStructure)[1])
    for r in 1:size(knowledgeStructure)[1]
        probabilities(r) = rand(1:10)
    end
    norm = sum(probabilities)
    probabilities./norm

    #generate numTrials simulated responses by randomly choosing response patterns from K
    responseSet = zeros(numTrials,numberOfQuestions)
    upperBound = size(knowledgeStructure)[1]
    for subject in 1:numTrials
        roll = rand()
        choice = 0
        while roll > 0
            choice += 1
            roll -= probabilities(choice)
        end
        knowledge = knowledgeStructure[choice,:]
        for q in 1:numberOfQuestions
            if knowledge[q] == 1
                if rand()<alpha
                    knowledge[q] = 0
                end
            else
                if rand()< beta
                    knowledge[q] = 1
                end
            end
        end
        responseSet[subject,:] = knowledge[:]
    end
    #analyse D with algorithm
    frequencies = Array{Any, 2}(UndefInitializer(),2^numberOfQuestions,6)
    for i in 1:2^numberOfQuestions
        frequencies[i,1:end-1] = powerSet[i]
    end
    frequencies[:,end] = zeros(2^numberOfQuestions)

    for subject in 1:size(responseSet)[1]
        for subset in 1:2^numberOfQuestions
            if responseSet[subject,:] == frequencies[subset,1:end-1]
                frequencies[subset,end] += 1
            end
        end
    end

    sortedFrequencies = frequencies[sortperm(frequencies[:, end]), :]


    Lopt = transpose(calculateLopt(sortedFrequencies))

    #determine diff(K,K_L)
    return diff(knowledgeStructure,K_L(sortedFrequencies, Lopt[1]))


end


end
