#make function to compute barcodes in dimension 0-2 for given connectivity matrix
module cElegansFunctions

function computeBarcodes(matrix, timeList)
    values = zeros(279,3)
    eireneMatrix = eirene(matrix, model = "vr", minrad = 0.5, maxrad = 3001, maxdim = 2)

    barcodes = barcode(eireneMatrix, dim=0)
    amount = size(barcodes,1)
    bettiNumber = 0


    #calculate and store dimension 0 betti values
    for i in 1:279
        for j in 1:amount
            if barcodes[j,1] == timeList[i]
                bettiNumber += 1
            end
            if barcodes[j,2] == timeList[i]
                bettiNumber -= 1
            end
        end
        values[i,1,k] = bettiNumber
    end

    #calculate and store dimension 1 betti values (resetting variables first)
    barcodes = barcode(eireneMatrix, dim=1)
    amount = size(barcodes,1)
    bettiNumber = 0

    for i in 1:279
        for j in 1:amount
            if barcodes[j,1] == timeList[i]
                bettiNumber += 1
            end
            if barcodes[j,2] == timeList[i]
                bettiNumber -= 1
            end
        end
        values[i,2,k] = bettiNumber
    end

    #calculate and store dimension 2 betti values (resetting variables first)
    barcodes = barcode(eireneMatrix, dim=2)
    amount = size(barcodes,1)
    bettiNumber = 0

    for i in 1:279
        for j in 1:amount
            if barcodes[j,1] == timeList[i]
                bettiNumber += 1
            end
            if barcodes[j,2] == timeList[i]
                bettiNumber -= 1
            end
        end
        values[i,3,k] = bettiNumber
    end

    return values
end

end



function computeRates(matrix, timeList)
    values = zeros(279,3,2)
    eireneMatrix = eirene(matrix, model = "vr", minrad = 0.5, maxrad = 3001, maxdim = 2)

    barcodes = barcode(eireneMatrix, dim=0)
    amount = size(barcodes,1)


    #calculate and store dimension 0 betti values
    for i in 1:279
        for j in 1:amount
            if barcodes[j,1] == timeList[i]
                values[i,1,1] += 1
            end
            if barcodes[j,2] == timeList[i]
                values[i,1,2] += 1
            end
        end
    end

    #calculate and store dimension 1 betti values (resetting variables first)
    barcodes = barcode(eireneMatrix, dim=1)
    amount = size(barcodes,1)

    for i in 1:279
        for j in 1:amount
            if barcodes[j,1] == timeList[i]
                values[i,2,1] += 1
            end
            if barcodes[j,2] == timeList[i]
                values[i,2,2] += 1
            end
        end
    end

    #calculate and store dimension 2 betti values (resetting variables first)
    barcodes = barcode(eireneMatrix, dim=2)
    amount = size(barcodes,1)

    for i in 1:279
        for j in 1:amount
            if barcodes[j,1] == timeList[i]
                values[i,3,1] += 1
            end
            if barcodes[j,2] == timeList[i]
                values[i,3,2] += 1
            end
        end
    end

    return values
end

end
