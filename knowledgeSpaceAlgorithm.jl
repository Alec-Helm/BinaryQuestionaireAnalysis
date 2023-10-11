#example 1 from paper
problemSet = [1,2,3,4,5]
knowledgeSpace = Array{Array,1}(undef,8)
knowledgeSpace[1] = [0,0,0,0,0]
knowledgeSpace[2] = [0,0,0,0,1]
knowledgeSpace[3] = [0,0,1,0,0]
knowledgeSpace[4] = [0,0,0,1,1]
knowledgeSpace[5] = [0,0,1,0,1]
knowledgeSpace[6] = [0,0,1,1,1]
knowledgeSpace[7] = [0,1,1,1,1]
knowledgeSpace[8] = [1,1,1,1,1]
println("output: knowledgeSpace")
