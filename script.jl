include("functions.jl");


##successful sample run 
config = [(1,1), (2,3), (3,2)]
n = 3;
lambda = 3; 
point = (1,1);
runCollisionChain(n, point, lambda, "right", config, printName = "samplerun", printBool = false);


###
#all possible configurations for n = 2, k = 2
generateConfiguration(2,2,true, "allpossiblen=2k=2")

##
#all possible configurations for n = 3, k = 3
generateConfiguration(3,3, true, "allpossiblen=3k=3")



# configurations and moves which result in invalid chains
# when dist = 1 and lambda is fixed 
#when n = 4 
n = 4
config1 = [(1,1), (2,2), (3,3), (4,4)];
point1 = (4,4);
lambda1 = 3; 
runCollisionChain(n, point1, lambda1, "right", config1, printName = "failure1", printBool = true, lambdaFixed = true);


point2 = (1,1);
lambda2 = 4; 
config2 = [(1,1), (2,2), (3,1), (4,2)];
runCollisionChain(n, point2, lambda2, "up", config1, printName = "failure2", printBool = true, lambdaFixed = true);


#when n = 5
n = 5
config1 = [(1,1), (2,2), (3,3), (4,4), (5,5)];
point1 = (5,5);
lambda1 = 3; 
runCollisionChain(n, point1, lambda1, "right", config1, printName = "failure1", printBool = true, lambdaFixed = true);


config2 = [(1,1), (2,2), (3,1), (4,2), (5,3)];
point2 = (1,1);
lambda2 = 5; 
runCollisionChain(n, point2, lambda2, "up", config2, printName = "failure2", printBool = true, lambdaFixed = true);



#check magnitude of second eigenvalue against vanilla metropolis
for n in 2:5
	for k in 2:(n*n-1) 
		K1 = generateTransitionMatrix(n,k,lambda)[1];
		e1 = eig(K1)[1];
		K2 = generateTransitionMatrixMetropolis(n,k);
		e2 = eig(K2)[1];
		println("n:", n, " k:",k, " RFM: ", norm(e1[2]), " Vanilla: ", e2[(length(e2)-1)]);
	end
end












