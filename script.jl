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
		e1 = sort(abs(e1));
		K2 = generateTransitionMatrixMetropolis(n,k);
		e2 = eig(K2)[1];
		println("n:", n, " k:",k, " RFM: ", norm(e1[length(e1)-1]), " Vanilla: ", e2[(length(e2)-1)]);
	end
end


vm = [0.7500000000000017,
 0.8333333333333336, 
 0.8750000000000003, 
 0.9000000000000001,
 0.916666666666667,
 0.9285714285714286,
 0.9375000000000001,
 0.9444444444444444,
 0.9500000000000004,
 0.9545454545454546,
 0.9583333333333333,
 0.9615384615384617,
 0.9642857142857151,
 0.9666666666666669];
lambdarfm2 = [0.9114378277661477,
0.8204438753025938,
0.834864412946574,
0.8546846189506956,
0.8666435490085984, 
0.8739731905132966, 
0.8782704750934184, 
0.8804071586580813, 
0.8808953226700512, 
0.8800542414365945, 
0.878093522355144, 
0.8751558682725222,
0.8713350888225718, 
0.866666666666667];
lambdarfm1 = [0.9378772832046405,
0.8750766698089074,
0.8426858055293812,
0.8587955394617385,
0.8687299185773811, 
0.8748238119471867, 
0.8783086092008078, 
0.879877376220958, 
0.8799400218374752, 
0.8787545023948367,
0.8765119837729034,
0.873427566013724, 
0.8698943992303472, 
0.866666666666667];

k = 2:15

clf();
plot(k, lambdarfm2, label = "RFM3", linestyle = ":");
scatter(k, lambdarfm2);
plot(k, lambdarfm1, label = "RFM2", linestyle = "--");
scatter(k, lambdarfm1);
plot(k, vm, label = "VM", linestyle = "-", color = "black");
scatter(k, vm);
xlabel("k");
ylabel("second largest eigenvalue")
legend(loc = "lower right");
savefig("figures/eigenvaluecomp.eps")











