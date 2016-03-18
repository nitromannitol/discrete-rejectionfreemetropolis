using PyPlot

#functions used in script.jl

###plotting functions 

#simple command line output print
#note that this prints out backwards 
#from what pyplot prints
function prettyprint(comb, n)
	for i in 1:n 
		for j in 1:n 
			if (i,j) in comb 
				print(" * ")
			else
				print(" - ")
			end
		end
		println("");
	end
end

#plots a given configuration of k points
#on a n by n grid and saves the graph
function plotConfig(config, n, hasTitle = false, titleName = "")
	clf();
	ax = gca();
	ax[:yaxis][:set_ticklabels]([])
	ax[:xaxis][:set_ticklabels]([])
	grid("on");
	for i in 1:n
		for j in 1:n
			if ((i,j) in config)
				scatter(i,j, marker = "o", s = 5000, c= "red");
			else
				scatter(i,j, s = 50)
			end
		end
	end
	if(hasTitle)
		savefig(string("test/", titleName, ".eps"))
		#savefig(string("figures/", titleName, ".eps"))
	end
end
	


#### math functions

#returns periodic manhattan distance of the two points
function manhattanDistance(p1, p2, n)
	i1 = min(p1[1], p2[1]);
	i2 = max(p1[1], p2[1]);

	dist1 = min(abs(i1 - i2), abs((i1+n - i2)));

	j1 = min(p1[2], p2[2]);
	j2 = max(p1[2], p2[2]);

	dist2 = min(abs(j1 - j2), abs((j1+n - j2)));

	return dist1 + dist2; 
end



#check a possible configuration to see if 
#it satisfies the distance greater than 1 apart
#criteria
function checkConfiguration(config, dist, n)
	for point1 in config 
		for point2 in config
			if point1 != point2 && manhattanDistance(point1, point2,n) <= dist
				println(point1, point2, manhattanDistance(point1,point2,n))
				return false;
			end
		end
	end
	return true; 
end

#checks if point collides with any other points
#in the configuration 
#if it does collide, return the point it collides with 
#otherwise return the original point 
function checkForCollisions(point1, config,n, dist)
	for point2 in config
		if(manhattanDistance(point1, point2, n) == dist)
			return point2 
		end
	end
	return point1; 
end


##generate all possible configurations of k particles in an n by n grid 
##constrained to be at distance strictly greater than dist apart 
##do this by looking at all possible configurations and then just checking
##to see if the constraint distance greater than dist apart is satisfies 
function generateConfiguration(n, k, printBool = false, titleName = "", dist = 1)
	#first generate all possible points 
	#on the grid
	tuples = (Int64,Int64)[]; 
	for i in 1:n
		for j in 1:n 
			push!(tuples, (i,j))
		end
	end

	#use built in package to generate all k possible 
	#combinations from tuples 
	possibleCombs = combinations(tuples, k);


	#now check each possible combination 
	#to see if the distance constraint is satisfied
	counter = 0; 
	allowedConfigs = Array{(Int64,Int64),1}[];  
	for (index, comb) in enumerate(possibleCombs)
		if(checkConfiguration(comb, dist, n))
			if(printBool)
				#println("---------------")
				#prettyprint(comb, n);
				plotConfig(comb,n, printBool, string(counter, titleName))
			end
			counter = counter+1;
			comb = sort(comb); #just so it's in a fixed order
			push!(allowedConfigs, comb);
		end
	end
	if(printBool) 
		println("We have generated all ", counter, " possible configurations.");
	end
	return allowedConfigs;
end

#helper function for runCollisionChain
#helps moving up and to the right
#on a periodic graph
function movePoint(point, direction, n)
	if(direction == "right")
		p = point[1];
		if(p == n)
			p = 1; 
		else
			p = p+1; 
		end
		return (p, point[2]);
	elseif(direction == "up")
		p = point[2];
		if(p == n)
			p = 1;
		else
			p = p+1;
		end
		return (point[1], p);
	elseif(direction == "left")
		p = point[1];
		if(p == 1)
			p = n; 
		else
			p = p-1; 
		end
		return (p, point[2]);
	elseif(direction == "down")
		p = point[2];
		if(p == 1)
			p = n;
		else
			p = p-1;
		end
		return (point[1], p);
	end
end


#runs the event chain algorithm starting at the point 
#in the given direction with lambda displacement parameter  
#on the given configuration
#direction can be right or up
#dist - how far away apart the particles are constrained to be before
#a colliison occurs 
#returns the configuration after the chain 
function runCollisionChain(n, point, lambda, direction, config; printName = "", printBool = false, dist = 1, lambdaFixed = false)
	k = length(config);
	if(printBool)
		println("Initial point: ", point, " Direction: ", direction, " Lambda: ", lambda);
	end
	initLambda = lambda; 
	while(true)
			index = findfirst(config, point);
			if(printBool)
				pause(0.25);
				titleName = string(abs(initLambda - lambda),"n=",n,"k=", k, printName);
				plotConfig(config,n, true, titleName);
				#prettyprint(config,n);
				#println("--------");
			end
			newPoint = movePoint(point, direction, n);
			newConfig = [config[1:(index-1)],newPoint, config[(index+1):end]]
			point = checkForCollisions(newPoint, [config[1:(index-1)], config[(index+1):end]],n, dist);
			lambda = lambda-1;
			config = newConfig;
			if(lambda == 0)
				if( (lambdaFixed == false) && ((newPoint != point) || (countnz(newConfig.==(newPoint)) >1)))
																		#check if point is on top of other point
					lambda = lambda+1;
				else
					break;
				end
				#break;
			end
			if(printBool)
				println("Point collided with: ", point);
			end

	end
	if(printBool)
		titleName = string(abs(initLambda - lambda),"n=",n,"k=", k, printName);
		plotConfig(config,n, true, titleName);
		#prettyprint(config,n);
		pause(0.25);
		#println("--------");
	end
	return sort(config); 
end



function generateTransitionMatrix(n,k, lambda, dist = 0, lambdaFixed = false)
	posConfigs = generateConfiguration(n,k, false, "", dist);
	l = length(posConfigs);
	countMatrix = zeros(l,l);

	numInvalid = 0; #keep track of number of invalid configurations
	for (originalIndex, config) in enumerate(posConfigs)
		for point in config
			#run chain up or right
			index1 = findfirst(posConfigs, runCollisionChain(n,point, lambda, "up", config, dist = dist, lambdaFixed = lambdaFixed));
			index2 = findfirst(posConfigs, runCollisionChain(n,point, lambda, "right", config, dist = dist, lambdaFixed = lambdaFixed));


			#there should not be invalid configurations with the updated algorithm 
			numInvalid = numInvalid + (index1==0) + (index2==0);  #increment number of invalid increments
			if( index1 == 0 )
				println("Invalid configuration!");
				#runCollisionChain(n,point, lambda, "up", config, string("invalid", point, "up", lambda), true, dist)
				#rejects = index is same position 
				#index1 = originalIndex
			end
			if( index2 == 0 )
				println("Invalid configuration!");
				#runCollisionChain(n,point, lambda, "right", config, string("invalid", point, "right", lambda), true, dist)

				#rejects = index is same position 
				#index2 = originalIndex
			end


			countMatrix[originalIndex, index1] = countMatrix[originalIndex, index1] +1;
			countMatrix[originalIndex, index2] = countMatrix[originalIndex, index2] +1;
		end

	end
	normalizationConstant = 2*k; 
	stochasticMatrix = countMatrix./normalizationConstant;
	percent = numInvalid/(normalizationConstant*l);

	return stochasticMatrix, numInvalid, percent
end



#simulate the rejection free algorithm 
#starting with initial configuration, initialConfig 
#on an n by n grid with k particles constrained to be at distance dist
#and displacement parameter lambda
#the chain is run for niters and then the expected value of the function
#fn is returned 
function simulateChain(initialConfig, n,k,lambda,niters, fn, dist)

	config = initialConfig;
	chain = fn(config);

	for i in 1:niters
		#pick a point randomly 
		point = initialConfig[rand(1:k)];

		#flip a coin and pick a direction randomly
		direction = "up";
		if(rand(1) > 0.5) 
			direction = "right"
		end

		#run collision chain 
 		config = runCollisionChain(n, point, lambda, direction, config,dist = 0, lambdaFixed = false);
 		chain = chain + fn(config);
	end

	return chain./niters;

end


function runMetropolis(n, point, config, direction, posConfigs, originalIndex)
	index = findfirst(config, point);
	newPoint = movePoint(point, direction, n);
	newConfig = [config[1:(index-1)],newPoint, config[(index+1):end]]
	point = checkForCollisions(newPoint, [config[1:(index-1)], config[(index+1):end]],n, 0);
	if( (newPoint != point) || (countnz(newConfig.==(newPoint)) >1))
		index = originalIndex
	else
		newConfig = sort(newConfig);
		index = findfirst(posConfigs, newConfig);
		if(index == 0)
			println(config);
			println(newConfig);
		end
	end
	return index;
end


function generateTransitionMatrixMetropolis(n,k)
	posConfigs = generateConfiguration(n,k, false, "", 0);
	l = length(posConfigs);
	countMatrix = zeros(l,l);

	for (originalIndex, config) in enumerate(posConfigs)
		for point in config

			#run meropolis up, down, left, or right
			ind1 = runMetropolis(n,point,config,"up", posConfigs, originalIndex);
			ind2 = runMetropolis(n,point,config,"down", posConfigs, originalIndex);
			ind3 = runMetropolis(n,point,config,"left", posConfigs, originalIndex);
			ind4 = runMetropolis(n,point,config,"right", posConfigs, originalIndex);

			countMatrix[originalIndex, ind1] = countMatrix[originalIndex, ind1] +1;
			countMatrix[originalIndex, ind2] = countMatrix[originalIndex, ind2] +1;
			countMatrix[originalIndex, ind3] = countMatrix[originalIndex, ind3] +1;
			countMatrix[originalIndex, ind4] = countMatrix[originalIndex, ind4] +1;
		end

	end
	normalizationConstant = 4*k; 
	stochasticMatrix = countMatrix./normalizationConstant;

	return stochasticMatrix

end




