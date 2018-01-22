Tutorial based on notes by Mikkel Schierup and Ida Moltke

# Effect	of	population	size	on	coalescence	trees

The	time	scale	in	the	trees	we	have	been	looking	at	so	far	is	2N	(assuming	we	are	looking	at	a	
diploid	organism),	where	N	is	the	population	size.	With	this	in	mind	answer	the	following	
questions:

1) How	would	you	expect	a	tree	for	five	gene	copies	from	a	population	of	size	1000	to	differ	from
a	tree	for	five	gene	copies	from	population	of	size	2000? (if	we	drew	them	in	real	time	scale)

2) The	tree	below	is	simulated	assuming	a	constant	population	size:

![alt text](https://github.com/FerRacimo/DemographicCourseAdelaide2018/blob/master/Tree1.png)

How	would	the	shape	of	this tree	be	different if	the	population	size	had	not	been	constant, but	
instead	to	population	size	had	been	10	times	bigger	from the	time	point	indicated	with	^ to the	
present?	(i.e. the	population	had	increased	dramatically	in	size	recently)

3) And	what	if	the	population	size	change	had	not	been	an	increase	but	instead	a	decrease?

4) And	by	the	same	logic, which	of	the	following	trees	do you	think	are	simulated	under	
decreasing	population	size,	constant	population	size	and	increasing	population	size?

![alt text](https://github.com/FerRacimo/DemographicCourseAdelaide2018/blob/master/Tree2.png)

# The	connection	between	a	mutation’s	location	on	a	tree	and	the	SFS

1) Which	categories	would	a	Site	Frequency	Spectrum	(SFS)	for	five	gene	copies have?
(e.g. singletons,	doubletons,…?)

2) Which	of	these	categories	would	the	mutation	in	this	tree	(marked	by	a	red	star)	contribute	to:

![alt text](https://github.com/FerRacimo/DemographicCourseAdelaide2018/blob/master/Tree3.png)

3) Which	category	of	the	SFS	would	this	mutation	contribute	to	in	an	SFS:	

![alt text](https://github.com/FerRacimo/DemographicCourseAdelaide2018/blob/master/Tree4.png)

4) Draw	a	SFS	based	on	this	tree:	

![alt text](https://github.com/FerRacimo/DemographicCourseAdelaide2018/blob/master/Tree5.png)

# Effect	of	population	size	changes	on	the	SFS
Usually	the	number	of	mutations	that	happen on	a	branch	is	proportional	to	the	length	of	the	
branch. With	this	and	your	answers	in exercise	1C	and	1D	in	mind answer	the	following	questions:

1) How	do	you	think	population	growth affects the	frequency	spectrum?	

2) How	do	you	think	population	size decrease affects the	frequency	spectrum?


# Recombination and the Coalescent

Let's go back to using the Hudson animator. The recombination rate is determined by rho=2Nc. In the animation, recombination events are marked as blue nodes (in contrast to the green nodes of coalescent events).
Set n=5, and rho=1. Press recalc and study the animation. Does any recombination events occur? If not, try again. How many recombination events do you expect to see? Try a total of 25 times and write down the number of recombination events in each case. What does the distribution look like?

Look at a couple of simulations in more detail. Study where in the sequence, recombination events occur. Can you find examples that different part of the sequence have different most recent common ancestors (marked in green) at different time points?

Try to press the trees window pane. Here it is possible to study each of the different trees over the sequence. How many different trees are there and how does this relate to the number of recombination events? Try to find examples of each of the following recombination types:

1) Recombination changing the topology of the tree

2) Recombination changing the branch length but not the topology of the tree

3) Recombination that does not change the tree

Try setting rho=2 and n=10
How many recombinations occur now?
At which time do different parts of the sequence find a most recent common ancestor?
What is the time until the first part of the sequence finds a most recent common ancestor (calculate the average over 5 replicates). At which time have all the different parts of the sequence found a most recent common ancestor (average over 5 replicates)? Compare this time with the time to most recent common ancestor without recombination


Try setting rho=5 and n=20, how many recombination events occur now? How many different most recent common ancestors are there over the sequence now?
