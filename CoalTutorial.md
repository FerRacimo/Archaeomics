Exercises in coalescent theory
===============

#Exercise	1A:	Simulating	a	coalescence	tree	assuming	a	constant	population	size

The	purpose	of	this first	exercise	is	to	make	sure	it	is	clear	how	a	coalescence tree	is	simulated.	To	
do	this	let	us	try to	simulate	a	coalescence tree	for	five gene	copies by	hand (+ a	little	help	from	R):
1. Start	by	drawing	a	node	for	each	of	the	five	gene	copies on	an	invisible line	(with	space	for	
drawing	a	tree	above them).	Name	these	nodes 1,2,3,4,5
2. Also,	make	a	list	of	the	node	names.	You	can	either	do	this	by	hand	or	you	can	do	it	in	R	by	
simply	writing	
nodes = c(1,2,3,4,5) # make the list and call it nodes
nodes # print the list
3. Sample	which	two	nodes	will	coalesce	first	(going	back	in	time)	by	randomly	picking	two	of	the	
nodes.	You	can	either	do	this	by	hand	or	you	can	do	it	in	R	by	typing
nodecount = length(nodes) # save the number of nodes in the variable nodecount
tocoalesce = sample(1:nodecount, size=2) # sample 2 different nodes in node list
nodes[tocoalesce[1]] # print the first node sampled
nodes[tocoalesce[2]] # print the second node sampled
If	you	used	R	then	make	sure	you	understand	what	the	R	code	does	before	moving	on.
4. Sample	the	time	it	takes	before	these	two	nodes	coalesce	(measured	from	previous	
coalescence	event in	units	of	2N)	by	sampling	from	an	exponential	distribution	with	rate	equal	
to	nodecount*(nodecount-1)/2	where	nodecount	is	the	number	of	nodes	in	your node	list.	Do	
this	in	R	by	typing:
nodecount = length(nodes) # save the number of nodes in the variable nodecount
coalescencerate = nodecount*(nodecount-1)/2 # calculate the coalescent rate
coalescencetime = rexp(1, rate=coalescencerate) # sample from exponential w. that rate
coalescencetime
Make sure	you	understand	what	the	R	code	does	before	moving	on.
5. Now	draw	a	node	that	is	the	sampled amount	of	time	further	up	in	the	tree	than	the currently	
highest node	(so	if	the	currently	highest	node	is	drawn	at	height	T	then	draw	the	new	one	at	
height	T plus the	sampled	coalescence	time)	and	draw	a	branch	from	each	of	the	nodes	you	
sampled	in	step	3	to	this	new	node	indicating	that	these	two	nodes	coalesce at	this	time.	
6. Next,	make	an	updated	list	of	the	nodes	that	are	left	by	removing	the two	nodes	that	
coalesced	and	instead	adding	the	newly	drawn	node	that represents	their	common	ancestor.	
You	can	call	the	new	node	the	next	number	not	used	as	a	name	yet	(e.g. if	this	is	the	first	
coalescence event you	can	call	it	6, if	it	is	the	second	coalescence	event you	can	call	it	7	etc.).	
You	can	either	do	this	by	hand	or	in	R.	If	you	want	to	do	it	R	you	can	do	it	as	follows:
nodes <- nodes[-tocoalesce] # remove the two nodes that coalesced
nodes <- c(nodes,2*5-length(nodes)-1) # add the new node
nodes # print the new list
If	you	used	R	then	make	sure	you	understand	what	the	R	code	does	before	moving	on.
