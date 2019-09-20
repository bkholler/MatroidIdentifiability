(* ::Package:: *)

(*This should be the directory that all of the files are in*)
SetDirectory["/Users/bkholler/Downloads/TreeIDFor_K3P_CFN/FilesForPaper"];


jacobianMatrix[map_, vars_] := Module[{curRow,outMat},
	outMat = {};
	
	Do[
		curRow = {}; 
		
		Do[
			AppendTo[curRow, D[map[[i]],vars[[j]]]]
			
		,{j,1,Length[vars]}];
		
		AppendTo[outMat,curRow];
		
	,{i,1,Length[map]}];
	
outMat
];


fallFact[x_,k_]:= Product[x-i,{i,0,k-1}];


(* Permutes the leaf labels of the tree defined by splits with permutation perm *)
permuteTree[splits_,perm_] := Sort[Map[Sort,Map[Sort,splits/. perm,2]]];


(* pairs, a pair of pairs of trees given by their splits,
 * perm, a permutation of the leaves,
 * treeSet, the set of trees with the number of leaves that the trees in pairs have,
 * Returns a the position of each tree in treeSet after permuting the leaves of each tree by perm as a pair of pairs
*)
permutePairs[pairs_, perm_,treeSet_] := Module[{t1,t2,t3,t4},
t1 = Position[treeSet,   permuteTree[treeSet[[pairs[[1,1]]]],perm]][[1,1]];
t2 = Position[treeSet,   permuteTree[treeSet[[pairs[[1,2]]]],perm]][[1,1]];
t3 = Position[treeSet,   permuteTree[treeSet[[pairs[[2,1]]]],perm]][[1,1]];
t4 = Position[treeSet,   permuteTree[treeSet[[pairs[[2,2]]]],perm]][[1,1]];

Sort[{Sort[{t1,t2}],Sort[{t3,t4}]}]
];


(* Converts a permutation to a list of rules that can be used with /. *)
permToRules[perm_] := Module[{outRules},
outRules= {};
Do[
AppendTo[outRules,i -> perm[[i]]]; 
,{i,1,Length[perm]}];

outRules
];


(* inSet, the original set of all pairs of pairs, given as pairs of integers encoding the position of the trees in treeSet,
 * m, the number of leaves that the trees in the pairs have,
 * treeSet, the set of trees with m leaves given by their splits,
 * Returns a list of representatives for each orbit given. Each pair of pairs is an integer pair encoding the position of the trees in treeSet
*)
findOrbitReps[inSet_,m_,treeSet_] := Module[{remainingPairs,orbReps,curOrb,i},
remainingPairs= inSet;

orbReps = {};

While[Length[remainingPairs] > 0,

curOrb = {remainingPairs[[1]]};

Do[
AppendTo[curOrb, permutePairs[curOrb[[1]],permToRules[Permutations[Range[m]][[i]]], treeSet]];
,{i,1,m!}];

curOrb = DeleteDuplicates[curOrb];

AppendTo[orbReps,curOrb[[1]]];

remainingPairs = Complement[remainingPairs,curOrb];
];

orbReps
];


(* tree, a leaf labelled binary tree given by its splits
* Returns the A set in each split {A,B}
*)
splitSetToA[tree_] := Module[{outSplits},
outSplits = {};

Do[
AppendTo[outSplits,tree[[i,1]]] 
,{i,1,Length[tree]}];

outSplits
];


(* Converts a list of integers {n_1,n_2,...,n_k} to the string "n_1n_2...n_k" *)
intVectToString[intVect_] := Module[{outString},
outString = "";

Do[
outString = outString<>ToString[intVect[[i]]];
,{i,1,Length[intVect]}];

outString
];



(* vars, a list of variables,
 * vals, a list of values those variables can take,
 * Returns a list of rules of the form vars[[i]] \[Rule] vals[[i]] which can be used to substitue values with /.
 *)
varSubs[vars_, vals_] := Module[{outRules},
outRules = {};

Do[
AppendTo[outRules,vars[[i]]-> vals[[i]]];
,{i,1,Length[vars]}];

outRules
];


randomRealParameters[symbolicFunction_] := symbolicFunction /. varSubs[Variables[symbolicFunction], RandomReal[{-1, 1}, Length[Variables[symbolicFunction]]]];


randomIntParameters[poly_,p_] := poly /. varSubs[Variables[poly],  RandomInteger[{0,p-1}, Length[Variables[poly]]]];





(* m, the number of leaves of the trees that intPairs represent,
 * intPairs, a pair of integers representing trees by their position in the set of unrooted m leaf trees,
 * Returns the actual set of trees as a set of their A-splits
*)
intPairsToTrees[m_,intPairs_] := Map[splitSetToA, {{unrootedTrees[[m,intPairs[[1, 1]]]],unrootedTrees[[m,intPairs[[1, 2]]]]},
{unrootedTrees[[m,intPairs[[2, 1]]]], unrootedTrees[[m,intPairs[[2, 2]]]]}},{2}];


(* Unrooted [n]-trees given by their splits A|B *)
unrooted4 = Get["unrooted4"];
unrooted5 = Get["unrooted5"];
unrooted6 = Get["unrooted6"];
unrootedTrees = {0,0,0,unrooted4,unrooted5,unrooted6};


(*integer pairs of pairs. The integers correspond to the position of each tree in the set of unrooted trees.*)
fourLeafOrbs = Get["fourLeafOrbs"];
fiveLeafOrbs = Get["fiveLeafOrbs"];
sixLeafOrbs = Get["sixLeafOrbs"];
certsCFN = Get["certsCFN"];


remove0Coords[list_] := list[[Complement[Range[Length[list]],Select[Range[Length[list]],list[[#]]=== 0&]]]];


(* puts a total ordering on the Klein group Z_2 x Z_2 *)
kleinPos[el_] := Position[{{0,0},{1,0},{0,1},{1,1}},el][[1,1]];


(* m, the number of leaves of tree
 * tree, a tree given by the A-sets of its splits A|B
 * paramSym, a string to be used as a formal variable in the parameterization
 * model, a string specifying the model (CFN, JC, K2P, K3P)
 * IDParamTo1, a boolean set to True by default. This automatically sets the parameters associated to the identity element to 1 if true
 * Outputs a matrix of parameters to be used in the parameterization map of a group-based phylogenetic tree model in the Fourier coordinates
*)
phyloTreeParameterMatrix[m_, tree_, paramSym_, model_, IDParamTo1_: True] := Module[{parameterMatrix},
parameterMatrix = {};

If[model == "CFN",

	Do[
		AppendTo[parameterMatrix, {paramSym<>"A"<>ToString[i], paramSym<>"C"<>ToString[i]}];
	
	,{i, 1, m}];

	Do[
		AppendTo[parameterMatrix, {paramSym<>"A"<>intVectToString[tree[[i]]],paramSym<>"C"<>intVectToString[tree[[i]]]}];
	
	,{i, 1, Length[tree]}];	

(* else if *),

If[model == "K3P",

	Do[
		AppendTo[parameterMatrix,{paramSym<>"A"<>ToString[i],paramSym<>"C"<>ToString[i],paramSym<>"G"<>ToString[i],paramSym<>"T"<>ToString[i]}];
	
	,{i,1,m}];

	Do[
		AppendTo[parameterMatrix,{paramSym<>"A"<>intVectToString[tree[[i]]],paramSym<>"C"<>intVectToString[tree[[i]]],paramSym<>"G"<>intVectToString[tree[[i]]],paramSym<>"T"<>intVectToString[tree[[i]]]}];
	
	,{i,1,Length[tree]}];	
	
(* else if *),

If[model == "K2P",

	Do[
		AppendTo[parameterMatrix,{paramSym<>"A"<>ToString[i],paramSym<>"C"<>ToString[i],paramSym<>"G"<>ToString[i],paramSym<>"C"<>ToString[i]}];
	
	,{i,1,m}];

	Do[
		AppendTo[parameterMatrix,{paramSym<>"A"<>intVectToString[tree[[i]]],paramSym<>"C"<>intVectToString[tree[[i]]],paramSym<>"G"<>intVectToString[tree[[i]]],paramSym<>"C"<>intVectToString[tree[[i]]]}];
	
	,{i,1,Length[tree]}];	

(* else if *),

If[model == "JC",

	Do[
		AppendTo[parameterMatrix,{paramSym<>"A"<>ToString[i],paramSym<>"C"<>ToString[i],paramSym<>"C"<>ToString[i],paramSym<>"C"<>ToString[i]}];
	
	,{i,1,m}];

	Do[
		AppendTo[parameterMatrix,{paramSym<>"A"<>intVectToString[tree[[i]]],paramSym<>"C"<>intVectToString[tree[[i]]],paramSym<>"C"<>intVectToString[tree[[i]]],paramSym<>"C"<>intVectToString[tree[[i]]]}];
	
	,{i,1,Length[tree]}];
];
];
];
];

parameterMatrix = Map[ToExpression, parameterMatrix,2];

If[IDParamTo1,

	parameterMatrix = parameterMatrix /. varSubs[parameterMatrix[[All,1]], ConstantArray[1, Length[parameterMatrix[[All,1]]]]]
];

parameterMatrix
];


(* m, the number of leaves of tree
 * tree, a tree given by the A-sets of its splits A|B
 * paramSym, a string to be used as a formal variable in the parameterization
 * model, a string specifying the model (CFN, JC, K2P, K3P)
 * IDParamTo1, a boolean set to True by default. This automatically sets the parameters associated to the identity element to 1 if true
 * Outputs the Jacobian of the parameterization map of a group-based phylogenetic tree model in the Fourier coordinates
*)
phyloTreeParameterization[m_, tree_, paramSym_, model_, IDParamTo1_: True] := Module[{outMap, G, Gtuples, tupleProb, tuple, parameterMatrix},
outMap = {};

If[model == "CFN",

	G = {0,1};
	Gtuples = Tuples[G,m];
	parameterMatrix = phyloTreeParameterMatrix[m, tree, paramSym, model, IDParamTo1];

	Do[
	
		tupleProb = 0;
		tuple = Gtuples[[i]];

		If[Mod[Total[tuple],2] == 0,

			tupleProb = Product[parameterMatrix[[j]][[tuple[[j]] + 1]], {j, 1, m}];

			Do[
				tupleProb *= parameterMatrix[[j+m]][[Mod[Total[tuple[[tree[[j]]]]], 2] + 1]];
			,{j, 1, Length[tree]}];
		];

		AppendTo[outMap, tupleProb];

	,{i, 1, 2^m}];
	
	(* else *),

	G = {{0,0},{1,0},{0,1},{1,1}}; (*A,C,G,T*)
	Gtuples = Tuples[G,m];
	parameterMatrix = phyloTreeParameterMatrix[m, tree, paramSym, model, IDParamTo1];

	Do[
	
		tuple = Gtuples[[i]];
		tupleProb = 0;
	
		If[Mod[Total[tuple], 2] == {0,0},

			tupleProb = Product[parameterMatrix[[j]][[kleinPos[tuple[[j]]]]], {j, 1, m}];

			Do[
				tupleProb *= parameterMatrix[[j+m]][[kleinPos[Mod[Total[tuple[[tree[[j]]]]], 2]]]]

			,{j,1,Length[tree]}];
		];

		AppendTo[outMap, tupleProb];

	,{i, 1, 4^m}];
];

outMap
];


(* m, the number of leaves of tree
 * tree, a tree given by the A-sets of its splits A|B
 * paramSym, a string to be used as a formal variable in the parameterization
 * homogenizeParam, a formal parameter that can be used to rehomogenize the parameterization if the parameters associated to the identity element were set to 0
 * model, a string specifying the model (CFN, JC, K2P, K3P)
 * IDParamTo1, a boolean set to True by default. This automatically sets the parameters associated to the identity element to 1 if true
 * Outputs the Jacobian of the parameterization map of a group-based phylogenetic tree model in the Fourier coordinates
*)
phyloTreeJac[m_, tree_ , paramSym_, homogenizeParam_, model_, IDParamTo1_: True] := Module[{map, vars},
map = homogenizeParam*phyloTreeParameterization[m, tree, paramSym, model, IDParamTo1];
vars = Flatten[phyloTreeParameterMatrix[m, tree, paramSym, model, IDParamTo1]];
vars = Join[{homogenizeParam}, vars];
vars = DeleteDuplicates[DeleteCases[vars , 1]];

jacobianMatrix[map, vars]
];


(* m, the number of leaves of tree
 * pair, a pair of trees each given by the A-sets of its splits A|B
 * paramSym1, a string to be used as a formal variable in the parameterization for the first tree
 * paramSym2, a string to be used as a formal variable in the parameterization for the first tree
 * mixParam, a formal parameter that can be used to form the mixture of the two parameterizations
 * model, a string specifying the model (CFN, JC, K2P, K3P)
 * IDParamTo1, a boolean set to True by default. This automatically sets the parameters associated to the identity element to 1 if true
 * Outputs the parameterization map of a 2-tree group-based mixture model in the Fourier coordinates
*)
phyloTreeMixParameterization[m_, pair_, paramSym1_,paramSym2_, mixParam_, model_, IDParamTo1_: True] := mixParam*phyloTreeParameterization[m, pair[[1]], paramSym1, model, IDParamTo1] + (1-mixParam)*phyloTreeParameterization[m, pair[[2]], paramSym2, model, IDParamTo1]


(* m, the number of leaves of tree
 * pair, a pair of trees each given by the A-sets of its splits A|B
 * paramSym1, a string to be used as a formal variable in the parameterization for the first tree
 * paramSym2, a string to be used as a formal variable in the parameterization for the first tree
 * mixParam, a formal parameter that can be used to form the mixture of the two parameterizations
 * model, a string specifying the model (CFN, JC, K2P, K3P)
 * IDParamTo1, a boolean set to True by default. This automatically sets the parameters associated to the identity element to 1 if true
 * Outputs the Jacobian of the parameterization map of a 2-tree group-based mixture model in the Fourier coordinates 
*)
phyloTreeMixJac[m_, pair_, paramSym1_, paramSym2_, mixParam_, model_, IDParamTo1_: True] := Module[{map, vars},
map = phyloTreeMixParameterization[m, pair, paramSym1, paramSym2, mixParam, model, IDParamTo1];
vars = {mixParam};

If[model == "CFN",
	vars = Join[vars, Flatten[phyloTreeParameterMatrix[m, pair[[1]], paramSym1, model, IDParamTo1][[All,{1,2}]]], Flatten[phyloTreeParameterMatrix[m, pair[[2]], paramSym2, model, IDParamTo1][[All,{1,2}]]]];
];

If[model == "K3P" || model == "K2P" || model == "JC",
	vars = Join[vars, Flatten[phyloTreeParameterMatrix[m, pair[[1]], paramSym1, model, IDParamTo1]], Flatten[phyloTreeParameterMatrix[m, pair[[2]], paramSym2, model, IDParamTo1]]];
];

vars = DeleteDuplicates[DeleteCases[vars, 1]];

jacobianMatrix[map, vars]
];


(* m, the number of leaves of the network
 * edgeSet, a list with elements of the form {edgeLabel, A-set of split} corresponding to the tree that one gets when a reticulation edge is deleted from the network
 * paramSym, a string to be used as a formal variable in the parameterization
 * model, a string specifying the model (CFN, JC, K2P, K3P)
 * IDParamTo1, a boolean set to True by default. This automatically sets the parameters associated to the identity element to 1 if true
 * Outputs a matrix of parameters to be used in the parameterization map of the k-cycle group-based phylogenetic network model in the Fourier coordinates
*)
phyloNetParameterMatrix[m_, edgeSet_, paramSym_, model_, IDParamTo1_: True] := Module[{parameterMatrix},

parameterMatrix = {};

If[model == "CFN",

	Do[
		AppendTo[parameterMatrix, {paramSym<>"A"<>ToString[edgeSet[[i,1]]], paramSym<>"C"<>ToString[edgeSet[[i,1]]]}];
	
	,{i, 1, Length[edgeSet]}];	

(* else if *),

If[model == "K3P",

	Do[
		AppendTo[parameterMatrix,{paramSym<>"A"<>ToString[edgeSet[[i,1]]], paramSym<>"C"<>ToString[edgeSet[[i,1]]], paramSym<>"G"<>ToString[edgeSet[[i,1]]], paramSym<>"T"<>ToString[edgeSet[[i,1]]]}];
	
	,{i, 1, Length[edgeSet]}];	
	
(* else if *),

If[model == "K2P",

	Do[
		AppendTo[parameterMatrix,{paramSym<>"A"<>ToString[edgeSet[[i,1]]], paramSym<>"C"<>ToString[edgeSet[[i,1]]], paramSym<>"G"<>ToString[edgeSet[[i,1]]], paramSym<>"C"<>ToString[edgeSet[[i,1]]]}];
	
	,{i, 1, Length[edgeSet]}];	

(* else if *),

If[model == "JC",

	Do[
		AppendTo[parameterMatrix,{paramSym<>"A"<>ToString[edgeSet[[i,1]]], paramSym<>"C"<>ToString[edgeSet[[i,1]]], paramSym<>"C"<>ToString[edgeSet[[i,1]]], paramSym<>"C"<>ToString[edgeSet[[i,1]]]}];
	
	,{i,1,Length[edgeSet]}];
];
];
];
];

parameterMatrix = Map[ToExpression, parameterMatrix, 2];

If[IDParamTo1,

	parameterMatrix = parameterMatrix /. varSubs[parameterMatrix[[All,1]], ConstantArray[1, Length[parameterMatrix[[All,1]]]]]
];

parameterMatrix
];


(* m, the number of leaves of the network
 * network, a k-cycle phylogenetic network in the form {tree1, tree2} where tree1 and tree2 are the trees obtained by deleting a reticulation edge from the network
 * paramSym, a string to be used as a formal variable in the parameterization
 * model, a string specifying the model (CFN, JC, K2P, K3P)
 * IDParamTo1, a boolean set to True by default. This automatically sets the parameters associated to the identity element to 1 if true
 * Outputs the parameterization map of the k-cycle group-based phylogenetic network model in the Fourier coordinates
*)
phyloNetParameterization[m_, network_, paramSym_, mixParam_, model_, IDParamTo1_: True] := Module[{outMap, G, Gtuples, parameterMatrix1, parameterMatrix2, tuple, tupleProb},

outMap = {};

If[model == "CFN",

	G = {0,1};
	Gtuples = Tuples[G,m];
	parameterMatrix1 = phyloNetParameterMatrix[m, network[[1]], paramSym, model, IDParamTo1];
	parameterMatrix2 = phyloNetParameterMatrix[m, network[[2]], paramSym, model, IDParamTo1];
	
	Do[
		tuple = Gtuples[[i]];
		tupleProb = 0;
	
		If[Mod[Total[tuple],2] == 0,
			
			tupleProb = mixParam*(Product[parameterMatrix1[[j]][[Mod[Total[tuple[[network[[1, j, 2]]]]],2] + 1]], {j, 1, Length[network[[1]]]}] + Product[parameterMatrix2[[j]][[Mod[Total[tuple[[network[[2, j, 2]]]]], 2] + 1]], {j, 1, Length[network[[2]]]}])
		];
		
		AppendTo[outMap, tupleProb];
		
	,{i, 1, 2^m}];
];

If[model !=  "CFN",

	G = {{0,0},{1,0},{0,1},{1,1}}; (*A,C,G,T*)
	Gtuples = Tuples[G,m];
	parameterMatrix1 = phyloNetParameterMatrix[m, network[[1]], paramSym, model, IDParamTo1];
	parameterMatrix2 = phyloNetParameterMatrix[m, network[[2]], paramSym, model, IDParamTo1];
	
	Do[
		tuple = Gtuples[[i]];
		tupleProb = 0;
	
		If[Mod[Total[tuple],2] == {0,0},
			
			tupleProb = mixParam*(Product[parameterMatrix1[[j]][[kleinPos[Mod[Total[tuple[[network[[1,j,2]]]]],2]]]],{j,1,Length[network[[1]]]}] + Product[parameterMatrix2[[j]][[kleinPos[Mod[Total[tuple[[network[[2,j,2]]]]],2]]]],{j,1,Length[network[[2]]]}])
		];
		
		AppendTo[outMap, tupleProb];
		
	,{i, 1, 4^m}];
];

remove0Coords[outMap]
];


(* m, the number of leaves of the network
 * network, a k-cycle phylogenetic network in the form {tree1, tree2} where tree1 and tree2 are the trees obtained by deleting a reticulation edge from the network
 * paramSym, a string to be used as a formal variable in the parameterization
 * model, a string specifying the model (CFN, JC, K2P, K3P)
 * IDParamTo1, a boolean set to True by default. This automatically sets the parameters associated to the identity element to 1 if true
 * Ouputs the Jacobian of the parameterization of a k-cycle group-based phylogenetic network model in the Fourier coordinates
*)
phyloNetJac[m_, network_, paramSym_, mixParam_, model_, IDParamTo1_: True] := Module[{map, vars},

map = phyloNetParameterization[m, network, paramSym, mixParam, model, IDParamTo1];
vars = {mixParam};

If[model == "CFN",
	vars = Join[vars, Flatten[phyloNetParameterMatrix[m, network[[1]], paramSym, model, IDParamTo1]], Flatten[phyloNetParameterMatrix[m, network[[2]], paramSym, model, IDParamTo1]]];
];

If[model == "K3P" || model == "K2P" || model == "JC",
	vars = Join[vars, Flatten[phyloNetParameterMatrix[m, network[[1]], paramSym, model, IDParamTo1]], Flatten[phyloNetParameterMatrix[m, network[[2]], paramSym, model, IDParamTo1]]];
];

vars = DeleteDuplicates[DeleteCases[vars, 1]];

jacobianMatrix[map, vars]
];


(* jac1, the jacobian of a polynomial map phi_1 parameterizing a statistical model M_1
 * jac2, the jacobian of a polynomial map phi_2 parameterizing a statistical model M_2
 * dim1, the dimension of M_1
 * dim2, the dimension of M_2
 * startSize, the smallest set size to use in Algorithm 3.3.
 * trials, for j such that startSize \[LessEqual] j \[LessEqual] min{dim1, dim2}, this algorithm will randomly sample trials number of sets of size j and see if they are certificates
 * Outputs a certificate that dim(M_1 \cap M_2) < min(dim1, dim2)
*)
matroidSeparate[jac1_, jac2_, dim1_, dim2_, startSize_, trials_] := Module[{numJac1, numJac2, trialSet, foundCert, certSet},

foundCert = False;
numJac1 = randomRealParameters[jac1];
numJac2 = randomRealParameters[jac2];

If[dim1 === dim2,
		
	Do[
		Do[
			trialSet = RandomSample[Range[Length[jac1]], j];
				
			If[MatrixRank[numJac1[[trialSet]]] != MatrixRank[numJac2[[trialSet]]],
				
				If[MatrixRank[jac1[[trialSet]]] != MatrixRank[jac2[[trialSet]]],
					
					foundCert = True;
					certSet = trialSet;
					Break[];
				];
			];
			
			, {i, 1, trials}];
			
		If[foundCert,
			
			Break[]
		];
		
	, {j, startSize, dim1}];
];

If[dim1 > dim2,
		
	Do[
		Do[
			trialSet = RandomSample[Range[Length[jac1]], j];
				
			If[MatrixRank[numJac1[[trialSet]]] < j && MatrixRank[numJac2[[trialSet]]] == j,
				
				If[MatrixRank[jac1[[trialSet]]] < j && MatrixRank[jac2[[trialSet]]] == j,
					
					foundCert = True;
					certSet = trialSet;
					Break[];
				];
			];
			
			, {i, 1, trials}];
			
		If[foundCert,
			
			Break[]
		];
		
	, {j, startSize, dim2}];
];

If[!foundCert,
	
	certSet = {-1}
];

certSet
]



matroidSeparateSZ[jac1_, jac2_, dim1_, dim2_,maxDegree_, startSize_, trials_, epsilon_] := Module[{p, l, numJac1, numJac2,trialSet, foundCert, certSet, alpha},

foundCert = False;

If[dim1 === dim2,

	Do[
	alpha  = j*maxDegree;
 p = NextPrime[alpha,2];
 l = Ceiling[Log[epsilon]/Log[alpha/(p-1)]];

 numJac1 =randomIntParameters[jac1,p];
 numJac2 = randomIntParameters[jac2, p];

	       Do[
			trialSet = RandomSample[Range[Length[jac1]], j];

			If[MatrixRank[numJac1[[trialSet]],Modulus->p] != MatrixRank[numJac2[[trialSet]], Modulus -> p],
		
		Do[
		
		numJac1 = Mod[randomIntParameters[jac1,p],p];
		numJac2 = Mod[randomIntParameters[jac2, p],p];

		If[ MatrixRank[numJac1[[trialSet]],Modulus->p] === MatrixRank[numJac2[[trialSet]], Modulus -> p],

			Break[];
		];

	    ,{k,1,l}];

		certSet = trialSet;
		foundCert = True;
		Break [];
			];
			
			, {i, 1, trials}];
			
		If[foundCert,
			
			Break[]
		];
		
	, {j, startSize, dim1}];
];

If[dim1 > dim2,
		
	Do[
		Do[
			trialSet = RandomSample[Range[Length[jac1]], j];
				
			If[MatrixRank[numJac1[[trialSet]],Modulus->p] < j &&  MatrixRank[numJac2[[trialSet]], Modulus -> p] == j,
		
		Do[
		
		numJac1 = Mod[randomIntParameters[jac1,p],p];
		numJac2 = Mod[randomIntParameters[jac2, p],p];

		If[MatrixRank[numJac1[[trialSet]],Modulus->p] < j &&  MatrixRank[numJac2[[trialSet]], Modulus -> p] == j,

			Break[];
		];

	    ,{k,1,l}];

		certSet = trialSet;
		foundCert = True;
		Break [];
			];
			
			, {i, 1, trials}];
			
		If[foundCert,
			
			Break[]
		];
		
	, {j, startSize, dim1}];
];

If[!foundCert,
	
	certSet = {-1}
];

certSet
];


(* jac1, the jacobian of a polynomial map phi_1 parameterizing a statistical model M_1
 * jac2, the jacobian of a polynomial map phi_2 parameterizing a statistical model M_2
 * dim1, the dimension of M_1
 * dim2, the dimension of M_2
 * certSet, a certificate set confirming that dim(M_1 \cap M_2) < min(dim1, dim2)
 * Confirms that certSet is an independent set for the matroid associated to M_1 and not M_2 or vice versa symbolically. This proves that dim(M_1 \cap M_2) < min(dim1, dim2). 
*)
confirmMatroidCertSym[jac1_, jac2_, dim1_, dim2_, certSet_] := Module[{isCert},

isCert = False;

If[dim1 == dim2,

	If[MatrixRank[jac1[[certSet]]] != MatrixRank[jac2[[certSet]]],
		
		isCert = True
	];
];

If[dim1 > dim2,

	If[MatrixRank[jac1[[certSet]]] < Length[certSet] && MatrixRank[jac2[[certSet]]] == Length[certSet],
		
		isCert = True
	];
];

isCert
];


(* jac1, the jacobian of a polynomial map phi_1 parameterizing a statistical model M_1
 * jac2, the jacobian of a polynomial map phi_2 parameterizing a statistical model M_2
 * dim1, the dimension of M_1
 * dim2, the dimension of M_2
 * certSet, a certificate set confirming that dim(M_1 \cap M_2) < min(dim1, dim2)
 * Confirms that certSet is an independent set for the matroid associated to M_1 and not M_2 or vice versa numerically. This does not prove that dim(M_1 \cap M_2) < min(dim1, dim2). 
*)
confirmMatroidCertNum[jac1_, jac2_, dim1_, dim2_, certSet_] := Module[{isCert, numJac1, numJac2},

isCert = False;
numJac1 = randomRealParameters[jac1];
numJac2 = randomRealParameters[jac2];

If[dim1 == dim2,

	If[MatrixRank[numJac1[[certSet]]] != MatrixRank[numJac2[[certSet]]],
		
		isCert = True
	];
];

If[dim1 > dim2,

	If[MatrixRank[numJac1[[certSet]]] < Length[certSet] && MatrixRank[numJac2[[certSet]]] == Length[certSet],
		
		isCert = True
	];
];

isCert
];


(* m, the number of leaves of network
 * network, a k-cycle network given as a list of directed edges with loops at each leaf
 * retEdges, the reticulation edges of the network
 * Outputs a representation of the network as {tree1, tree2} where tree1 and tree2 are the trees obtained by deleting the edges in retEdges. These trees are given as a list with elements of the form {i, split associated to i} 
*)
netToEdgeSplits[m_, network_, retEdges_] := Module[{trees, edgeSplits, curSplitSet},
trees = {Delete[network, retEdges[[1]]], Delete[network, retEdges[[2]]]};
edgeSplits = {};

Do[
curSplitSet = {};

Do[
AppendTo[curSplitSet,{Position[network,trees[[i,j]]][[1,1]],Sort[Map[Sort,Map[Select[#<= m&],ConnectedComponents[Graph[Delete[trees[[i]],j],DirectedEdges->False]]]]][[1]]}];
,{j,1,2m-1}];

AppendTo[edgeSplits,curSplitSet];
,{i,1,2}];

edgeSplits
];


(* This is a small library of 4 and 5 leaf networks. *)

edgesL4C4 = {1->5, 2->6, 3->7, 4->8, 5->6, 6->7, 7->8, 8->5};
leafSet4 = {1->1, 2->2, 3->3, 4->4};

permsL4C4 = {{1->1,2->2,3->3,4->4},{1->1,2->4,3->3,4->2},{1->1,2->2,3->4,4->3},{1->1,2->3,3->4,4->2},{1->4,2->2,3->3,4->1},{1->4,2->1,3->3,4->2},{1->2,2->1,3->4,4->3},{1->2,2->3,3->4,4->1},{1->2,2->1,3->3,4->4},{1->2,2->4,3->3,4->1},{1->1,2->3,3->2,4->4},{1->1,2->4,3->2,4->3}};
netsL4C4 = Table[netToEdgeSplits[4,Join[Sort[edgesL4C4/.permsL4C4[[i]]],leafSet4],{5,6}],{i,1,12}];

edgesL4C3 = {1->5,2->5,3->7,4->8,5->6,6->7,7->8,8->6};

permsL4C3Ret67 = {{1->1,2->2,3->3,4->4},{1->1,2->3,3->2,4->4},{1->1,2->4,3->3,4->2},{1->3,2->2,3->1,4->4},{1->4,2->2,3->3,4->1},{1->3,2->4,3->1,4->2}};
netsL4C3Ret67 = Table[netToEdgeSplits[4,Join[Sort[edgesL4C3/.permsL4C3Ret67[[i]]],leafSet4],{6,7}],{i,1,6}];
permsL4C3Ret68 = {{1->1,2->2,3->3,4->4},{1->1,2->3,3->2,4->4}, {1->1,2->4, 3->3, 4->2}};
netsL4C3Ret68 =Table[netToEdgeSplits[4,Join[Sort[edgesL4C3/.permsL4C3Ret68[[i]]],leafSet4],{6,8}],{i,1,3}];
netsL4C3 = Join[netsL4C3Ret67, netsL4C3Ret68];


leafSet5 = {1->1,2->2,3->3,4->4,5->5};
edgesL5C3Top1 = {1->6,2->6,3->7,4->7,5->8,6->10,7->9,8->9,9->10,10->8};

permsL5C3Top1Ret810=Map[permToRules,{{1,2,3,4,5},{1,2,3,5,4},{1,2,4,5,3},{1,3,2,4,5},{1,3,2,5,4},{1,4,2,3,5},{1,5,2,3,4},{1,4,2,5,3},{1,5,2,4,3},{1,3,4,5,2},{1,4,3,5,2},{1,5,3,4,2},{2,3,4,5,1},{2,4,3,5,1},{2,5,3,4,1}}];
permsL5C3Top1Ret89 = Map[permToRules,{{1,2,3,4,5},{1,2,3,5,4},{1,2,4,5,3},{1,3,2,4,5},{1,3,2,5,4},{1,4,2,3,5},{1,5,2,3,4},{1,4,2,5,3},{1,5,2,4,3},{1,3,4,5,2},{1,4,3,5,2},{1,5,3,4,2},{2,3,1,4,5},{2,3,1,5,4},{2,4,1,3,5},{2,5,1,3,4},{2,4,1,5,3},{2,5,1,4,3},{3,4,1,2,5},{3,5,1,2,4},{4,5,1,2,3},{3,4,1,5,2},{3,5,1,4,2},{4,5,1,3,2},{2,3,4,5,1},{2,4,3,5,1},{2,5,3,4,1},{3,4,2,5,1},{3,5,2,4,1},{4,5,2,3,1}}];

netsL5C3Top1Ret810 = Table[netToEdgeSplits[5,Join[Sort[edgesL5C3Top1 /. permsL5C3Top1Ret810[[i]]],leafSet5],{8,10}],{i,1,Length[permsL5C3Top1Ret810]}];
netsL5C3Top1Ret89 = Table[netToEdgeSplits[5,Join[Sort[edgesL5C3Top1 /. permsL5C3Top1Ret89[[i]]],leafSet5],{8,9}],{i,1,Length[permsL5C3Top1Ret89]}];

edgesL5C3Top2 = {1->6,2->6,3->7,4->10,5->9,6->7,7->8,8->9,9->10,10->8};
permsL5C3Top2Ret810= Map[permToRules,{{1,2,3,4,5},{1,2,4,3,5},{1,2,5,3,4},{1,3,2,4,5},{1,4,2,3,5},{1,5,2,3,4},{1,3,4,2,5},{1,3,5,2,4},{1,4,3,2,5},{1,5,3,2,4},{1,4,5,2,3},{1,5,4,2,3},{2,3,1,4,5},{2,4,1,3,5},{2,5,1,3,4},{3,4,1,2,5},{3,5,1,2,4},{4,5,1,2,3},{2,3,4,1,5},{2,3,5,1,4},{2,4,3,1,5},{2,5,3,1,4},{2,4,5,1,3},{2,5,4,1,3},{3,4,2,1,5},{3,5,2,1,4},{4,5,2,1,3},{3,4,5,1,2},{3,5,4,1,2},{4,5,3,1,2}}];
permsL5C3Top2Ret89 = Map[permToRules,{{1,2,3,4,5},{1,2,3,5,4},{1,2,4,3,5},{1,2,5,3,4},{1,2,4,5,3},{1,2,5,4,3},{1,3,2,4,5},{1,3,2,5,4},{1,4,2,3,5},{1,5,2,3,4},{1,4,2,5,3},{1,5,2,4,3},{1,3,4,2,5},{1,3,5,2,4},{1,4,3,2,5},{1,5,3,2,4},{1,4,5,2,3},{1,5,4,2,3},{1,3,4,5,2},{1,3,5,4,2},{1,4,3,5,2},{1,5,3,4,2},{1,4,5,3,2},{1,5,4,3,2},{2,3,1,4,5},{2,3,1,5,4},{2,4,1,3,5},{2,5,1,3,4},{2,4,1,5,3},{2,5,1,4,3},{3,4,1,2,5},{3,5,1,2,4},{4,5,1,2,3},{3,4,1,5,2},{3,5,1,4,2},{4,5,1,3,2},{2,3,4,1,5},{2,3,5,1,4},{2,4,3,1,5},{2,5,3,1,4},{2,4,5,1,3},{2,5,4,1,3},{3,4,2,1,5},{3,5,2,1,4},{4,5,2,1,3},{3,4,5,1,2},{3,5,4,1,2},{4,5,3,1,2},{2,3,4,5,1},{2,3,5,4,1},{2,4,3,5,1},{2,5,3,4,1},{2,4,5,3,1},{2,5,4,3,1},{3,4,2,5,1},{3,5,2,4,1},{4,5,2,3,1},{3,4,5,2,1},{3,5,4,2,1},{4,5,3,2,1}}];

netsL5C3Top2Ret810 = Table[netToEdgeSplits[5,Join[Sort[edgesL5C3Top2 /. permsL5C3Top2Ret810[[i]]],leafSet5],{8,10}],{i,1,Length[permsL5C3Top2Ret810]}];
netsL5C3Top2Ret89 = Table[netToEdgeSplits[5,Join[Sort[edgesL5C3Top2 /. permsL5C3Top2Ret89[[i]]],leafSet5],{8,9}],{i,1,Length[permsL5C3Top2Ret89]}];

edgesL5C4 = {1->6,2->6,3->8,4->9,5->10,6->7,7->8,8->9,9->10,10->7};

permsL5C4Ret710 = Map[permToRules,{{1,2,3,4,5},{1,2,3,5,4},{1,2,4,3,5},{1,3,2,4,5},{1,3,2,5,4},{1,4,2,3,5},{1,5,2,3,4},{1,4,2,5,3},{1,5,2,4,3},{1,3,4,2,5},{1,4,3,2,5},{1,5,3,2,4},{2,3,1,4,5},{2,3,1,5,4},{2,4,1,3,5},{2,5,1,3,4},{2,4,1,5,3},{2,5,1,4,3},{3,4,1,2,5},{3,5,1,2,4},{4,5,1,2,3},{3,4,1,5,2},{3,5,1,4,2},{4,5,1,3,2},{2,3,4,1,5},{2,4,3,1,5},{2,5,3,1,4},{3,4,2,1,5},{3,5,2,1,4},{4,5,2,1,3}}];
permsL5C4Ret89 =  Map[permToRules,{{1,2,3,4,5},{1,2,3,5,4},{1,2,4,3,5},{1,3,2,4,5},{1,3,2,5,4},{1,4,2,3,5},{1,5,2,3,4},{1,4,2,5,3},{1,5,2,4,3},{1,3,4,2,5},{1,4,3,2,5},{1,5,3,2,4},{2,3,1,4,5},{2,3,1,5,4},{2,4,1,3,5},{2,5,1,3,4},{2,4,1,5,3},{2,5,1,4,3},{3,4,1,2,5},{3,5,1,2,4},{4,5,1,2,3},{3,4,1,5,2},{3,5,1,4,2},{4,5,1,3,2},{2,3,4,1,5},{2,4,3,1,5},{2,5,3,1,4},{3,4,2,1,5},{3,5,2,1,4},{4,5,2,1,3}}];
permsL5C4Ret910 = Map[permToRules,{{1,2,3,4,5},{1,2,3,5,4},{1,2,4,3,5},{1,2,5,3,4},{1,2,4,5,3},{1,2,5,4,3},{1,3,2,4,5},{1,3,2,5,4},{1,4,2,3,5},{1,5,2,3,4},{1,4,2,5,3},{1,5,2,4,3},{1,3,4,2,5},{1,3,5,2,4},{1,4,3,2,5},{1,5,3,2,4},{1,4,5,2,3},{1,5,4,2,3},{1,3,4,5,2},{1,3,5,4,2},{1,4,3,5,2},{1,5,3,4,2},{1,4,5,3,2},{1,5,4,3,2},{2,3,1,4,5},{2,3,1,5,4},{2,4,1,3,5},{2,5,1,3,4},{2,4,1,5,3},{2,5,1,4,3},{3,4,1,2,5},{3,5,1,2,4},{4,5,1,2,3},{3,4,1,5,2},{3,5,1,4,2},{4,5,1,3,2},{2,3,4,1,5},{2,3,5,1,4},{2,4,3,1,5},{2,5,3,1,4},{2,4,5,1,3},{2,5,4,1,3},{3,4,2,1,5},{3,5,2,1,4},{4,5,2,1,3},{3,4,5,1,2},{3,5,4,1,2},{4,5,3,1,2},{2,3,4,5,1},{2,3,5,4,1},{2,4,3,5,1},{2,5,3,4,1},{2,4,5,3,1},{2,5,4,3,1},{3,4,2,5,1},{3,5,2,4,1},{4,5,2,3,1},{3,4,5,2,1},{3,5,4,2,1},{4,5,3,2,1}}];

netsL5C4Ret710 = Table[netToEdgeSplits[5,Join[Sort[edgesL5C4/. permsL5C4Ret710[[i]]],leafSet5],{7,10}],{i,1,Length[permsL5C4Ret710]}];
netsL5C4Ret89 = Table[netToEdgeSplits[5,Join[Sort[edgesL5C4/. permsL5C4Ret89 [[i]]],leafSet5],{8,9}],{i,1,Length[permsL5C4Ret89 ]}];
netsL5C4Ret910 = Table[netToEdgeSplits[5,Join[Sort[edgesL5C4/. permsL5C4Ret910 [[i]]],leafSet5],{9,10}],{i,1,Length[permsL5C4Ret910 ]}];


edgesL5C5 = {1->6,2->7,3->8,4->9,5->10,6->7,7->8,8->9,9->10,10->6};

permsL5C5 = Map[permToRules,{{1,2,3,4,5},{1,2,3,5,4},{1,2,4,5,3},{1,3,2,4,5},{1,3,2,5,4},{1,4,2,3,5},{2,1,3,4,5},{2,1,3,5,4},{2,1,4,5,3},{3,1,2,4,5},{3,1,2,5,4},{4,1,2,3,5},{5,1,2,3,4},{4,1,2,5,3},{5,1,2,4,3},{3,1,4,5,2},{4,1,3,5,2},{5,1,3,4,2},{2,3,1,4,5},{2,3,1,5,4},{2,4,1,3,5},{3,2,1,4,5},{3,2,1,5,4},{4,2,1,3,5},{5,2,1,3,4},{4,2,1,5,3},{5,2,1,4,3},{3,4,1,2,5},{4,3,1,2,5},{5,3,1,2,4}}];
netsL5C5 = Table[netToEdgeSplits[5,Join[Sort[edgesL5C5/. permsL5C5[[i]]],leafSet5],{6,10}],{i,1,Length[permsL5C5]}];
