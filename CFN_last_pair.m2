-- Coordinates for the first variety
qVars = {q000000, q000011, q000101, q000110, q001001, q001010, q001100, q001111, q010001, q010010, q010100, q010111, q011000, q011011, q011101, q011110, q100001, q100010, 
 q100100, q100111, q101000, q101011, q101101, q101110, q110000, q110011, q110101, q110110, q111001, q111010, q111100, q111111};

-- Coordinates for the second variety
pVars = {p000000, p000011, p000101, p000110, p001001, p001010, p001100, p001111, p010001, p010010, p010100, p010111, p011000, p011011, p011101, p011110, p100001, p100010, 
 p100100, p100111, p101000, p101011, p101101, p101110, p110000, p110011, p110101, p110110, p111001, p111010, p111100, p111111};

-- Coordinates for the join varieties
rVars = {r000000, r000011, r000101, r000110, r001001, r001010, r001100, r001111, r010001, r010010, r010100, r010111, r011000, r011011, r011101, r011110, r100001, r100010, 
 r100100, r100111, r101000, r101011, r101101, r101110, r110000, r110011, r110101, r110110, r111001, r111010, r111100, r111111};


-- Ring For Elimination Ideal
R = QQ[(qVars|pVars|rVars), MonomialOrder => Eliminate 64]

-- Matrices corresponding to the splits in the trees.

qSplitMatrix1230 = matrix {{q000000, q000011, q000101, q000110}, {q011000, q011011, q011101, q011110}, {q101000, q101011, q101101, q101110}, {q110000, q110011, q110101, q110110}}
qSplitMatrix1231 = matrix {{q001001, q001010, q001100, q001111}, {q010001, q010010, q010100, q010111}, {q100001, q100010, q100100, q100111}, {q111001, q111010, q111100, q111111}}

pSplitMatrix1230 = matrix {{p000000, p000011, p000101, p000110}, {p011000, p011011, p011101, p011110}, {p101000, p101011, p101101, p101110}, {p110000, p110011, p110101, p110110}}
pSplitMatrix1231 = matrix {{p001001, p001010, p001100, p001111}, {p010001, p010010, p010100, p010111}, {p100001, p100010, p100100, p100111}, {p111001, p111010, p111100, p111111}}

qSplitMatrix120 = matrix {{q000000, q000011, q000101, q000110, q001001, q001010, q001100, q001111}, {q110000, q110011, q110101, q110110, q111001, q111010, q111100, q111111}}
qSplitMatrix121 = matrix {{q010001, q010010, q010100, q010111, q011000, q011011, q011101, q011110}, {q100001, q100010, q100100, q100111, q101000, q101011, q101101, q101110}}

pSplitMatrix120 = matrix {{p000000, p000011, p000101, p000110, p001001, p001010, p001100, p001111}, {p110000, p110011, p110101, p110110, p111001, p111010, p111100, p111111}}
pSplitMatrix121 = matrix {{p010001, p010010, p010100, p010111, p011000, p011011, p011101, p011110}, {p100001, p100010, p100100, p100111, p101000, p101011, p101101, p101110}}


qSplitMatrix12340 = matrix {{q000000, q000011}, {q001100, q001111}, {q010100, q010111}, {q011000, q011011}, {q100100, q100111}, {q101000, q101011}, {q110000, q110011}, {q111100, q111111}}
qSplitMatrix12341 = matrix {{q000101, q000110}, {q001001, q001010}, {q010001, q010010}, {q011101, q011110}, {q100001, q100010}, {q101101, q101110}, {q110101, q110110}, {q111001, q111010}}

pSplitMatrix12340 = matrix {{p000000, p000011}, {p001100, p001111}, {p010100, p010111}, {p011000, p011011}, {p100100, p100111}, {p101000, p101011}, {p110000, p110011}, {p111100, p111111}}
pSplitMatrix12341 = matrix {{p000101, p000110}, {p001001, p001010}, {p010001, p010010}, {p011101, p011110}, {p100001, p100010}, {p101101, p101110}, {p110101, p110110}, {p111001, p111010}}

qSplitMatrix230 = matrix {{q000000, q000011, q000101, q000110, q100001, q100010, q100100, q100111}, {q011000, q011011, q011101, q011110, q111001, q111010, q111100, q111111}}
qSplitMatrix231 = matrix {{q001001, q001010, q001100, q001111, q101000, q101011, q101101, q101110}, {q010001, q010010, q010100, q010111, q110000, q110011, q110101, q110110}}

pSplitMatrix230 = matrix {{p000000, p000011, p000101, p000110, p100001, p100010, p100100, p100111}, {p011000, p011011, p011101, p011110, p111001, p111010, p111100, p111111}}
pSplitMatrix231 = matrix {{p001001, p001010, p001100, p001111, p101000, p101011, p101101, p101110}, {p010001, p010010, p010100, p010111, p110000, p110011, p110101, p110110}}

qSplitMatrix12360 = matrix {{q000000, q000110}, {q001001, q001111}, {q010001, q010111}, {q011000, q011110}, {q100001, q100111}, {q101000, q101110}, {q110000, q110110}, {q111001, q111111}}
qSplitMatrix12361 = matrix {{q000011, q000101}, {q001010, q001100}, {q010010, q010100}, {q011011, q011101}, {q100010, q100100}, {q101011, q101101}, {q110011, q110101}, {q111010, q111100}}

pSplitMatrix12360 = matrix {{p000000, p000110}, {p001001, p001111}, {p010001, p010111}, {p011000, p011110}, {p100001, p100111}, {p101000, p101110}, {p110000, p110110}, {p111001, p111111}}
pSplitMatrix12361 = matrix {{p000011, p000101}, {p001010, p001100}, {p010010, p010100}, {p011011, p011101}, {p100010, p100100}, {p101011, p101101}, {p110011, p110101}, {p111010, p111100}}


-- Equations given by the map parameterizing the join, i.e. r = p + q
sumEqns = {p000000 + q000000 - r000000, p000011 + q000011 - r000011, p000101 + q000101 - r000101, p000110 + q000110 - r000110, p001001 + q001001 - r001001, p001010 + q001010 - r001010, p001100 + q001100 - r001100, 
 p001111 + q001111 - r001111, p010001 + q010001 - r010001, p010010 + q010010 - r010010, p010100 + q010100 - r010100, p010111 + q010111 - r010111, p011000 + q011000 - r011000, p011011 + q011011 - r011011, 
 p011101 + q011101 - r011101, p011110 + q011110 - r011110, p100001 + q100001 - r100001, p100010 + q100010 - r100010, p100100 + q100100 - r100100, p100111 + q100111 - r100111, p101000 + q101000 - r101000, 
 p101011 + q101011 - r101011, p101101 + q101101 - r101101, p101110 + q101110 - r101110, p110000 + q110000 - r110000, p110011 + q110011 - r110011, p110101 + q110101 - r110101, p110110 + q110110 - r110110, 
 p111001 + q111001 - r111001, p111010 + q111010 - r111010, p111100 + q111100 - r111100, p111111 + q111111 - r111111}




-- Create elimination ideal for the first pair {T_1,T_2} where T_1 corresponds to the q coords and T_2 to the p coords. Projecting onto the r coords gives us the vanishing ideal of the join

-- ideal of eqns we get from the splits of T_1 = {{1,2,3},{1,2},{1,2,3,4}}
I_1 = minors(2, qSplitMatrix1230) + minors(2, qSplitMatrix1231) + minors(2, qSplitMatrix120) + minors(2, qSplitMatrix121) + minors(2, qSplitMatrix12340) + minors(2, qSplitMatrix12341)

-- ideal of eqns we get from the splits of T_2 = {{1,2,3},{2,3},{1,2,3,6}}
I_2 = minors(2, pSplitMatrix1230) + minors(2, pSplitMatrix1231) + minors(2, pSplitMatrix230) + minors(2, pSplitMatrix231) + minors(2, pSplitMatrix12360) + minors(2, pSplitMatrix12361)

elimIdealT = I_1 + I_2 + ideal(sumEqns)

gbT = selectInSubring(1, gens gb(elimIdealT, DegreeLimit => 4))

I = ideal(gbT)


-- Create elimination ideal for the second pair {S_1,S_2} where S_1 corresponds to the q coords and S_2 to the p coords. Projecting onto the r coords gives us the vanishing ideal of the join

-- ideal of eqns we get from the splits of S_1 = {{1,2,3},{1,2},{1,2,3,6}}
J_1 = minors(2, qSplitMatrix1230) + minors(2, qSplitMatrix1231) + minors(2, qSplitMatrix120) + minors(2, qSplitMatrix121) + minors(2, qSplitMatrix12360) + minors(2, qSplitMatrix12361)

-- ideal of eqns we get from the splits of S_2 = {{1,2,3},{2,3},{1,2,3,4}}
J_2 = minors(2, pSplitMatrix1230) + minors(2, pSplitMatrix1231) + minors(2, pSplitMatrix230) + minors(2, pSplitMatrix231) + minors(2, pSplitMatrix12340) + minors(2, pSplitMatrix12341)

elimIdealS = J_1 + J_2 + ideal(sumEqns)

gbS = selectInSubring(1, gens gb(elimIdealS, DegreeLimit => 4))

J = ideal(gbS)


I == J






















