
#####################################################################################
#####################################################################################
## COMPLETE THE FUNCTIONS BELOW
#####################################################################################
#####################################################################################

import numpy as np

def InitializeEquation(N_NODE,N_PRE_DISP,DISP_NODE, DOF_NODE, N_ELEM, NNODE_ELE, ELEM_NODE):

    # define your function HERE

    # Return EQ_NUM as 2D numpy array of intergers of appropriate shape.
    EQ_NUM = np.zeros((N_NODE, DOF_NODE), dtype=int)
    for i in range(N_PRE_DISP):
        NODE = DISP_NODE[i,0]-1
        IDOF = DISP_NODE[i,1]-1
        EQ_NUM[NODE,IDOF]= -(i+1)
        
    ROW=0
    
    for i in range(N_NODE):
        for j in range(DOF_NODE):
            if EQ_NUM[i,j]==0:
                ROW+=1
                EQ_NUM[i,j]=ROW

    EDOF=4

    LM = np.zeros((N_ELEM, EDOF), dtype=int)

    for i in range(N_ELEM):
        for j in range(NNODE_ELE):
            NODE = ELEM_NODE[i,j]-1
            for k in range(DOF_NODE):
                eq=EQ_NUM[NODE,k]
                col_idx = k +(j)*DOF_NODE
                LM[i,col_idx]=eq




    # Return N_DOF, TOTAL_DOF, EDOF as scalars.

    N_DOF = ROW
    TOTAL_DOF = 2* N_ELEM
    # Return LM as 2D numpy array of integers of appropriate shape.

    

    return EQ_NUM, N_DOF, TOTAL_DOF, EDOF, LM


def AssemblePrescribedDisplacement(N_PRE_DISP, DISP_NODE, DISP_VAL, EQ_NUM):

    # Initialize arrays that need to be returned (UP)
    # Return UP as 1D numpy array of appropriate shape.

    # define your function HERE
    
    UP= np.zeros(N_PRE_DISP)
    
    for i in range(N_PRE_DISP):  # Python uses 0-based indexing
        node = DISP_NODE[i, 0] - 1  # Adjust for 0-based indexing (node number)
        idof = DISP_NODE[i, 1] - 1  # Adjust for 0-based indexing (degree of freedom, 1 or 2)
        u = DISP_VAL[i]  # Prescribed displacement value
        row = -EQ_NUM[node, idof] - 1  # Adjust for 0-based indexing, assuming EQ_NUM is negative-based
        UP[row] = u  # Assign the displacement value to the correct position


    return UP


def AssembleConcentratedNodalForces(N_LOAD, FORCE_NODE, FORCE_VAL, EQ_NUM, N_DOF):

    # Initialize arrays that need to be returned (PF)
    # Return PF as 1D numpy array of appropriate shape.

    # define your function HERE
    
    PF = np.zeros(N_DOF)
    for i in range(N_LOAD):
        node = FORCE_NODE[i,0]-1
        idof = FORCE_NODE[i,1]-1
        f = FORCE_VAL[i]
        row = EQ_NUM[node, idof] - 1 
        if row>=0:
            PF[row]+=f
            
        
    

    return PF


def ElementArrays(I_ELEM, ELEM_STIFF, COORDS, ELEM_NODE, ELEM_LOAD, ELEM_AREA):

    # define your function HERE
    # Return KEL as 2D numpy array of appropriate shape.
    node_a, node_b = ELEM_NODE[I_ELEM]  # This assumes ELEM_NODE is zero-indexed and I_ELEM is given
    
    E = ELEM_STIFF[I_ELEM,0]
    A = ELEM_AREA[I_ELEM]
    
    coord_a = COORDS[node_a-1]
    coord_b = COORDS[node_b-1]
    
    a = coord_b[0]- coord_a[0]
    b = coord_b[1]- coord_a[1]

    
    l = np.sqrt(a**2 + b**2)
    
    c = a/l
    s = b/l
   
    KEL = (E * A / l) * np.array([[ c**2,  c*s, -c**2, -c*s],
                                   [ c*s,   s**2, -c*s,  -s**2],
                                   [-c**2, -c*s,  c**2,  c*s],
                                   [-c*s,  -s**2, c*s,   s**2]])
    
        
    # Return PEL as 1D numpy array of appropriate shape.
    
    PG1 = ELEM_LOAD[I_ELEM,0]
    PG = (A*l)* np.array([0.0, -0.5 *PG1, 0.0, -0.5 * PG1])
    
    
    alpha = ELEM_STIFF[I_ELEM,1]
    deltaT = ELEM_LOAD[I_ELEM,1]
    
    PT = (E*A*alpha*deltaT)* np.array([-(c),
                                       -(s),  # The load is distributed along the y-direction
                                       (c),
                                       (s)])
    
    PEL = PG+PT
    return KEL,PEL


def AssembleGlobalArrays(KEL, PEL, ELEM_NODE, EQ_NUM, I_ELEM, PP, PF, KPP, KPF, KFF, KFP, NNODE_ELE, DOF_NODE, LM, EDOF):
    # Loop over the element degrees of freedom
    for i in range(EDOF):
        ROW = LM[I_ELEM, i]
        if ROW > 0:  # Free degree of freedom
            PF[ROW - 1] += PEL[i]  # Adjust for zero-based indexing
        else:  # Prescribed degree of freedom
            PP[abs(ROW) - 1] += PEL[i]  # Convert to positive index and adjust for zero-based indexing

        for j in range(EDOF):
            COL = LM[I_ELEM, j]
            if ROW > 0:  # Free row
                if COL > 0:  # Free column
                    KFF[ROW - 1, COL - 1] += KEL[i, j]  # Adjust for zero-based indexing
                else:  # Prescribed column
                    KFP[ROW - 1, abs(COL) - 1] += KEL[i, j]  # Convert to positive index and adjust for zero-based indexing
            else:  # Prescribed row
                if COL > 0:  # Free column
                    KPF[abs(ROW) - 1, COL - 1] += KEL[i, j]  # Convert to positive index and adjust for zero-based indexing
                else:  # Prescribed column
                    KPP[abs(ROW) - 1, abs(COL) - 1] += KEL[i, j]  # Convert to 

    return PP, PF, KPP, KPF, KFF, KFP



def AssembleModule(N_ELEM, N_LOAD, N_PRE_DISP, ELEM_NODE, ELEM_STIFF, FORCE_NODE, FORCE_VAL, DISP_NODE, DISP_VAL, EQ_NUM, N_DOF, NNODE_ELE, COORDS, ELEM_AREA, ELEM_LOAD, DOF_NODE, LM, EDOF):

    ## Initialize arrays th)at need to be returned (KPP, KPF, KFF, KFP, PP)
    # KPP, KPF, KFF, KFP are 2D numpy arrays of appropriate shapes.
    KPP = np.zeros((N_PRE_DISP, N_PRE_DISP))

    KFF = np.zeros((N_DOF, N_DOF))
    
    KPF = np.zeros((N_PRE_DISP,N_DOF))
    
    KFP = np.zeros((N_DOF,N_PRE_DISP))
    
    # PP is a 1D numpy array of appropriate shape.
    PP = np.zeros((N_PRE_DISP))

    ## Prescribe boundary conditions
    UP = AssemblePrescribedDisplacement(N_PRE_DISP, DISP_NODE, DISP_VAL, EQ_NUM)

    ## apply concentrated forces
    PF = AssembleConcentratedNodalForces(N_LOAD, FORCE_NODE, FORCE_VAL, EQ_NUM, N_DOF)

    ## assemble global stiffness and force
    for I_ELEM in range(N_ELEM):

        ##get element stiffness matrix and load vector
         KEL,PEL = ElementArrays(I_ELEM, ELEM_STIFF, COORDS, ELEM_NODE, ELEM_LOAD, ELEM_AREA)

        ##assemble global stiffness and load vector
         PP,PF,KPP,KPF,KFF,KFP=AssembleGlobalArrays(KEL, PEL, ELEM_NODE, EQ_NUM, I_ELEM, PP, PF, KPP, KPF, KFF, KFP, NNODE_ELE, DOF_NODE, LM, EDOF)

    return KPP,KPF,KFF,KFP,PF,PP,UP

def SolveModule(N_NODE, KPP, KPF, KFF, KFP, PF, UP, EQ_NUM, DOF_NODE, PP):
    # Solve for UF: the unknown displacements
    # We use np.linalg.solve to solve the system KFF * UF = PF - KFP * UP
    UF = np.linalg.solve(KFF, PF - np.dot(KFP, UP))
    
    # Solve for R (reactions):
    # R = KPF * UF + KPP * UP - PP
    R = np.dot(KPF, UF) + np.dot(KPP, UP) - PP
    
    # The function AssembleGlobalResponses is assumed to be implemented.
    # It will compile the global displacements into a single array.
    UUR = AssembleGlobalResponses(EQ_NUM, UF, UP)
    
    return UUR, UF, R

def AssembleGlobalResponses(EQ_NUM, UF, UP):
    # Initialize the global displacement vector with zeros
    UUR = np.zeros((EQ_NUM.shape[0], EQ_NUM.shape[1]))
    
    # Loop through each node and degree of freedom
    for i in range(EQ_NUM.shape[0]):
        for j in range(EQ_NUM.shape[1]):
            eq_num = EQ_NUM[i, j]
            if eq_num > 0:  # Free DOF
                # -1 for zero-based indexing
                UUR[i, j] = UF[eq_num - 1]
            elif eq_num < 0:  # Prescribed DOF
                # -1 for zero-based indexing and negative sign for prescribed values
                UUR[i, j] = UP[-eq_num - 1]
    
    # Flatten UUR to 1D if necessary or keep it 2D based on how you want to return it
    # UUR = UUR.flatten()
    
    return UUR


def GetStrainStress(N_ELEM, NNODE_ELE, NDIM, ELEM_NODE, ELEM_STIFF, ELEM_LOAD, COORDS, UUR):

    # define your function HERE
    # Return strain, stress as 1D numpy arrays of appropriate shapes.

    strain = np.zeros(N_ELEM)
    stress = np.zeros(N_ELEM)

    for i in range(N_ELEM):
        # Get the node numbers for the current element
        node_a, node_b = ELEM_NODE[i]
        
        # Get coordinates for each node
        x1, y1 = COORDS[node_a - 1]  # Adjust for zero-based indexing
        x2, y2 = COORDS[node_b - 1]  # Adjust for zero-based indexing
        
        # Calculate the element's length and directional cosines
        a = x2 - x1
        b = y2 - y1
        l = np.sqrt(a**2 + b**2)
        c = a / l
        s = b / l
        
        # Transformation matrix
        T = np.array([[c, s, 0, 0],
                      [0, 0, c, s]])
        
        # Local displacements
        Ua = UUR[node_a - 1]  # Adjust for zero-based indexing
        Ub = UUR[node_b - 1]  # Adjust for zero-based indexing
        UE = np.dot(T, np.concatenate((Ua, Ub)))
        
        # Strain and stress calculation
        epsilon = (UE[1] - UE[0]) / l
        alpha = ELEM_STIFF[i, 1]
        deltaT = ELEM_LOAD[i, 1]
        E = ELEM_STIFF[i, 0]
        
        # Save computed strain and stress
        strain[i] = epsilon
        stress[i] = E * (epsilon - alpha * deltaT)
    
    return strain, stress

#####################################################################################
#####################################################################################
## DO NOT MODIFY THE REMAINING OF THIS FILE
## The following functions are provided to you - do not make any changes
#####################################################################################
#####################################################################################

import numpy as np
import csv
import matplotlib.pyplot as plt

def ReadInput(filename, ndim, nnode_ele):
# % INPUT:
# % filename: name of the file with input values
# % ndim: problem dimension
# % nnode_ele: number of nodes per element
#
# % OUTPUT:
# % N_NODE, the (integer) number of nodes
# % N_ELEM, the (integer) number of elements
# % N_LOAD, the (integer) number of nonzero nodal forces, i.e. Pf
# % N_PRE_DISP, the (integer) number of nodes with prescribed displacement, i.e. Up
# % COORDS(N_NODE, ndim), matrix that has the nodal coordinates
# % ELEM_NODE(N_ELEM, nnode_ele), matrix that contains node numbers (integers) for each element
# % ELEM_STIFF(N_ELEM,2), 2d array that contains stiffness value (real) and thermal expansion coefficient for each element
# % ELEM_AREA(N_ELEM,), 1d array that contains cross sectional area for each element
# % ELEM_LOAD(N_ELEM,2), 2d array that contains specific weight (real) and change in temperature for each element
# % FORCE_NODE(N_LOAD,2), 2d array that contains the node numbers (integer) where forces are applied and the direction
# % FORCE_VAL(N_LOAD,), 1d array that contains the value of the forces (real) corresponding to the numbers in FORCE_NODE
# % DISP_NODE(N_PRE_DISP,2), 2d array that contains the node numbers (integer) where boundary conditions are imposed and the direction
# % DISP_VAL(N_PRE_DISP,),  1d array that contains the value of the prescribed displacements (real) corresponding to the numbers in DISP_NODE

    with open(filename, newline = '') as fin:
        mesh_data = csv.reader(fin, delimiter=',')
        line = next(mesh_data)
        [n_node, n_elem, n_load, n_pre_disp] = list(map(int, line))

        # Get coordinates
        coords = np.zeros((n_node,ndim))
        for i in range(n_node):
            line = next(mesh_data)
            coords[i] = list(map(float, line))

        # Get element connectivity
        elem_node = np.zeros((n_elem,nnode_ele),dtype=int)
        # stores the Elasticity modulus and thermal expansion coefficient for each truss member
        elem_stiff = np.zeros((n_elem,2))
        # cross sectional area of the truss element
        elem_area = np.zeros(n_elem)
        # this will store the specific weight (gamma) of the material, used to calculate the load due to gravity,
        # and the change in temperature, for the thermal loading
        elem_load = np.zeros((n_elem,2))
        for i in range(n_elem):
            line = next(mesh_data)
            line_list = list(map(float, line))
            elem_node[i] = line_list[:2]
            elem_stiff[i,0] = line_list[2]
            elem_load[i,0] = line_list[3]
            elem_area[i] = line_list[4]
            elem_stiff[i,1] = line_list[5]
            elem_load[i,1] = line_list[6]

        # Get external forces at nodes
        force_node = np.zeros((n_load,2), dtype=int)
        force_val = np.zeros(n_load)
        for i in range(n_load):
            line = next(mesh_data)
            line_list = list(map(float, line))
            force_node[i] = line_list[:2]
            force_val[i] = line_list[2]


        # Get boundary conditions
        disp_node = np.zeros((n_pre_disp,2), dtype=int)
        disp_val = np.zeros(n_pre_disp)
        for i in range(n_pre_disp):
            line = next(mesh_data)
            line_list = list(map(float, line))
            disp_node[i] = line_list[:2]
            disp_val[i] = line_list[2]


    return(n_node, n_elem, n_load, n_pre_disp,coords,elem_node, elem_stiff, elem_area, elem_load, force_node, force_val, disp_node, disp_val)


def AssembleGlobalResponses(eq_num, uf, up):

    uur = np.zeros(eq_num.shape)
    for i,val in enumerate(eq_num):
        for j, row in enumerate(val):
            if row > 0:
                uur[i,j]=uf[row-1]
            else:
                uur[i,j]=up[-row-1]
    return uur


def plotTruss(uur, n_elem, elem_node, coords):

    plottruss = plt.figure()
    bar_lengthv = np.zeros(n_elem)
    for iel in range(n_elem):
        xel = coords[elem_node[iel,:]-1,:]
        a = xel[1,0] - xel[0,0]
        b = xel[1,1] - xel[0,1]
        bar_lengthv[iel] = np.sqrt( a**2 +  b**2 )
    bar_length_max=np.amax(bar_lengthv)

    disp_max=np.sum(np.amax(np.abs(uur)))
    ampli=0.2*bar_length_max/disp_max
    plottitle = ' undeformed and deformed shapes \n disp. amplification factor =' + str(ampli)
    plt.axes().set_aspect('equal')
    plt.title(plottitle)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')

    for iel in range(n_elem):
        xel_undef = coords[elem_node[iel,:]-1,:]
        uel = uur[elem_node[iel,:]-1,:]
        xel_def = xel_undef + ampli*uel
        plt.plot(xel_undef[:,0],xel_undef[:,1],'bo-',xel_def[:,0],xel_def[:,1],'ro--',linewidth=2)
    plt.show()

    return(plottruss)
