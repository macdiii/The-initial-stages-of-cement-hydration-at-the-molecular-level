#密度重写嘤嘤嘤

import numpy as np
import pandas as pd
def data_init(natoms,dumplines,i):
    atom=np.zeros((natoms,3))   
    idtype=np.zeros((natoms,2))    
    xlo,xhi_boundary=dumplines[i+5].split()                                                                
    ylo,yhi=dumplines[i+6].split()                                                                            
    xlo,xhi_boundary,ylo,yhi=float(xlo),float(xhi_boundary),float(ylo),float(yhi)
    xhi=xhi_boundary
    for j in range(natoms) :
        line0=dumplines[i+9+j].split() 
        ID = int(line0[3])
        idtype[ID-1][0]=ID 
        idtype[ID-1][1]=int(line0[4]) 
        atom[ID-1][0]=float(line0[0]) 
        atom[ID-1][1]=float(line0[1]) 
        atom[ID-1][2]=float(line0[2]) 
    return atom,idtype,  xlo,xhi_boundary,ylo,yhi,xhi
def huge_atoms(atom,xlo,xhi_boundary,ylo,yhi,xhi):
    a1=atom+[-xhi+xlo,yhi-ylo,0]                    
    a2=atom+[0,yhi-ylo,0]                      
    a3=atom+[xhi-xlo,yhi-ylo,0]                     
    a4=atom+[-xhi+xlo,0,0]
    a5=atom+[0,0,0]
    a6=atom+[xhi-xlo,0,0]
    a7=atom+[-xhi+xlo,-yhi+ylo,0]
    a8=atom+[0,-yhi+ylo,0] 
    a9=atom+[xhi-xlo,-yhi+ylo,0]
    ATOMS=np.concatenate((a1,a2),axis=0)
    ATOMS=np.concatenate((ATOMS,a3),axis=0)
    ATOMS=np.concatenate((ATOMS,a4),axis=0)
    ATOMS=np.concatenate((ATOMS,a5),axis=0)
    ATOMS=np.concatenate((ATOMS,a6),axis=0)
    ATOMS=np.concatenate((ATOMS,a7),axis=0)
    ATOMS=np.concatenate((ATOMS,a8),axis=0)
    ATOMS=np.concatenate((ATOMS,a9),axis=0)
    return ATOMS# ,a1,a2,a3,a4,a5,a6,a7,a8,a9
def huge_idtype(idtype):
    IDTYPE=idtype
    for i in range(8):            
        IDTYPE=np.concatenate((IDTYPE,idtype),axis=0)
    return IDTYPE
def if_surface(j,atom,surfacelo,surfacehi):
    if atom[j][2] < surfacelo or atom[j][2] > surfacehi :
        return False
    else:
        return True   
def if_soild(j,idtype):
    if (idtype[j,0]-1)%1472 >= 560 :
        return False
    else:
        return True   
def if_type(j,idtype,target_type):
    if idtype[j][1] == target_type:
        return True
    else:
        return False  
def distence(x1,y1,z1,x2,y2,z2):
    result=np.sqrt(np.square(x1-x2)+np.square(y1-y2)+np.square(z1-z2))
    return result
def position(atom,j):
    x=atom[j][0]
    y=atom[j][1]
    z=atom[j][2]
    return x,y,z 
def boundry(x,y,z,value):
    x_L=x-value
    x_R=x+value
    y_L=y-value
    y_R=y+value
    z_L=z-value
    z_R=z+value
    return x_L,x_R,y_L,y_R,z_L,z_R
def judge_boumdry(x,y,z,x_L,x_R,y_L,y_R,z_L,z_R):
    if x>=x_L and x<=x_R and y>=y_L and y<=y_R and z>=z_L and z<=z_R :
        return True
    else:
        return False
def cut_atom_area(distance,xlo,xhi_boundary,yhi,ylo,surfacelo,surfacehi):                
    L_boundary=xlo-distance                              
    R_boundary=xhi_boundary+distance                   
    U_boundary=yhi+distance                            
    D_boundary=ylo-distance                            
    U_Z_boundary=surfacehi+distance                      
    D_Z_boundary=surfacelo-distance                     
    return L_boundary,R_boundary,U_boundary,D_boundary,U_Z_boundary,D_Z_boundary
def judge_cut_area(big_atom,big_idtype,L_boundary,R_boundary,U_boundary,D_boundary,U_Z_boundary,D_Z_boundary):     
    list_delete=[]
    
    for i in range(big_atom.shape[0]):
        if big_atom[i][0]<L_boundary or big_atom[i][0]>R_boundary or big_atom[i][1]<D_boundary or big_atom[i][1]>U_boundary or big_atom[i][2]<D_Z_boundary or big_atom[i][2]>U_Z_boundary:
            list_delete.append(i)    
    atom_update=np.delete(big_atom,list_delete,0)
    idtype_update=np.delete(big_idtype,list_delete,0)
    sizeee=atom_update.shape[0]
    return atom_update , idtype_update ,sizeee
Ca_O=2.83
H_O=1.2
O_Si=2
surfacelo=130 
surfacehi=150
distance=10
n_file=1
lo=-5      
hi=305     

for p in range(1,n_file+1):

    with open("dump.trj","r") as f1:             
        dumplines = f1.readlines()
    natoms=int(dumplines[3].split()[0])
    zhen=0
    result=np.zeros ( (500,6) )
    for i in range(0,len(dumplines),natoms+9) :
        print("zhen",zhen)
        time=int(dumplines[i+1].split()[0]) 
        idtype=np.zeros((natoms,2))
        atom=np.zeros((natoms,3))
        time_step=int(dumplines[i+1].split()[0])
        atom,idtype,  xlo,xhi_boundary,ylo,yhi,xhi=data_init(natoms,dumplines,i)                                
        big_atom=huge_atoms(atom,xlo,xhi_boundary,ylo,yhi,xhi)
        big_idtype=huge_idtype(idtype)
        L_boundary,R_boundary,U_boundary,D_boundary,U_Z_boundary,D_Z_boundary=cut_atom_area(distance,xlo,xhi_boundary,yhi,ylo,surfacelo,surfacehi)########
        U_Z_boundary,D_Z_boundary=350,-10
        big_atom,big_idtype ,sizeee=judge_cut_area(big_atom,big_idtype,L_boundary,R_boundary,U_boundary,D_boundary,U_Z_boundary,D_Z_boundary)######
        S=(xhi-xlo)*(yhi-ylo)
        part=np.zeros ( (500,6) )
        for j in range(natoms):
            hight=atom[j,2]
            if if_type(j,idtype,1):
                part[int(hight)+10,1]+=1
                
            elif if_type(j,idtype,2):
                part[int(hight)+10,2]+=1
            elif if_type(j,idtype,3):
                if if_soild(j,idtype):
                    part[int(hight)+10,3]+=1
                else:
                    part[int(hight)+10,4]+=1
            elif if_type(j,idtype,4):
                part[int(hight)+10,5]+=1
        result+=part/S
    result=result/11
    
    df=pd.DataFrame(result)
    df.to_csv("atom.csv")
