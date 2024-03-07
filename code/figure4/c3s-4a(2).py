import numpy as np
import pandas as pd
def data_init(natoms,dumplines,i):
    atom=np.zeros((natoms,3))   
    idtype=np.zeros((natoms,2))
    xlo,xhi_boundary,xy=dumplines[i+5].split()                                                                     
    ylo,yhi,xz=dumplines[i+6].split()                                                                         
    xlo,xhi_boundary,xy,ylo,yhi,xz=float(xlo),float(xhi_boundary),float(xy),float(ylo),float(yhi),float(xz)
    xhi=xhi_boundary-xy
    for j in range(natoms) :
        line0=dumplines[i+9+j].split()
        ID = int(line0[3])
        idtype[ID-1][0]=ID 
        idtype[ID-1][1]=int(line0[4])
        atom[ID-1][0]=float(line0[0]) 
        atom[ID-1][1]=float(line0[1]) 
        atom[ID-1][2]=float(line0[2]) 
    return atom,idtype,  xlo,xhi_boundary,xy,ylo,yhi,xz,xhi
def huge_atoms(atom,xlo,xhi_boundary,xy,ylo,yhi,xz,xhi):
    a1=atom+[-xhi+xlo+xy,yhi-ylo,0]                    
    a2=atom+[xy,yhi-ylo,0]                      
    a3=atom+[xy+xhi-xlo,yhi-ylo,0]                    
    a4=atom+[-xhi+xlo,0,0]
    a5=atom+[0,0,0]
    a6=atom+[xhi-xlo,0,0]
    a7=atom+[-xy-xhi+xlo,-yhi+ylo,0]
    a8=atom+[-xy,-yhi+ylo,0] 
    a9=atom+[xhi-xlo-xy,-yhi+ylo,0]
    ATOMS=np.concatenate((a1,a2),axis=0)
    ATOMS=np.concatenate((ATOMS,a3),axis=0)
    ATOMS=np.concatenate((ATOMS,a4),axis=0)
    ATOMS=np.concatenate((ATOMS,a5),axis=0)
    ATOMS=np.concatenate((ATOMS,a6),axis=0)
    ATOMS=np.concatenate((ATOMS,a7),axis=0)
    ATOMS=np.concatenate((ATOMS,a8),axis=0)
    ATOMS=np.concatenate((ATOMS,a9),axis=0)
    return ATOMS 
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
    if (idtype[j,0]-1)%2544 >= 1080:
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

result_z=np.zeros((1000,4))

for p in range(1,n_file+1):
    with open("dump.trj","r") as f1:             
        dumplines = f1.readlines()
    natoms=int(dumplines[3].split()[0])
    zhen=0
    for i in range(0,len(dumplines),natoms+9) :
        time=int(dumplines[i+1].split()[0]) 
        idtype=np.zeros((natoms,2))
        atom=np.zeros((natoms,3))
        time_step=int(dumplines[i+1].split()[0])
        atom,idtype,  xlo,xhi_boundary,xy,ylo,yhi,xz,xhi=data_init(natoms,dumplines,i)                                
        big_atom=huge_atoms(atom,xlo,xhi_boundary,xy,ylo,yhi,xz,xhi)
        big_idtype=huge_idtype(idtype)
        L_boundary,R_boundary,U_boundary,D_boundary,U_Z_boundary,D_Z_boundary=cut_atom_area(distance,xlo,xhi_boundary,yhi,ylo,surfacelo,surfacehi)########
        U_Z_boundary,D_Z_boundary=350,-10
        big_atom,big_idtype ,sizeee=judge_cut_area(big_atom,big_idtype,L_boundary,R_boundary,U_boundary,D_boundary,U_Z_boundary,D_Z_boundary)######
        S=(xhi-xlo)*(yhi-ylo)
        count_zhen=np.zeros((1000,4))
        for j in range (natoms):                                                       
            if if_type(j,idtype,3):                                                
                x,y,z=position(atom,j)                                                 
                if if_soild(j,idtype):        
                    continue
                else:                               
                    count_zhen[int(z)+10,0]+=1                              
                                         
                    x_L,x_R,y_L,y_R,z_L,z_R=boundry(x,y,z,Ca_O)                        
                    x_L_H,x_R_H,y_L_H,y_R_H,z_L_H,z_R_H=boundry(x,y,z,H_O)          
                    flag_Ca=0                                                                                                           
                    for k in range(sizeee):                                   
                        if if_type(k,big_idtype,1):                                     
                            x_n,y_n,z_n=position(big_atom,k)                                
                            if judge_boumdry(x_n,y_n,z_n,x_L,x_R,y_L,y_R,z_L,z_R):          
                                if distence(x_n,y_n,z_n,x,y,z)<=Ca_O :    
                                    count_Ca_H=0
                                    flag_Ca=1
                                    for u in range(sizeee): 
                                        if if_type(u,big_idtype,4):
                                            x_H,y_H,z_H=position(big_atom,u)
                                            if judge_boumdry(x_H,y_H,z_H,x_L_H,x_R_H,y_L_H,y_R_H,z_L_H,z_R_H):                                    
                                                if distence(x_H,y_H,z_H,x,y,z)<=H_O :
                                                    count_Ca_H+=1
                                                else:continue
                                            else:continue
                                        else:continue
                                    if count_Ca_H==1:                              
                                        count_zhen[int(z)+10,2]+=1
                                        break
                                    elif count_Ca_H==2:
                                        count_zhen[int(z)+10,3]+=1
                                        break
                                    else:continue
                                else:continue
                            else:continue
                        else:continue
                    if flag_Ca==0:            
                        count_H2O=0
                        for zz in range(sizeee):  #H
                            if if_type(zz,big_idtype,4):
                                x_H2,y_H2,z_H2=position(big_atom,zz)
                                if judge_boumdry(x_H2,y_H2,z_H2,x_L_H,x_R_H,y_L_H,y_R_H,z_L_H,z_R_H):                                    
                                    if distence(x_H2,y_H2,z_H2,x,y,z)<=H_O :
                                        count_H2O+=1
                                    else:continue
                                else:continue
                            else:continue
                        if count_H2O==2:
                            count_zhen[int(z)+10,1]+=1                             
                        else:continue
            else:continue        
        result_z+=count_zhen/S
df=pd.DataFrame(result_z/10)
df.to_csv("all.csv")