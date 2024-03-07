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
        ID = int(line0[0])
        idtype[ID-1][0]=ID 
        idtype[ID-1][1]=int(line0[1]) 
        atom[ID-1][0]=float(line0[3]) 
        atom[ID-1][1]=float(line0[4]) 
        atom[ID-1][2]=float(line0[5])
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
def if_soild(j):
    if j%2544 >= 1080 :
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
Ca_O=2.83
H_O=1.2
surfacelo=130 
surfacehi=150
n_file=161
for p in range(157,n_file+1):
    Ca=np.zeros((10000,10))   
    with open("{}/dump2.trj".format(p),"r") as f1:             
        dumplines = f1.readlines()
    natoms=int(dumplines[3].split()[0])
    zhen=0
    Ca_structure=np.zeros((10000,13))
    for i in range(0,len(dumplines),natoms+9) :
        time=int(dumplines[i+1].split()[0]) 
        idtype=np.zeros((natoms,2))
        atom=np.zeros((natoms,3))
        time_step=int(dumplines[i+1].split()[0])
        atom,idtype,  xlo,xhi_boundary,xy,ylo,yhi,xz,xhi=data_init(natoms,dumplines,i)                                
        big_atom=huge_atoms(atom,xlo,xhi_boundary,xy,ylo,yhi,xz,xhi)
        big_idtype=huge_idtype(idtype)
        for j in range (natoms):                                                  
            if if_soild(j):                                                          
                if if_surface(j,atom,surfacelo,surfacehi):                                
                    if if_type(j,idtype,1):                                                 
                        x_Ca,y_Ca,z_Ca=position(atom,j)                                     
                        flag=1
                        x_L,x_R,y_L,y_R,z_L,z_R=boundry(x_Ca,y_Ca,z_Ca,Ca_O)                      
                        id_Ca_O=np.zeros(500)                                                     
                        id_Ca_O[0]=j                                                            
                        num_Ca_O=1
                        Ca_O_structure=0                                                                
                        for k in  range(natoms*9):                                              
                            xx,yy,zz=position(big_atom,k)
                            if judge_boumdry(xx,yy,zz,x_L,x_R,y_L,y_R,z_L,z_R):                           
                                if if_type(k,big_idtype,3):                                       
                                    if if_soild(k):                                                                      
                                        x_O,y_O,z_O=position(big_atom,k)
                                        if distence(x_O,y_O,z_O,x_Ca,y_Ca,z_Ca)<=Ca_O:                    
                                            flag=0
                                            break
                                    else:
                                        x_O_w,y_O_w,z_O_w=position(big_atom,k)              
                                        if distence(x_O_w,y_O_w,z_O_w,x_Ca,y_Ca,z_Ca)<=Ca_O:
                                            id_Ca_O[num_Ca_O]=k
                                            num_Ca_O+=1
                                            Ca_O_structure+=1
                                        continue
                                else: continue
                            else:continue
                        if flag==1:
                            Ca[zhen,1]+=1                                       
                            if Ca_O_structure==5:
                                Ca[zhen,2]+=1
                            elif Ca_O_structure==6:                             
                                Ca[zhen,3]+=1
                            elif Ca_O_structure==7:
                                Ca[zhen,4]+=1
                            elif Ca_O_structure==4:
                                Ca[zhen,5]+=1
                            elif Ca_O_structure==3:
                                Ca[zhen,6]+=1
                            elif Ca_O_structure==2:
                                Ca[zhen,7]+=1
                            elif Ca_O_structure==1:
                                Ca[zhen,8]+=1
                            count_structure_Ca_O_H=0
                            count_structure_Ca_O_H2=0
                            for r in range(1,num_Ca_O+1):                                                              
                                x_O_ww,y_O_ww,z_O_ww=position(big_atom,int(id_Ca_O[r]))
                                x_L_ww,x_R_ww,y_L_ww,y_R_ww,z_L_ww,z_R_ww=boundry(x_O_ww,y_O_ww,z_O_ww,H_O)
                                count_OH=0
                                for q in range(natoms*9):                                                             
                                    xx_H,yy_H,zz_H=position(big_atom,q)
                                    if judge_boumdry(xx_H,yy_H,zz_H,    x_L_ww,x_R_ww,y_L_ww,y_R_ww,z_L_ww,z_R_ww):
                                        if if_type(q,big_idtype,4):
                                            if distence(xx_H,yy_H,zz_H,x_O_ww,y_O_ww,z_O_ww)<=H_O:
                                                count_OH+=1
                                            else:
                                                continue   
                                        else:continue
                                    else:continue     
                                if count_OH==1:
                                    count_structure_Ca_O_H+=1                 
                                if count_OH==2:
                                    count_structure_Ca_O_H2+=1                                
                            if count_structure_Ca_O_H==3:
                                Ca_structure[zhen,1]+=1
                            elif count_structure_Ca_O_H==2:
                                Ca_structure[zhen,2]+=1
                            elif count_structure_Ca_O_H==1:
                                Ca_structure[zhen,3]+=1
                            elif count_structure_Ca_O_H==0:
                                Ca_structure[zhen,4]+=1
                            if count_structure_Ca_O_H2==3:
                                Ca_structure[zhen,5]+=1
                            elif count_structure_Ca_O_H2==4:
                                Ca_structure[zhen,6]+=1
                            elif count_structure_Ca_O_H2==5:
                                Ca_structure[zhen,7]+=1
                            elif count_structure_Ca_O_H2==6:
                                Ca_structure[zhen,8]+=1       
                            if count_structure_Ca_O_H==3 and count_structure_Ca_O_H2==3:
                                Ca_structure[zhen,9]+=1
                            if count_structure_Ca_O_H==2 and count_structure_Ca_O_H2==4:
                                Ca_structure[zhen,10]+=1
                            if count_structure_Ca_O_H==1 and count_structure_Ca_O_H2==5:
                                Ca_structure[zhen,11]+=1
                            if count_structure_Ca_O_H==0 and count_structure_Ca_O_H2==6:
                                Ca_structure[zhen,12]+=1
        Ca_structure[zhen,0]= time_step                       
        Ca[zhen,0]=time_step
        zhen+=1
    df=pd.DataFrame(Ca)
    df.to_csv("Ca_{}.csv".format(p))     
    df=pd.DataFrame(Ca_structure)
    df.to_csv("Ca_structure_{}.csv".format(p))       
