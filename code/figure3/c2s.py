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
        ID = int(line0[0])
        idtype[ID-1][0]=ID 
        idtype[ID-1][1]=int(line0[1]) 
        atom[ID-1][0]=float(line0[3]) 
        atom[ID-1][1]=float(line0[4]) 
        atom[ID-1][2]=float(line0[5]) 
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
surfacelo=100
surfacehi=120
distance=10
n_file=210
Ca_7=[]
Ca_8=[]
for p in range(210,n_file+1):
    Ca=np.zeros((10000,80)) 
    with open("{}/dump2.trj".format(p),"r") as f1:      
        dumplines = f1.readlines()
    natoms=int(dumplines[3].split()[0])
    zhen=0
    for i in range(0,len(dumplines),natoms+9) :
        pin=2 
        time=int(dumplines[i+1].split()[0]) 
        idtype=np.zeros((natoms,2))
        atom=np.zeros((natoms,3))
        time_step=int(dumplines[i+1].split()[0])
        Ca[zhen,0]=time_step
        atom,idtype,  xlo,xhi_boundary,ylo,yhi,xhi=data_init(natoms,dumplines,i)                                
        big_atom=huge_atoms(atom,xlo,xhi_boundary,ylo,yhi,xhi)
        big_idtype=huge_idtype(idtype)
        L_boundary,R_boundary,U_boundary,D_boundary,U_Z_boundary,D_Z_boundary=cut_atom_area(distance,xlo,xhi_boundary,yhi,ylo,surfacelo,surfacehi)
        big_atom,big_idtype ,sizeee=judge_cut_area(big_atom,big_idtype,L_boundary,R_boundary,U_boundary,D_boundary,U_Z_boundary,D_Z_boundary)     
        for j in range (natoms):                                                  
            flag=1
            peiwei_sum=0
            peiwei_Os=0
            peiwei_Ow=0
            if if_soild(j,idtype):                                                       
                if if_surface(j,atom,surfacelo,surfacehi):                                  
                    if if_type(j,idtype,1):                                                  
                        x_Ca,y_Ca,z_Ca=position(atom,j)                                         
                        x_L,x_R,y_L,y_R,z_L,z_R=boundry(x_Ca,y_Ca,z_Ca,Ca_O)                        
                        Os=[]
                        Ow=[]                                                                       
                        for k in  range(sizeee):                                                
                            xx,yy,zz=position(big_atom,k)
                            if judge_boumdry(xx,yy,zz,x_L,x_R,y_L,y_R,z_L,z_R):                            
                                if if_type(k,big_idtype,3):                                        
                                    if if_soild(k,big_idtype):                                                                    
                                        if distence(xx,yy,zz,x_Ca,y_Ca,z_Ca)<=Ca_O:                      
                                            Os.append(k)
                                            peiwei_Os+=1
                                        else:continue
                                    else:
                                        if distence(xx,yy,zz,x_Ca,y_Ca,z_Ca)<=Ca_O:                   
                                            peiwei_Ow+=1
                                            Ow.append(k)                                            
                                        else:continue
                                else: continue
                            else:continue
                        peiwei_sum=peiwei_Os+peiwei_Ow
                        if peiwei_Os==5:
                            Ca[zhen,1]+=1
                            if peiwei_sum==5:
                                Ca[zhen,2]+=1
                            elif peiwei_sum==6:
                                Ca[zhen,3]+=1
                            elif peiwei_sum==7:
                                Ca[zhen,4]+=1
                                Ca_7.append(j)
                                O_H1=[]
                                O_H2=[]
                                for t in range(len(Ow)):                         
                                    countt=0
                                    x_Ow,y_Ow,z_Ow=position(big_atom,Ow[t])      
                                    x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow=boundry(x_Ow,y_Ow,z_Ow,H_O)
                                    for z in range(sizeee):                     
                                        xx_H,yy_H,zz_H=position(big_atom,z)       
                                        if judge_boumdry(xx_H,yy_H,zz_H,x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow):     
                                            if if_type(z,big_idtype,4):                          
                                                if distence(x_Ow,y_Ow,z_Ow,xx_H,yy_H,zz_H)<=H_O:     
                                                    countt+=1                                              
                                                else:continue
                                            else :continue
                                        else:continue
                                    if countt==1:
                                        O_H1.append(big_idtype[t,0])
                                    elif countt==2:
                                        O_H2.append(big_idtype[t,0])
                                    else: continue
                                if len(O_H1)==2 and len(O_H2)==0:
                                    Ca[zhen,5]+=1
                                if len(O_H1)==1 and len(O_H2)==1:
                                    Ca[zhen,6]+=1  
                                if len(O_H1)==0 and len(O_H2)==2:
                                    Ca[zhen,7]+=1  
                            elif peiwei_sum==8:
                                Ca[zhen,8]+=1
                                Ca_8.append(j)
                                O_H1=[]
                                O_H2=[]
                                for t in range(len(Ow)):                         
                                    countt=0
                                    x_Ow,y_Ow,z_Ow=position(big_atom,Ow[t])      
                                    x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow=boundry(x_Ow,y_Ow,z_Ow,H_O)
                                    for z in range(sizeee):                     
                                        xx_H,yy_H,zz_H=position(big_atom,z)       
                                        if judge_boumdry(xx_H,yy_H,zz_H,x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow):     
                                            if if_type(z,big_idtype,4):                          
                                                if distence(x_Ow,y_Ow,z_Ow,xx_H,yy_H,zz_H)<=H_O:     
                                                    countt+=1                                              
                                                else:continue
                                            else :continue
                                        else:continue
                                    if countt==1:
                                        O_H1.append(big_idtype[t,0])
                                    elif countt==2:
                                        O_H2.append(big_idtype[t,0])
                                    else: continue
                                if len(O_H1)==3 and len(O_H2)==0:
                                    Ca[zhen,9]+=1
                                if len(O_H1)==2 and len(O_H2)==1:
                                    Ca[zhen,10]+=1  
                                if len(O_H1)==1 and len(O_H2)==2:
                                    Ca[zhen,11]+=1  
                                if len(O_H1)==0 and len(O_H2)==3:
                                    Ca[zhen,12]+=1  
                            O_H1=[]
                            O_H2=[]
                            for t in range(len(Ow)):                         
                                countt=0
                                x_Ow,y_Ow,z_Ow=position(big_atom,Ow[t])      
                                x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow=boundry(x_Ow,y_Ow,z_Ow,H_O)
                                for z in range(sizeee):                     
                                    xx_H,yy_H,zz_H=position(big_atom,z)       
                                    if judge_boumdry(xx_H,yy_H,zz_H,x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow):     
                                        if if_type(z,big_idtype,4):                          
                                            if distence(x_Ow,y_Ow,z_Ow,xx_H,yy_H,zz_H)<=H_O:     
                                                countt+=1                                              
                                            else:continue
                                        else :continue
                                    else:continue
                                if countt==1:
                                    O_H1.append(big_idtype[t,0])
                                elif countt==2:
                                    O_H2.append(big_idtype[t,0])
                                else: continue
                            if len(O_H1)==1 and len(O_H2)==0:
                                Ca[zhen,13]+=1
                            if len(O_H2)==1 and len(O_H1)==0:
                                Ca[zhen,14]+=1  
                        elif peiwei_Os==4:
                            Ca[zhen,15]+=1
                            if peiwei_sum==4:
                                Ca[zhen,16]+=1
                            elif peiwei_sum==5:
                                Ca[zhen,17]+=1
                            elif peiwei_sum==6:
                                Ca[zhen,18]+=1
                            elif peiwei_sum==7:
                                Ca[zhen,19]+=1
                            elif peiwei_sum==8:
                                Ca[zhen,20]+=1   
                            O_H1=[]
                            O_H2=[]
                            for t in range(len(Ow)):                          
                                countt=0
                                x_Ow,y_Ow,z_Ow=position(big_atom,Ow[t])       
                                x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow=boundry(x_Ow,y_Ow,z_Ow,H_O)
                                for z in range(sizeee):                    
                                    xx_H,yy_H,zz_H=position(big_atom,z)       
                                    if judge_boumdry(xx_H,yy_H,zz_H,x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow):     
                                        if if_type(z,big_idtype,4):                       
                                            if distence(x_Ow,y_Ow,z_Ow,xx_H,yy_H,zz_H)<=H_O:       
                                                countt+=1                                       
                                            else:continue
                                        else :continue
                                    else:continue
                                if countt==1:
                                    O_H1.append(big_idtype[t,0])
                                elif countt==2:
                                    O_H2.append(big_idtype[t,0])
                                else: continue
                            if len(O_H1)==1 and len(O_H2)==1:
                                Ca[zhen,21]+=1         
                            if len(O_H2)==2 and len(O_H1)==0:
                                Ca[zhen,22]+=1
                            if len(O_H2)==0 and len(O_H1)==2:
                                Ca[zhen,23]+=1
                        elif peiwei_Os==3:
                            Ca[zhen,24]+=1
                            if peiwei_sum==3:
                                Ca[zhen,25]+=1
                            elif peiwei_sum==4:
                                Ca[zhen,26]+=1
                            elif peiwei_sum==5:
                                Ca[zhen,27]+=1
                            elif peiwei_sum==6:
                                Ca[zhen,28]+=1
                            elif peiwei_sum==7:
                                Ca[zhen,29]+=1   
                            elif peiwei_sum==8:
                                Ca[zhen,30]+=1   
                            O_H1=[]
                            O_H2=[]
                            for t in range(len(Ow)):                          
                                countt=0
                                x_Ow,y_Ow,z_Ow=position(big_atom,Ow[t])       
                                x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow=boundry(x_Ow,y_Ow,z_Ow,H_O)
                                for z in range(sizeee):                        
                                    xx_H,yy_H,zz_H=position(big_atom,z)        
                                    if judge_boumdry(xx_H,yy_H,zz_H,x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow):    
                                        if if_type(z,big_idtype,4):                            
                                            if distence(x_Ow,y_Ow,z_Ow,xx_H,yy_H,zz_H)<=H_O:      
                                                countt+=1                                             
                                            else:continue
                                        else :continue
                                    else:continue
                                if countt==1:
                                    O_H1.append(big_idtype[t,0])
                                elif countt==2:
                                    O_H2.append(big_idtype[t,0])
                                else: continue
                            if len(O_H1)==1 and len(O_H2)==2:
                                Ca[zhen,31]+=1
                            if len(O_H1)==0 and len(O_H2)==3:
                                Ca[zhen,32]+=1                                
                            if len(O_H2)==1 and len(O_H1)==2:
                                Ca[zhen,33]+=1
                            if len(O_H2)==0 and len(O_H1)==3: 
                                Ca[zhen,34]+=1
                                
                        elif peiwei_Os==2:
                            Ca[zhen,35]+=1
                            if peiwei_sum==2:
                                Ca[zhen,36]+=1
                            elif peiwei_sum==3:
                                Ca[zhen,37]+=1
                            elif peiwei_sum==4:
                                Ca[zhen,38]+=1
                            elif peiwei_sum==5:
                                Ca[zhen,39]+=1
                            elif peiwei_sum==6:
                                Ca[zhen,40]+=1   
                            elif peiwei_sum==7:
                                Ca[zhen,41]+=1  
                            elif peiwei_sum==8:
                                Ca[zhen,42]+=1  
                            O_H1=[]
                            O_H2=[]
                            for t in range(len(Ow)):                          
                                countt=0
                                x_Ow,y_Ow,z_Ow=position(big_atom,Ow[t])        
                                x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow=boundry(x_Ow,y_Ow,z_Ow,H_O)
                                for z in range(sizeee):                        
                                    xx_H,yy_H,zz_H=position(big_atom,z)       
                                    if judge_boumdry(xx_H,yy_H,zz_H,x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow):    
                                        if if_type(z,big_idtype,4):                            
                                            if distence(x_Ow,y_Ow,z_Ow,xx_H,yy_H,zz_H)<=H_O:      
                                                countt+=1                                            
                                            else:continue
                                        else :continue
                                    else:continue
                                if countt==1:
                                    O_H1.append(big_idtype[t,0])
                                elif countt==2:
                                    O_H2.append(big_idtype[t,0])
                                else: continue
                            if len(O_H1)==1 and len(O_H2)==3:
                                Ca[zhen,43]+=1
                            if len(O_H1)==0 and len(O_H2)==4:
                                Ca[zhen,44]+=1
                            if len(O_H1)==2 and len(O_H2)==2:
                                Ca[zhen,45]+=1
                            if len(O_H1)==3 and len(O_H2)==1:
                                Ca[zhen,46]+=1
                            if len(O_H1)==4 and len(O_H2)==0:
                                Ca[zhen,47]+=1
                        elif peiwei_Os==1:
                            Ca[zhen,48]+=1
                            if peiwei_sum==1:
                                Ca[zhen,49]+=1
                            elif peiwei_sum==2:
                                Ca[zhen,50]+=1
                            elif peiwei_sum==3:
                                Ca[zhen,51]+=1
                            elif peiwei_sum==4:
                                Ca[zhen,52]+=1
                            elif peiwei_sum==5:
                                Ca[zhen,53]+=1   
                            elif peiwei_sum==6:
                                Ca[zhen,54]+=1  
                            elif peiwei_sum==7:
                                Ca[zhen,55]+=1  
                            elif peiwei_sum==8:
                                Ca[zhen,56]+=1  
                            O_H1=[]
                            O_H2=[]
                            for t in range(len(Ow)):                          
                                countt=0
                                x_Ow,y_Ow,z_Ow=position(big_atom,Ow[t])        
                                x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow=boundry(x_Ow,y_Ow,z_Ow,H_O)
                                for z in range(sizeee):                      
                                    xx_H,yy_H,zz_H=position(big_atom,z)        
                                    if judge_boumdry(xx_H,yy_H,zz_H,x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow):   
                                        if if_type(z,big_idtype,4):                          
                                            if distence(x_Ow,y_Ow,z_Ow,xx_H,yy_H,zz_H)<=H_O:      
                                                countt+=1                                              
                                            else:continue
                                        else :continue
                                    else:continue
                                if countt==1:
                                    O_H1.append(big_idtype[t,0])
                                elif countt==2:
                                    O_H2.append(big_idtype[t,0])
                                else: continue
                            if len(O_H1)==0 and len(O_H2)==5:
                                Ca[zhen,57]+=1
                            if len(O_H1)==1 and len(O_H2)==4:
                                Ca[zhen,58]+=1       
                            if len(O_H1)==2 and len(O_H2)==3:
                                Ca[zhen,59]+=1
                            if len(O_H1)==3 and len(O_H2)==2:
                                Ca[zhen,60]+=1
                            if len(O_H1)==4 and len(O_H2)==1:
                                Ca[zhen,61]+=1
                            if len(O_H1)==5 and len(O_H2)==0:
                                Ca[zhen,62]+=1
                        if peiwei_Os<=4:
                            Ca[zhen,63]+=1
                            O_H1=[]
                            O_H2=[]
                            for t in range(len(Ow)):                        
                                countt=0
                                x_Ow,y_Ow,z_Ow=position(big_atom,Ow[t])        
                                x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow=boundry(x_Ow,y_Ow,z_Ow,H_O)
                                for z in range(sizeee):                         
                                    xx_H,yy_H,zz_H=position(big_atom,z)   
                                    if judge_boumdry(xx_H,yy_H,zz_H,x_L_Ow,x_R_Ow,y_L_Ow,y_R_Ow,z_L_Ow,z_R_Ow):   
                                        if if_type(z,big_idtype,4):                        
                                            if distence(x_Ow,y_Ow,z_Ow,xx_H,yy_H,zz_H)<=H_O:      
                                                countt+=1                                             
                                            else:continue
                                        else :continue
                                    else:continue
                                if countt==1:
                                    O_H1.append(big_idtype[t,0])
                                elif countt==2:
                                    O_H2.append(big_idtype[t,0])
                                else: continue
                            if len(O_H1)==1 and len(O_H2)==0:
                                Ca[zhen,64]+=1
                            if len(O_H1)==2 and len(O_H2)==0:
                                Ca[zhen,65]+=1
                            if len(O_H1)==3 and len(O_H2)==0:
                                Ca[zhen,66]+=1
                        else:continue
                    else:continue    
                else:continue        
            else:continue            
        zhen+=1                
    df=pd.DataFrame(Ca)
    df.to_csv("Ca_peiwei{}.csv".format(p)) 
    Ca_7=list(set(Ca_7))
    Ca_8=list(set(Ca_8))
    df=pd.DataFrame(Ca_7)
    df.to_csv("Ca_7_id{}.csv".format(p))
    df=pd.DataFrame(Ca_8)
    df.to_csv("Ca_8_id{}.csv".format(p))



