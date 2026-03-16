import numpy as np
from sympy.physics.quantum.cg import CG

#%%%Class for constructing operators, stores the operator in different bases to allow for easy switching.
class Operator():
    
    def __init__(self,I,J,U_IJtoFmF=''):      
        totalstates = int((2*I+1)*(2*J+1))
        Fmin        = np.absolute(I-J)
        Fmax        = np.absolute(I+J)

        # Make some lists which indicate what state corresponds to what index in each base.
        mImJindexarray = np.zeros((totalstates,2),dtype=np.float64)
        index = 0
        for mI in np.arange(-I,I+1,1):
            for mJ in np.arange(-J,J+1,1):
                mImJindexarray[index] = [mI,mJ]
                index += 1
                
        FmFindexarray  = np.zeros((totalstates,2),dtype=np.float64)
        index = 0
        for f in np.arange(Fmin,Fmax+1,1):
            for mf in np.arange(-f,f+1,1):
                FmFindexarray[index] = [f,mf]
                index += 1
                
        self.mImJ = np.zeros((totalstates,totalstates),dtype=np.complex128)
        self.FmF  = np.zeros((totalstates,totalstates),dtype=np.complex128)
        self.mImJindex = mImJindexarray
        self.FmFindex = FmFindexarray

        # Basis transformations between mI,mJ and F,mF
        if type(U_IJtoFmF) == str:
            self.mI_mJ_to_F_mF = np.zeros((totalstates,totalstates),dtype=np.complex128)
            self.F_mF_to_mI_mJ = np.zeros((totalstates,totalstates),dtype=np.complex128)
            for x in np.arange(totalstates):
                for y in np.arange(totalstates):
                    tmp1 = mImJindexarray[x]
                    tmp2 = FmFindexarray[y]

                    mi,mj   = tmp1[0], tmp1[1]
                    f,mf    = tmp2[0], tmp2[1]

                    self.mI_mJ_to_F_mF[y,x] = CG(I,mi,J,mj,f,mf).doit()
                    self.F_mF_to_mI_mJ[x,y] = CG(I,mi,J,mj,f,mf).doit()
        else:
            self.mI_mJ_to_F_mF = U_IJtoFmF
            self.F_mF_to_mI_mJ = np.transpose(U_IJtoFmF)
		
			 
        
    # Functions for setting the operator 
    def set_FmF(self,oprt):
        self.FmF = oprt
        self.mImJ = np.dot(self.F_mF_to_mI_mJ,np.dot(oprt,self.mI_mJ_to_F_mF))
    
    def set_mImJ(self,oprt):
        self.mImJ = oprt
        self.FmF = np.dot(self.mI_mJ_to_F_mF,np.dot(oprt,self.F_mF_to_mI_mJ))
        
        
    
#%%%    
#Class that holds the angular momentum operators, all given in their natural basis.
class AngMomOperators():
    
    def __init__(self,I,J):
        
        
        totalstates = int((2*I+1)*(2*J+1))    
        self.totalstates = int((2*I+1)*(2*J+1))    

        dum  = Operator(I,J)
        basistransformation = dum.mI_mJ_to_F_mF
        #%%%% Create the objects for the angular momentum operators
        self.Isqrd = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Iplus = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Iminus = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Ix = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Iy = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Iz = Operator(I,J,U_IJtoFmF=basistransformation)

        self.Jsqrd = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Jplus = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Jminus = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Jx = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Jy = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Jz = Operator(I,J,U_IJtoFmF=basistransformation)

        self.IdotJ = Operator(I,J,U_IJtoFmF=basistransformation)

        self.Fsqrd = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Fplus = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Fminus = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Fx = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Fy = Operator(I,J,U_IJtoFmF=basistransformation)
        self.Fz = Operator(I,J,U_IJtoFmF=basistransformation)
        
        self.Id = Operator(I,J,U_IJtoFmF=basistransformation)

        
        #%%%% Temporary matrices for constructing the operators.
        matIsqrd       = np.zeros((totalstates,totalstates),dtype=np.complex128)
        matIz          = np.zeros((totalstates,totalstates),dtype=np.complex128)
        matIplus       = np.zeros((totalstates,totalstates),dtype=np.complex128)
        matIminus      = np.zeros((totalstates,totalstates),dtype=np.complex128)

        matJsqrd       = np.zeros((totalstates,totalstates),dtype=np.complex128)
        matJz          = np.zeros((totalstates,totalstates),dtype=np.complex128)
        matJplus       = np.zeros((totalstates,totalstates),dtype=np.complex128)
        matJminus      = np.zeros((totalstates,totalstates),dtype=np.complex128)

        matFsqrd       = np.zeros((totalstates,totalstates),dtype=np.complex128)
        matFz          = np.zeros((totalstates,totalstates),dtype=np.complex128)
        matFplus       = np.zeros((totalstates,totalstates),dtype=np.complex128)
        matFminus      = np.zeros((totalstates,totalstates),dtype=np.complex128)


        #Construct the diagonal operators
        matIsqrd = np.identity(totalstates)*I*(I+1)
        matJsqrd = np.identity(totalstates)*J*(J+1.0)
        matFsqrd = np.diag(dum.FmFindex[:,0]*(dum.FmFindex[:,0]+1.0))
        matIz    = np.diag(dum.mImJindex[:,0])
        matJz    = np.diag(dum.mImJindex[:,1])
        matFz    = np.diag(dum.FmFindex[:,1])

        
        # constructing the ladder operators
        Fmin        = np.absolute(I-J)
        Fmax        = np.absolute(I+J)
        skipindex = 0   
        for f in np.arange(Fmin,Fmax+1,1):
            mf = np.arange(-f,f,1)
            matFplus[skipindex:skipindex+int(2*f+1),skipindex:skipindex+int(2*f+1)] = np.diag(np.sqrt(f*(f+1) - (mf)*(mf+1)),k=-1)
            matFminus[skipindex:skipindex+int(2*f+1),skipindex:skipindex+int(2*f+1)]= np.diag(np.sqrt(f*(f+1) - mf*(mf+1)),k=+1)
            skipindex += int(2*f+1)
            
        skipindex = 0
        mJ = np.arange(-J,J,1)
        for mI in np.arange(-I,I+1,1):
            matJplus[skipindex:skipindex+int(2*J+1),skipindex:skipindex+int(2*J+1)] = np.diag(np.sqrt(J*(J+1) - mJ*(mJ+1)),k=-1)
            matJminus[skipindex:skipindex+int(2*J+1),skipindex:skipindex+int(2*J+1)]= np.diag(np.sqrt(J*(J+1) - mJ*(mJ+1)) ,k=+1)
            skipindex += int(2 * J +1)
            
        mI = np.arange(-I,I,1)
        for mJ in np.arange(-J,J+1,1):
            index = int(mJ+J)
            matIplus[index::int(2*J+1),index::int(2*J+1)] = np.diag(np.sqrt(I*(I+1) - mI*(mI+1)),k=-1)
            matIminus[index::int(2*J+1),index::int(2*J+1)] = np.diag(np.sqrt(I*(I+1) - mI*(mI+1)),k=+1)
        
        self.Isqrd.set_mImJ(matIsqrd)
        self.Iplus.set_mImJ(matIplus)
        self.Iminus.set_mImJ(matIminus)
        self.Ix.set_mImJ(0.5*(matIplus + matIminus))
        self.Iy.set_mImJ(0.5j*(matIminus - matIplus))
        self.Iz.set_mImJ(matIz)

        self.Jsqrd.set_mImJ(matJsqrd)
        self.Jplus.set_mImJ(matJplus)
        self.Jminus.set_mImJ(matJminus)
        self.Jx.set_mImJ(0.5*(matJplus + matJminus))
        self.Jy.set_mImJ(0.5j*(matJminus - matJplus))
        self.Jz.set_mImJ(matJz)

        self.Fsqrd.set_FmF(matFsqrd)
        self.Fplus.set_FmF(matFplus)
        self.Fminus.set_FmF(matFminus)
        self.Fx.set_FmF(0.5*(matFplus + matFminus))
        self.Fy.set_FmF(0.5j*(matFminus - matFplus) )
        self.Fz.set_FmF(matFz)

        self.IdotJ.set_FmF(0.5*self.Fsqrd.FmF)
        self.IdotJ.set_mImJ(self.IdotJ.mImJ - 0.5*self.Jsqrd.mImJ - 0.5*self.Isqrd.mImJ)
        
        
        self.Id.set_mImJ(np.identity(totalstates))
        
        return 
