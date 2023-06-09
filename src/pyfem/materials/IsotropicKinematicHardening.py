

from pyfem.materials.BaseMaterial import BaseMaterial
from pyfem.materials.MatUtils     import vonMisesStress,hydrostaticStress
from pyfem.materials.MatUtils     import transform3To2,transform2To3
from numpy import zeros, ones, dot, array, outer
from math import sqrt

class IsotropicKinematicHardening( BaseMaterial ):

  def __init__ ( self, props ):

    self.tolerance = 1.0e-6

    BaseMaterial.__init__( self, props )

    self.ebulk3 = self.E / ( 1.0 - 2.0*self.nu )
    self.eg2    = self.E / ( 1.0 + self.nu )
    self.eg     = 0.5*self.eg2
    self.eg3    = 3.0*self.eg
    self.elam   = ( self.ebulk3 - self.eg2 ) / 3.0

    self.ctang = zeros(shape=(6,6))

    self.ctang[:3,:3] = self.elam

    self.ctang[0,0] += self.eg2
    self.ctang[1,1] = self.ctang[0,0]
    self.ctang[2,2] = self.ctang[0,0]

    self.ctang[3,3] = self.eg
    self.ctang[4,4] = self.ctang[3,3]
    self.ctang[5,5] = self.ctang[3,3]
 
    self.setHistoryParameter( 'sigma', zeros(6) )
    self.setHistoryParameter( 'eelas', zeros(6) )
    self.setHistoryParameter( 'eplas', zeros(6) )
    self.setHistoryParameter( 'alpha', zeros(6) )

    self.commit_history()

    #Set the labels for the output data in this material model
    self.outLabels = [ "S11" , "S22" , "S33" , "S23" , "S13" , "S12" , "Epl" ]
    self.outData   = zeros(7)

#------------------------------------------------------------------------------
#  pre:  kinematics object containing current strain (kinemtics.strain)
#  post: stress vector and tangent matrix
#------------------------------------------------------------------------------

  def getStress( self, kinematics ):

    eelas = self.getHistoryParameter('eelas')   
    eplas = self.getHistoryParameter('eplas')   
    alpha = self.getHistoryParameter('alpha')   
    sigma = self.getHistoryParameter('sigma') 
     
    if len(kinematics.dstrain) == 6:
      dstrain = kinematics.dstrain
    else:
      dstrain = transform2To3(kinematics.dstrain)

    eelas += dstrain

    sigma += dot( self.ctang , dstrain )

    tang = self.ctang

    smises = vonMisesStress( sigma - alpha )

    deqpl = 0.

    if smises > ( 1.0 + self.tolerance ) * self.syield:
      shydro = hydrostaticStress( sigma )
   
      flow = sigma - alpha

      flow[:3] = flow[:3]-shydro*ones(3)
      flow *= 1.0/smises

      deqpl = ( smises - self.syield ) / ( self.eg3 + self.hard )

      alpha += self.hard * flow * deqpl
      eplas[:3] +=  1.5 * flow[:3] * deqpl
      eelas[:3] += -1.5 * flow[:3] * deqpl

      eplas[3:] +=  3.0 * flow[3:] * deqpl
      eelas[3:] += -3.0 * flow[3:] * deqpl

      sigma = alpha + flow * self.syield
      sigma[:3] += shydro * ones(3)
     
      effg   = self.eg*(self.syield + self.hard*deqpl ) / smises
      effg2  = 2.0*effg
      effg3  = 3.0*effg
      efflam = 1.0/3.0 * ( self.ebulk3-effg2 )
      effhdr = self.eg3 * self.hard/(self.eg3+self.hard)-effg3

      tang = zeros(shape=(6,6))
      tang[:3,:3] = efflam
    
      for i in range(3):
        tang[i,i]     += effg2
        tang[i+3,i+3] += effg

      tang += effhdr*outer(flow,flow)
 
    self.setHistoryParameter( 'eelas', eelas )
    self.setHistoryParameter( 'eplas', eplas )
    self.setHistoryParameter( 'alpha', alpha )
    self.setHistoryParameter( 'sigma', sigma )

    # Store output eplas

    self.outData[:6] = sigma
    self.outData[6]  = eplas[0]

    if len(kinematics.dstrain) == 6:
      return sigma , tang  
    else:
      return transform3To2(sigma,tang)
         
 
