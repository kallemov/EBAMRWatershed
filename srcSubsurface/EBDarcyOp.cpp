#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBArith.H"
#include "EBDarcyOp.H"
#include "EBDarcyOpF_F.H"
#include "NamespaceHeader.H"

bool EBDarcyOp::s_turnOffBCs = false; //REALLY needs to default to false
//bool EBDarcyOp::s_forceNoEBCF = false; //REALLY needs to default to false

//-----------------------------------------------------------------------
EBDarcyOp::
EBDarcyOp(const EBLevelGrid &                                  a_eblgFine,
	  const EBLevelGrid &                                  a_eblg,
	  const EBLevelGrid &                                  a_eblgCoar,
	  const RefCountedPtr<RichardsDomainBC>&               a_domainBC,
	  const RefCountedPtr<RichardsEBBC>&                   a_ebBC,
	  const RefCountedPtr<EBQuadCFInterp>&                 a_quadCFI,
	  const int&                                           a_refToFine,
	  const int&                                           a_refToCoar,
	  const RealVect&                                      a_vectDx,
	  const Real&                                          a_beta,
	  const RefCountedPtr<LevelData<EBFluxFAB> >&          a_bcoef,
	  const RefCountedPtr<LevelData<BaseIVFAB<Real> > >&   a_bcoIrreg,
	  const RefCountedPtr<LevelData<EBFluxFAB> >&          a_darcyFlux,
	  const IntVect&                                       a_ghostCellsPhi,
	  const IntVect&                                       a_ghostCellsRHS)
   :m_ghostCellsPhi(a_ghostCellsPhi),
    m_ghostCellsRHS(a_ghostCellsRHS),
    m_quadCFIWithCoar(a_quadCFI),
    m_eblg(a_eblg),
    m_domainBC(a_domainBC),
    m_ebBC(a_ebBC),
    m_vectDx(a_vectDx),
    m_bcoef(a_bcoef),
    m_bcoIrreg(a_bcoIrreg),
    m_DarcyFlux(a_darcyFlux),
    m_beta(a_beta),
    m_opEBStencil(),
    m_vofIterIrreg(),
    m_vofIterMulti(),
    m_vofIterDomLo(),
    m_vofIterDomHi(),
    m_fastFR(),
    m_refToFine(a_refToFine),
    m_refToCoar(a_refToCoar),
    m_hasFine(a_refToFine>0 ? true : false),
    m_hasCoar(a_refToCoar>0 ? true : false)
{
  CH_TIME("EBDarcyOp::DarcyOp");
  //  int ncomp = 1;

  if (m_hasFine)
    {
      m_eblgFine       = a_eblgFine;
      m_vectDxFine     = m_vectDx/a_refToFine;
    }

  if (m_hasCoar)
    {
      m_eblgCoar       = a_eblgCoar;
      m_vectDxCoar     = m_vectDx*a_refToCoar;
    }
  //define stencils for the operator
  // defineStencils();

}
//-----------------------------------------------------------------------
EBDarcyOp::
~EBDarcyOp()
{
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
fillGrad(const LevelData<EBCellFAB>& a_phi)
{
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
calculateAlphaWeight()
{
  DataIterator dit = m_eblg.getDBL().dataIterator(); 
  int nbox = dit.size();
//#pragma omp parallel
  {
//#pragma omp for
    for(int mybox=0; mybox<nbox; mybox++)
      {
	  m_alphaDiagWeight[dit[mybox]].setVal(0.0);
     }
  }
}
//-----------------------------------------------------------------------

void
EBDarcyOp::
getDivFStencil(VoFStencil&      a_vofStencil,
               const VolIndex&  a_vof,
               const DataIndex& a_dit)
{
  CH_TIME("EBDarcyOp::getDivFStencil");
  const EBISBox& ebisBox = m_eblg.getEBISL()[a_dit];
  a_vofStencil.clear();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int isign = sign(sit());
          Vector<FaceIndex> faces = ebisBox.getFaces(a_vof, idir, sit());

          for (int iface = 0; iface < faces.size(); iface++)
            {
              VoFStencil fluxStencil;
              getFluxStencil(fluxStencil, faces[iface], a_dit, idir);
              Real areaFrac = ebisBox.areaFrac(faces[iface]);
              fluxStencil *= Real(isign)*areaFrac/m_vectDx[idir];
              a_vofStencil += fluxStencil;

            }
        }
    }
}

void
EBDarcyOp::
getJacobianStencil(VoFStencil&      a_vofStencil,
		   const VolIndex&  a_vof,
		   const DataIndex& a_dit)
{
  //first let's get what is the terms from the regular stencil
  getDivFStencil(a_vofStencil,a_vof, a_dit);

  LayoutData<BaseIVFAB<VoFStencil> >* fluxStencil = m_ebBC->getFluxStencil(0);
  if (fluxStencil != NULL)
    {
      BaseIVFAB<VoFStencil>& fluxStencilBaseIVFAB = (*fluxStencil)[a_dit];
      //this fills the stencil with the gradient*beta*bcoef
      VoFStencil  fluxStencilPt = fluxStencilBaseIVFAB(a_vof,0);
      a_vofStencil += fluxStencilPt;
    }
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
getFluxStencil(VoFStencil&      a_fluxStencil,
               const FaceIndex& a_face,
               const DataIndex& a_dit,
	       const int a_dir)
{
  /// stencil for flux computation.   the truly ugly part of this computation
  /// beta and eta are multiplied in here

  CH_TIME("EBDarcyOp::getFluxStencil");
  //need to do this by interpolating to centroids
  //so get the stencil at each face center and add with
  //interpolation weights
  FaceStencil interpSten = EBArith::getInterpStencil(a_face,
                                                     (*(m_eblg.getCFIVS()))[a_dit],
                                                     m_eblg.getEBISL()[a_dit],
                                                     m_eblg.getDomain());

  a_fluxStencil.clear();

  for (int isten = 0; isten < interpSten.size(); isten++)
  {
      const FaceIndex& face = interpSten.face(isten);
      const Real&    weight = interpSten.weight(isten);
      VoFStencil faceCentSten;
      getFaceCenteredFluxStencil(faceCentSten, face, a_dit, a_dir);
      faceCentSten *= weight;
      a_fluxStencil += faceCentSten;
  
  }
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
                           const FaceIndex& a_face,
                           const DataIndex& a_dit,
			   const int a_dir)
{
  CH_TIME("EBDarcyOp::getFaceCenteredFluxStencil");
  //face centered gradient is just a centered diff
  //int faceDir= a_face.direction();
  a_fluxStencil.clear();

  if (!a_face.isBoundary())
    {
      a_fluxStencil.add(a_face.getVoF(Side::Hi),  1.0/m_vectDx[a_dir], 0);
      a_fluxStencil.add(a_face.getVoF(Side::Lo), -1.0/m_vectDx[a_dir], 0);
      a_fluxStencil *= (*m_bcoef)[a_dit][a_dir](a_face,0);
    }
  else
    {
      //the boundary condition handles this one.
    }
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
setBeta(const Real& a_beta)
{
  CH_TIME("EBDarcyOp::setAlphaAndBeta");
  m_beta  = a_beta;
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
setTime(Real a_time)
{
    m_time = a_time; //used in BC
    // Refine the stencils.
    //defineStencils();

}
//-----------------------------------------------------------------------
void
EBDarcyOp::
kappaScale(LevelData<EBCellFAB> & a_rhs)
{
  CH_TIME("EBDarcyOp::kappaScale");
  EBLevelDataOps::kappaWeight(a_rhs);
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
diagonalScale(LevelData<EBCellFAB> & a_rhs,
              bool a_kappaWeighted)
{

  CH_TIME("EBDarcyOp::diagonalScale");
  if (a_kappaWeighted)
    EBLevelDataOps::kappaWeight(a_rhs);
}
void
EBDarcyOp::
defineStencils()
{
  CH_TIME("EBDarcyOp::defineStencils");
  // create ebstencil for irregular applyOp
  m_opEBStencil.define(m_eblg.getDBL());
  // create vofstencils for applyOp and

  Real fakeBeta = 1.0;
  m_domainBC->setCoef(m_eblg,   fakeBeta ,      m_bcoef   );
  m_ebBC->setCoef(    m_eblg,   fakeBeta ,      m_bcoIrreg);

  //Real dxScale = 1.0/m_dx;
  //m_ebBC->define((*(m_eblg.getCFIVS())), dxScale); //has to happen AFTER coefs are set
  //you should assert this is NULL
  LayoutData<BaseIVFAB<VoFStencil> >* fluxStencil = m_ebBC->getFluxStencil(0);

  m_vofIterIrreg.define(     m_eblg.getDBL()); // vofiterator cache
  m_vofIterMulti.define(     m_eblg.getDBL()); // vofiterator cache
  m_alphaDiagWeight.define(  m_eblg.getDBL());
  Box sideBoxLo[SpaceDim];
  Box sideBoxHi[SpaceDim];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box domainBox = m_eblg.getDomain().domainBox();
      sideBoxLo[idir] = adjCellLo(domainBox, idir, 1);
      sideBoxLo[idir].shift(idir,  1);
      sideBoxHi[idir] = adjCellHi(domainBox, idir, 1);
      sideBoxHi[idir].shift(idir, -1);
      m_vofIterDomLo[idir].define( m_eblg.getDBL()); // vofiterator cache for domain lo
      m_vofIterDomHi[idir].define( m_eblg.getDBL()); // vofiterator cache for domain hi
    }

  DataIterator dit = m_eblg.getDBL().dataIterator(); 
  int nbox = dit.size();
  //#pragma omp parallel
  {
    // this breaks 
    //#pragma omp parallel for
    for(int mybox=0; mybox<nbox; mybox++)
      {

        const Box& curBox = m_eblg.getDBL().get(dit[mybox]);
        const EBISBox& ebisBox = m_eblg.getEBISL()[dit[mybox]];
        const EBGraph& ebgraph = ebisBox.getEBGraph();

        IntVectSet irregIVS = ebisBox.getIrregIVS(curBox);
        IntVectSet multiIVS = ebisBox.getMultiCells(curBox);

        BaseIVFAB<VoFStencil> opStencil(irregIVS,ebgraph, 1);

        //cache the vofIterators
        m_alphaDiagWeight[dit[mybox]].define(irregIVS,ebisBox.getEBGraph(), 1);
        m_vofIterIrreg   [dit[mybox]].define(irregIVS,ebisBox.getEBGraph());
        m_vofIterMulti   [dit[mybox]].define(multiIVS,ebisBox.getEBGraph());

        for (int idir = 0; idir < SpaceDim; idir++)
          {
            IntVectSet loIrreg = irregIVS;
            IntVectSet hiIrreg = irregIVS;
            loIrreg &= sideBoxLo[idir];
            hiIrreg &= sideBoxHi[idir];
            m_vofIterDomLo[idir][dit[mybox]].define(loIrreg,ebisBox.getEBGraph());
            m_vofIterDomHi[idir][dit[mybox]].define(hiIrreg,ebisBox.getEBGraph());
          }

        VoFIterator& vofit = m_vofIterIrreg[dit[mybox]];
        for (vofit.reset(); vofit.ok(); ++vofit)
          {
            const VolIndex& VoF = vofit();

            VoFStencil& curStencil = opStencil(VoF,0);

            //bcoef is included here in the flux consistent
            //with the regular

            getDivFStencil(curStencil,VoF, dit[mybox]);
            if (fluxStencil != NULL)
              {
                BaseIVFAB<VoFStencil>& fluxStencilBaseIVFAB = (*fluxStencil)[dit[mybox]];
                //this fills the stencil with the gradient*beta*bcoef
                VoFStencil  fluxStencilPt = fluxStencilBaseIVFAB(VoF,0);
                curStencil += fluxStencilPt;
              }
	  }

	m_alphaDiagWeight[dit[mybox]].setVal(0,0);

	//Operator ebstencil
        m_opEBStencil[dit[mybox]] = RefCountedPtr<EBStencil>
          (new EBStencil(m_vofIterIrreg[dit[mybox]].getVector(), opStencil, m_eblg.getDBL().get(dit[mybox]),
                         m_eblg.getEBISL()[dit[mybox]], m_ghostCellsPhi, m_ghostCellsRHS, 0, true));
      }//dit
  } //pragma

  if (m_hasFine)
    {

      int ncomp = 1;
      m_fastFR.define(m_eblgFine, m_eblg, m_refToFine, ncomp, false);//s_forceNoEBCF); set to false as default
      m_hasEBCF = m_fastFR.hasEBCF();
    }
  defineEBCFStencils();
}


void
EBDarcyOp::
defineEBCFStencils()
{
  if (m_hasFine && m_hasEBCF)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              //coarse fine stuff is between me and next finer level
              //fine stuff lives over m_eblgfine
              //coar stuff lives over m_eblg
              int index = m_fastFR.index(idir, sit());
              m_stencilCoar[index].define(m_eblg.getDBL());
              m_faceitCoar [index].define(m_eblg.getDBL());

              DataIterator dit = m_eblg.getDBL().dataIterator(); 
              int nbox = dit.size();

//#pragma omp parallel for
              for(int mybox=0; mybox<nbox; mybox++)
                {
                  Vector<FaceIndex>& facesEBCFCoar =  m_faceitCoar[index][dit[mybox]];
                  Vector<VoFStencil>& stencEBCFCoar= m_stencilCoar[index][dit[mybox]];
                  Vector<VoFIterator>& vofitlist = m_fastFR.getVoFItCoar(dit[mybox], idir, sit());
                  //first build up the list of the faces
                  for (int ivofit = 0; ivofit < vofitlist.size(); ivofit++)
                    {
                      VoFIterator& vofit = vofitlist[ivofit];
                      for (vofit.reset(); vofit.ok(); ++vofit)
                        {
                          //on the coarse side of the CF interface we are
                          //looking in the flip direction because we look
                          //back at the interface
                          Vector<FaceIndex> facespt = m_eblg.getEBISL()[dit[mybox]].getFaces(vofit(), idir, flip(sit()));
                          facesEBCFCoar.append(facespt);
                        }
                    }

                  stencEBCFCoar.resize(facesEBCFCoar.size());
                  for (int iface = 0; iface < stencEBCFCoar.size(); iface++)
                    {
                      IntVectSet cfivs; //does not apply here
                      getFluxStencil(stencEBCFCoar[iface], facesEBCFCoar[iface], dit[mybox], idir);
                    }
                }
            }
        }
    }
}

void
EBDarcyOp::
applyOp(LevelData<EBCellFAB>&             a_opPhi,
        const LevelData<EBCellFAB>&       a_phi,
        bool                              a_homogeneousPhysBC)
{
  //homogeneous CFBCs because that is all we can do.
  applyOp(a_opPhi, a_phi, NULL, a_homogeneousPhysBC, true);
}

void EBDarcyOp::
applyAMROp(LevelData<EBCellFAB>&       a_LofPhi,
	   const LevelData<EBCellFAB>& a_phiFine,
	   const LevelData<EBCellFAB>& a_phi,
	   const LevelData<EBCellFAB>& a_phiCoar,
	   const bool                  a_homogeneousPhysBC,
	   EBDarcyOp*                  a_finerOp)
{
  CH_TIMERS("EBDarcyOp::AMROperator");
  CH_TIMER("applyOp", t1);
  CH_TIMER("reflux", t2);
  CH_assert(a_LofPhi.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_LofPhi.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);

  //apply the operator between this and the next coarser level.
  CH_START(t1);
  applyOp(a_LofPhi, a_phi, &a_phiCoar,  a_homogeneousPhysBC, false);
  CH_STOP(t1);

  //  dumpLevelPoint(a_LofPhi, string("EBConductivityOp: AMROperator: pre-reflux lphi = "));
  //now reflux to enforce flux-matching from finer levels
  if (m_hasFine)
    {
      CH_assert(a_finerOp != NULL);
      CH_START(t2);

      reflux(a_LofPhi, a_phiFine, a_phi, a_finerOp);

      CH_STOP(t2);
    }
  //  dumpLevelPoint(a_LofPhi, string("EBConductivityOp: AMROperator: post-reflux lphi = "));
}
void
EBDarcyOp::reflux(LevelData<EBCellFAB>& a_residual,
		  const LevelData<EBCellFAB>& a_phiFine,
		  const LevelData<EBCellFAB>& a_phi,
		  EBDarcyOp*                  a_finerOp)
{
  CH_TIMERS("EBDarcyOp::fastReflux");
  CH_TIMER("setToZero",t2);
  CH_TIMER("incrementCoar",t3);
  CH_TIMER("incrementFine",t4);
  CH_TIMER("reflux_from_reg",t5);
  Interval interv(0,0);

  CH_START(t2);
  m_fastFR.setToZero();
  CH_STOP(t2);
  CH_START(t3);
  incrementFRCoar(m_fastFR, a_phiFine, a_phi);
  CH_STOP(t3);

  CH_START(t4);
  incrementFRFine(m_fastFR, a_phiFine, a_phi, a_finerOp);
  CH_STOP(t4);
  CH_START(t5);
  //make this 1/dx*dy*dz
  const Real cellvol = m_vectDx[0]*m_vectDx[1]*m_vectDx[2];
      Real scale = 1.0/cellvol;
  m_fastFR.reflux(a_residual, interv, scale);

  CH_STOP(t5);
}

void 
EBDarcyOp::incrementFRCoar(EBFastFR&             a_fluxReg,
			   const LevelData<EBCellFAB>& a_phiFine,
			   const LevelData<EBCellFAB>& a_phi)
{
  CH_TIME("EBDarcyOp::incrementFRCoar");
  CH_assert(a_phiFine.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);

  int ncomp = 1;
  Interval interv(0,0);

  DataIterator dit = m_eblg.getDBL().dataIterator(); 
  int nbox = dit.size();
//#pragma omp parallel for
  for(int mybox=0; mybox<nbox; mybox++)
    {

      const EBCellFAB& coarfab = a_phi[dit[mybox]];
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit[mybox]];
      const Box&  box = m_eblg.getDBL().get(dit[mybox]);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //no boundary faces here.

          Box ghostedBox = box;
          ghostedBox.grow(1);
          ghostedBox.grow(idir,-1);
          ghostedBox &= m_eblg.getDomain();

          EBFaceFAB coarflux(ebisBox, ghostedBox, idir, ncomp);

          //old way
          //getFlux(coarflux, coarfab, ghostedBox, box, m_eblg.getDomain(), ebisBox, m_dx, dit[mybox], idir);

          // new way
          getFluxRegOnly(coarflux, coarfab, ghostedBox, m_vectDx, dit[mybox], idir);
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Vector<FaceIndex>*  faceit=NULL;
              Vector<VoFStencil>* stencil=NULL;
              int index = EBFastFR::index(idir, sit());
              if (m_hasEBCF)
                {
                  faceit  = &( m_faceitCoar[index][dit[mybox]]);
                  stencil = &(m_stencilCoar[index][dit[mybox]]);
                }
              getFluxEBCF(coarflux, coarfab, ghostedBox, *faceit, *stencil);
            }

          //          dumpFlux(coarflux, idir,  string("incrementFRCoar: flux = "));
	  //make this dx*dy (exclude the dx_idir)
          Real scale = 1.0; //beta and bcoef already in flux
	  for (int iidir = 0; iidir < SpaceDim; iidir++)
	  {
	      if (iidir != idir) scale *=m_vectDx[iidir];
	  }
	  for (SideIterator sit; sit.ok(); ++sit)
	  {
              a_fluxReg.incrementCoarseBoth(coarflux, scale, dit[mybox], interv, idir, sit());
	  }
        }
    }
}

void
EBDarcyOp::incrementFRFine(EBFastFR&             a_fluxReg,
			   const LevelData<EBCellFAB>& a_phiFine,
			   const LevelData<EBCellFAB>& a_phi,
			   EBDarcyOp*                  a_finerOp)
{
  CH_TIME("EBConductivityOp::incrementFRFine");
  CH_assert(a_phiFine.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);
  CH_assert(m_hasFine);
  int ncomp = 1;
  Interval interv(0,0);
  EBDarcyOp& finerEBAMROp = *a_finerOp;

  //ghost cells of phiFine need to be filled
  LevelData<EBCellFAB>& phiFine = (LevelData<EBCellFAB>&) a_phiFine;
  finerEBAMROp.m_quadCFIWithCoar->interpolate(phiFine, a_phi, interv);
  phiFine.exchange(interv);
  
  DataIterator ditf = a_phiFine.dataIterator();
  int nbox = ditf.size();
//#pragma omp parallel for
 for(int mybox=0; mybox<nbox; mybox++)
   {
      const Box&     boxFine = m_eblgFine.getDBL().get(ditf[mybox]);
      const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[ditf[mybox]];
      const EBCellFAB& phiFine = a_phiFine[ditf[mybox]];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); sit.next())
            {
              Box fabBox = adjCellBox(boxFine, idir, sit(), 1);
              fabBox.shift(idir, -sign(sit()));

              Box ghostedBox = fabBox;
              ghostedBox.grow(1);
              ghostedBox.grow(idir,-1);
              ghostedBox &= m_eblgFine.getDomain();

              EBFaceFAB fluxFine(ebisBoxFine, ghostedBox, idir, ncomp);
              finerEBAMROp.getFlux(fluxFine, phiFine, ghostedBox, fabBox, m_eblgFine.getDomain(),
                                   ebisBoxFine, m_vectDxFine, ditf[mybox], idir);

	      //make this dx*dy (exclude the dx_idir)
              Real scale = 1.0; //beta and bcoef already in flux
	      for (int iidir = 0; iidir < SpaceDim; iidir++)
	      {
		  if (iidir != idir) scale *=m_vectDx[iidir];
	      }

              a_fluxReg.incrementFineBoth(fluxFine, scale, ditf[mybox], interv, idir, sit());
            }
        }
    }
}

//-----------------------------------------------------------------------
/***/
void
EBDarcyOp::
incrOpRegularAllDirs(Box * a_loBox,
                     Box * a_hiBox,
                     int * a_hasLo,
                     int * a_hasHi,
                     Box & a_curDblBox,
                     Box & a_curPhiBox,
                     int a_nComps,
                     BaseFab<Real> & a_curOpPhiFAB,
                     const BaseFab<Real> & a_curPhiFAB,
                     bool a_homogeneousPhysBC,
                     const DataIndex& a_dit)
{
  CH_TIME("EBDarcyOp::incrOpRegularAllDirs");
  CH_assert(m_domainBC != NULL);

  //need to monkey with the ghost cells to account for boundary conditions
  if (!s_turnOffBCs)
  {
      BaseFab<Real>& phiFAB = (BaseFab<Real>&) a_curPhiFAB;
      applyDomainFlux(a_loBox, a_hiBox, a_hasLo, a_hasHi,
                      a_curDblBox, a_nComps, phiFAB,
                      a_homogeneousPhysBC, a_dit);
  }

  for (int comp = 0; comp<a_nComps; comp++)
    {
      //data ptr fusses if it is truly zero size
      BaseFab<Real> dummy(Box(IntVect::Zero, IntVect::Zero), 1);

      BaseFab<Real>* bc[3];
      BaseFab<Real>* darcyflux[3];
      //need three coeffs because this has to work in 3d
      for (int iloc = 0; iloc < 3; iloc++)
        {
 	    bc[iloc] = &((*m_bcoef)[a_dit][iloc].getSingleValuedFAB());
	    darcyflux[iloc] = &((*m_DarcyFlux)[a_dit][iloc].getSingleValuedFAB());
	}

      FORT_DARCYCONDUCTIVITYINPLACE(CHF_FRA1(a_curOpPhiFAB,comp),
				    CHF_FRA1((*darcyflux[0]), comp),
				    CHF_FRA1((*darcyflux[1]), comp),
				    CHF_FRA1((*darcyflux[2]), comp),
				    CHF_CONST_FRA1(a_curPhiFAB,comp),
				    CHF_CONST_FRA1((*bc[0]),comp),
				    CHF_CONST_FRA1((*bc[1]),comp),
				    CHF_CONST_FRA1((*bc[2]),comp),
				    CHF_CONST_REAL(m_beta),
				    CHF_CONST_REALVECT(m_vectDx),
				    CHF_BOX(a_curDblBox));

    }
}

/***/
void
EBDarcyOp::
applyDomainFlux(Box * a_loBox,
                Box * a_hiBox,
                int * a_hasLo,
                int * a_hasHi,
                Box & a_dblBox,
                int a_nComps,
                BaseFab<Real> & a_phiFAB,
                bool a_homogeneousPhysBC,
                const DataIndex& a_dit)
{
  CH_TIME("EBDarcyOp::applyDomainFlux");
  CH_assert(m_domainBC != NULL);

  for (int idir=0; idir<SpaceDim; idir++)
  {

      EBArith::loHi(a_loBox[idir], a_hasLo[idir],
                    a_hiBox[idir], a_hasHi[idir],
                    m_eblg.getDomain(),a_dblBox, idir);

      for (int comp = 0; comp<a_nComps; comp++)
      {

          if (a_hasLo[idir] == 1)
	  {
              Box lbox=a_loBox[idir];
              lbox.shift(idir,-1);
              FArrayBox loFaceFlux(a_loBox[idir],a_nComps);
              int side = -1;
              m_domainBC->getFaceFlux(loFaceFlux,a_phiFAB,RealVect::Zero,m_vectDx, idir,Side::Lo,a_dit, m_time,a_homogeneousPhysBC);

              BaseFab<Real>& bc = ((*m_bcoef)[a_dit][idir].getSingleValuedFAB());
              //again, following the odd convention of EBAMRPoissonOp
              //(because I am reusing its BC classes),
              //the input flux here is CELL centered and the input box
              //is the box adjacent to the domain boundary on the valid side.
              //because I am not insane (yet) I will just shift the flux's box
              //over and multiply by the appropriate coefficient
              bc.shiftHalf(idir, 1);
              FORT_DARCYEBCOREGAPPLYDOMAINFLUX(CHF_FRA1(a_phiFAB,comp),
                                          CHF_CONST_FRA1(loFaceFlux,comp),
                                          CHF_CONST_FRA1(bc,comp),
                                          CHF_CONST_REAL(m_vectDx[idir]),
                                          CHF_CONST_INT(side),
                                          CHF_CONST_INT(idir),
                                          CHF_BOX(lbox));
              bc.shiftHalf(idir, -1);

               BaseFab<Real>& darcyflux = (*m_DarcyFlux)[a_dit][idir].getSingleValuedFAB();
	       darcyflux.shiftHalf(idir, 1);
	       darcyflux.copy(loFaceFlux);
	       darcyflux.shiftHalf(idir, -1);
	  }

          if (a_hasHi[idir] == 1)
	  {
              Box hbox=a_hiBox[idir];
              hbox.shift(idir,1);
              FArrayBox hiFaceFlux(a_hiBox[idir],a_nComps);
              int side = 1;
              m_domainBC->getFaceFlux(hiFaceFlux,a_phiFAB,RealVect::Zero,m_vectDx,idir,Side::Hi,a_dit,m_time,a_homogeneousPhysBC);

              BaseFab<Real>& bc = ((*m_bcoef)[a_dit][idir].getSingleValuedFAB());
              //again, following the odd convention of EBAMRPoissonOp
              //(because I am reusing its BC classes),
              //the input flux here is CELL centered and the input box
              //is the box adjacent to the domain boundary on the valid side.
              //because I am not insane (yet) I will just shift the flux's box
              //over and multiply by the appropriate coefficient
              bc.shiftHalf(idir, -1);
              FORT_DARCYEBCOREGAPPLYDOMAINFLUX(CHF_FRA1(a_phiFAB,comp),
                                          CHF_CONST_FRA1(hiFaceFlux,comp),
                                          CHF_CONST_FRA1(bc,comp),
                                          CHF_CONST_REAL(m_vectDx[idir]),
                                          CHF_CONST_INT(side),
                                          CHF_CONST_INT(idir),
                                          CHF_BOX(hbox));
              bc.shiftHalf(idir, 1);

               BaseFab<Real>& darcyflux = (*m_DarcyFlux)[a_dit][idir].getSingleValuedFAB();
	       darcyflux.shiftHalf(idir, -1);
	       darcyflux.copy(hiFaceFlux);
	       darcyflux.shiftHalf(idir, 1);

	  }
	  
      }
  }
}

//-----------------------------------------------------------------------
void
EBDarcyOp::
applyOp(LevelData<EBCellFAB>&                    a_lhs,
        const LevelData<EBCellFAB>&              a_phi,
        const LevelData<EBCellFAB>* const        a_phiCoar,
        const bool&                              a_homogeneousPhysBC,
        const bool&                              a_homogeneousCFBC)
{
  CH_TIME("ebco::applyOp");
  //we do not need it since it's interpolated in Newton-Krylov function

  // LevelData<EBCellFAB>& phi = const_cast<LevelData<EBCellFAB>&>(a_phi);
  // if (m_hasCoar && (!s_turnOffBCs))
  //   {
  //     //applyCFBCs(phi, a_phiCoar, a_homogeneousCFBC);
  //     Interval interv(0,0);
  //     m_quadCFIWithCoar->interpolate(phi, *a_phiCoar, interv);
  //   }
  // phi.exchange(phi.interval());

  EBLevelDataOps::setToZero(a_lhs);
  //incr( a_lhs, a_phi, m_alpha); //this multiplies by alpha
  DataIterator dit = m_eblg.getDBL().dataIterator(); 
  int nbox = dit.size();
//#pragma omp parallel for
  for(int mybox=0; mybox<nbox; mybox++)
    {

	//a_lhs[dit[mybox]].mult((*m_acoef)[dit[mybox]], 0, 0, 1);

      Box loBox[SpaceDim],hiBox[SpaceDim];
      int hasLo[SpaceDim],hasHi[SpaceDim];
      EBCellFAB      & phi = (EBCellFAB&)(a_phi[dit[mybox]]);
      EBCellFAB      & lph = a_lhs[dit[mybox]];
      //phi.setCoveredCellVal(0.0, 0);

      const BaseFab<Real>  & phiFAB = phi.getSingleValuedFAB();
      BaseFab<Real>        & lphFAB = lph.getSingleValuedFAB();
      Box dblBox = m_eblg.getDBL()[dit[mybox]];
      int nComps = 1;
      Box curPhiBox = phiFAB.box();

      if (!s_turnOffBCs)
        {
          incrOpRegularAllDirs( loBox, hiBox, hasLo, hasHi,
                                dblBox, curPhiBox, nComps,
                                lphFAB,
                                phiFAB,
                                a_homogeneousPhysBC,
                                dit[mybox]);
        }
      else
        {
          //the all dirs code is wrong for no bcs = true
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              incrOpRegularDir(a_lhs[dit[mybox]], a_phi[dit[mybox]], a_homogeneousPhysBC, idir, dit[mybox]);
            }
        }

      applyOpIrregular(a_lhs[dit[mybox]], a_phi[dit[mybox]], a_homogeneousPhysBC, dit[mybox]);

    }
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
incrOpRegularDir(EBCellFAB&             a_lhs,
                 const EBCellFAB&       a_phi,
                 const bool&            a_homogeneous,
                 const int&             a_dir,
                 const DataIndex&       a_datInd)
{
  CH_TIME("ebco::incrOpReg");
  const Box& grid = m_eblg.getDBL()[a_datInd];
  Box domainFaces = m_eblg.getDomain().domainBox();
  domainFaces.surroundingNodes(a_dir);
  Box interiorFaces = grid;
  interiorFaces.surroundingNodes(a_dir);
  interiorFaces.grow(a_dir, 1);
  interiorFaces &=  domainFaces;
  interiorFaces.grow( a_dir, -1);

  //do flux difference for interior points
  FArrayBox interiorFlux(interiorFaces, 1);
  const FArrayBox& phi  = (FArrayBox&)(a_phi.getSingleValuedFAB());
  getFlux(interiorFlux, phi,  interiorFaces, a_dir, m_vectDx, a_datInd);

  Box loBox, hiBox, centerBox;
  int hasLo, hasHi;
  EBArith::loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, m_eblg.getDomain(),grid, a_dir);

  //do the low high center thing
  BaseFab<Real>& reglhs        = a_lhs.getSingleValuedFAB();
  Box dummyBox(IntVect::Zero, IntVect::Unit);
  FArrayBox domainFluxLo(dummyBox, 1);
  FArrayBox domainFluxHi(dummyBox, 1);

  RealVect origin = RealVect::Zero;
  if (hasLo==1)
    {
      Box loBoxFace = loBox;
      loBoxFace.shiftHalf(a_dir, -1);
      domainFluxLo.resize(loBoxFace, 1);
      if (!s_turnOffBCs)
        {
          //using EBAMRPoissonOp BCs which require a cell centered box.
          domainFluxLo.shiftHalf(a_dir, 1);
          m_domainBC->getFaceFlux(domainFluxLo,phi,origin,m_vectDx,a_dir,Side::Lo,a_datInd,m_time,a_homogeneous);
          domainFluxLo *= m_beta;
          domainFluxLo.shiftHalf(a_dir,-1);
        }
      else
        {
          //extrapolate to domain flux if there is no bc
          for (BoxIterator boxit(loBoxFace); boxit.ok(); ++boxit)
            {
              const IntVect& iv = boxit();
              IntVect ivn= iv;
              ivn[a_dir]++;
              domainFluxLo(iv, 0) = interiorFlux(ivn, 0);
            }
        }
    }
  if (hasHi==1)
    {
      Box hiBoxFace = hiBox;
      hiBoxFace.shiftHalf(a_dir, 1);
      domainFluxHi.resize(hiBoxFace, 1);
      if (!s_turnOffBCs)
        {
          //using EBAMRPoissonOp BCs which require a cell centered box.
          domainFluxHi.shiftHalf(a_dir, -1);
          m_domainBC->getFaceFlux(domainFluxHi,phi,origin,m_vectDx,a_dir,Side::Hi,a_datInd,m_time,a_homogeneous);
          domainFluxHi *= m_beta;
          domainFluxHi.shiftHalf(a_dir,  1);
        }
      else
        {
          //extrapolate to domain flux if there is no bc
          for (BoxIterator boxit(hiBoxFace); boxit.ok(); ++boxit)
            {
              const IntVect& iv = boxit();
              IntVect ivn= iv;
              ivn[a_dir]--;
              domainFluxHi(iv, 0) = interiorFlux(ivn, 0);
            }
        }
    }
  Real unity = 1.0;  //beta already in the flux
  FORT_DARCYINCRAPPLYEBCO(CHF_FRA1(reglhs,0),
                     CHF_CONST_FRA1(interiorFlux, 0),
                     CHF_CONST_FRA1(domainFluxLo, 0),
                     CHF_CONST_FRA1(domainFluxHi, 0),
                     CHF_CONST_REAL(unity),
                     CHF_CONST_REAL(m_vectDx[a_dir]),
                     CHF_BOX(loBox),
                     CHF_BOX(hiBox),
                     CHF_BOX(centerBox),
                     CHF_CONST_INT(hasLo),
                     CHF_CONST_INT(hasHi),
                     CHF_CONST_INT(a_dir));
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
applyOpIrregular(EBCellFAB&             a_lhs,
                 const EBCellFAB&       a_phi,
                 const bool&            a_homogeneous,
                 const DataIndex&       a_datInd)
{
  CH_TIME("ebco::applyOpIrr");


  m_opEBStencil[a_datInd]->apply(a_lhs, a_phi, m_alphaDiagWeight[a_datInd], /*m_alpha*/ 0.0, m_beta, false);


  if (!a_homogeneous)
  {
    
    Real factor;
    factor = m_beta;//beta and bcoef handled within applyEBFlux
     
      /*
	NeumannConductivityEBBC* castbc = dynamic_cast<NeumannConductivityEBBC*> (&(*m_ebBC));
	if(castbc == NULL) MayDay::Error("cast failed");
	castbc->applyEBFluxVectDx(...);
      */
      m_ebBC->applyEBFlux(a_lhs, a_phi, m_vofIterIrreg[a_datInd], (*(m_eblg.getCFIVS())),
                          a_datInd, RealVect::Zero, m_vectDx, factor,
                          a_homogeneous, m_time);
    }
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int comp = 0;
      for (m_vofIterDomLo[idir][a_datInd].reset(); m_vofIterDomLo[idir][a_datInd].ok();  ++m_vofIterDomLo[idir][a_datInd])
        {
          Real flux;
          const VolIndex& vof = m_vofIterDomLo[idir][a_datInd]();
          m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                  RealVect::Zero, m_vectDx,idir,Side::Lo, a_datInd, m_time,
                                  a_homogeneous);
          //area gets multiplied in by bc operator
          a_lhs(vof,comp) -= flux*m_beta/m_vectDx[idir];
        }
      for (m_vofIterDomHi[idir][a_datInd].reset(); m_vofIterDomHi[idir][a_datInd].ok();  ++m_vofIterDomHi[idir][a_datInd])
        {
          Real flux;
          const VolIndex& vof = m_vofIterDomHi[idir][a_datInd]();
          m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                  RealVect::Zero, m_vectDx,idir,Side::Hi,a_datInd,0.0,
                                  a_homogeneous);
          //area gets multiplied in by bc operator
          a_lhs(vof,comp) += flux*m_beta/m_vectDx[idir];
        }
    }
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
applyOpNoBoundary(LevelData<EBCellFAB>&        a_opPhi,
                  const LevelData<EBCellFAB>&  a_phi)
{
  s_turnOffBCs = true;
  applyOp(a_opPhi, a_phi, true);
  s_turnOffBCs = false;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBDarcyOp::
create(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("ebco::create");
  int ncomp = a_rhs.nComp();
  EBCellFactory ebcellfact(m_eblg.getEBISL());
  a_lhs.define(m_eblg.getDBL(), ncomp, a_rhs.ghostVect(), ebcellfact);
}

void
EBDarcyOp::
assign(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("ebco::assign");
  EBLevelDataOps::assign(a_lhs,a_rhs);
}
//-----------------------------------------------------------------------
Real
EBDarcyOp::
dotProduct(const LevelData<EBCellFAB>& a_1,
           const LevelData<EBCellFAB>& a_2)
{
  CH_TIME("ebco::dotProd");
  ProblemDomain domain;
  Real volume;

  return EBLevelDataOps::kappaDotProduct(volume,a_1,a_2,EBLEVELDATAOPS_ALLVOFS,domain);
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
incr(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     Real                        a_scale)
{
  CH_TIME("ebco::incr");
  EBLevelDataOps::incr(a_lhs,a_x,a_scale);
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
axby(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     const LevelData<EBCellFAB>& a_y,
     Real                        a_a,
     Real                        a_b)
{
  CH_TIME("ebco::axby");
  EBLevelDataOps::axby(a_lhs,a_x,a_y,a_a,a_b);
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
scale(LevelData<EBCellFAB>& a_lhs,
      const Real&           a_scale)
{
  CH_TIME("ebco::scale");
  EBLevelDataOps::scale(a_lhs,a_scale);
}
//-------------------------------
Real 
EBDarcyOp::
norm(const LevelData<EBCellFAB>& a_rhs,
     int                         a_ord)
{
 CH_TIMERS("EBDarcyOp::norm");
 CH_TIMER("mpi_allreduce",t1);

 Real maxNorm = 0.0;

 maxNorm = localMaxNorm(a_rhs);

 CH_START(t1);
#ifdef CH_MPI
       Real tmp = 1.;
       int result = MPI_Allreduce(&maxNorm, &tmp, 1, MPI_CH_REAL,
                         MPI_MAX, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
         { //bark!!!
           MayDay::Error("sorry, but I had a communcation error on norm");
         }
       maxNorm = tmp;
#endif
 CH_STOP(t1);

 return maxNorm;
}

Real 
EBDarcyOp::
localMaxNorm(const LevelData<EBCellFAB>& a_rhs)
{
 CH_TIME("EBAMRPoissonOp::localMaxNorm");
 return  EBAMRPoissonOp::staticMaxNorm(a_rhs, m_eblg);
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
setToZero(LevelData<EBCellFAB>& a_lhs)
{
  CH_TIME("ebco::setToZero");
  EBLevelDataOps::setToZero(a_lhs);
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
setVal(LevelData<EBCellFAB>& a_lhs, const Real& a_value)
{
  CH_TIME("ebco::setVal");
  EBLevelDataOps::setVal(a_lhs, a_value);
}

void
EBDarcyOp::
getFlux(EBFluxFAB&                    a_flux,
        const LevelData<EBCellFAB>&   a_data,
        const Box&                    a_grid,
        const DataIndex&              a_dit,
        Real                          a_scale)
{
  CH_TIME("ebco::getflux1");
  a_flux.define(m_eblg.getEBISL()[a_dit], a_grid, 1);
  a_flux.setVal(0.);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box ghostedBox = a_grid;
      ghostedBox.grow(1);
      ghostedBox.grow(idir,-1);
      ghostedBox &= m_eblg.getDomain();


      getFlux(a_flux[idir], a_data[a_dit], ghostedBox, a_grid,
              m_eblg.getDomain(), m_eblg.getEBISL()[a_dit], m_vectDx, a_dit, idir);

    }
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
getFlux(EBFaceFAB&                    a_fluxCentroid,
        const EBCellFAB&              a_phi,
        const Box&                    a_ghostedBox,
        const Box&                    a_fabBox,
        const ProblemDomain&          a_domain,
        const EBISBox&                a_ebisBox,
        const RealVect&               a_vectDx,
        const DataIndex&              a_datInd,
        const int&                    a_idir)
{
  CH_TIME("ebco::getFlux2");
  //has some extra cells so...
  a_fluxCentroid.setVal(0.);
  int ncomp = a_phi.nComp();
  CH_assert(ncomp == a_fluxCentroid.nComp());
  Box cellBox = a_ghostedBox;
  //want only interior faces
  cellBox.grow(a_idir, 1);
  cellBox &= a_domain;
  cellBox.grow(a_idir,-1);

  Box faceBox = surroundingNodes(cellBox, a_idir);
  EBFaceFAB fluxCenter(a_ebisBox, a_ghostedBox, a_idir,1);

  //make a EBFaceFAB (including ghost cells) that will hold centered gradients
  BaseFab<Real>& regFlux  = fluxCenter.getSingleValuedFAB();
  const BaseFab<Real>& regPhi = a_phi.getSingleValuedFAB();
  const EBFaceFAB& bcoebff = (*m_bcoef)[a_datInd][a_idir];
  const FArrayBox& regBCo = (const FArrayBox&)(bcoebff.getSingleValuedFAB());

  FORT_DARCYGETFLUXEBCO(CHF_FRA1(regFlux,0),
                   CHF_CONST_FRA1(regBCo,0),
                   CHF_CONST_FRA1(regPhi, 0),
                   CHF_BOX(faceBox),
                   CHF_CONST_REAL(a_vectDx[a_idir]),
                   CHF_CONST_INT(a_idir));


  a_fluxCentroid.copy(fluxCenter);

  IntVectSet ivsCell = a_ebisBox.getIrregIVS(cellBox);
  if (!ivsCell.isEmpty())
    {
      FaceStop::WhichFaces stopCrit = FaceStop::SurroundingNoBoundary;

      for (FaceIterator faceit(ivsCell, a_ebisBox.getEBGraph(), a_idir,stopCrit);
          faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          Real phiHi = a_phi(face.getVoF(Side::Hi), 0);
          Real phiLo = a_phi(face.getVoF(Side::Lo), 0);
          Real fluxFace = bcoebff(face, 0)*(phiHi - phiLo)/a_vectDx[a_idir];
          fluxCenter(face, 0) = fluxFace;
        }
      //interpolate from face centers to face centroids
      Box cellBox = a_fluxCentroid.getCellRegion();
      EBArith::interpolateFluxToCentroids(a_fluxCentroid,
                                          fluxCenter,
                                          a_fabBox,
                                          a_ebisBox,
                                          a_domain,
                                          a_idir);
    }

  a_fluxCentroid *= m_beta;
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
getFluxRegOnly(EBFaceFAB&                    a_fluxCentroid,
               const EBCellFAB&              a_phi,
               const Box&                    a_ghostedBox,
               const RealVect&               a_vectDx,
               const DataIndex&              a_datInd,
               const int&                    a_idir)
{
  CH_TIME("ebco::getFluxRegOnly");
  const ProblemDomain& domain = m_eblg.getDomain();

  //has some extra cells so...
  a_fluxCentroid.setVal(0.);
  int ncomp = a_phi.nComp();
  CH_assert(ncomp == a_fluxCentroid.nComp());
  Box cellBox = a_ghostedBox;
  //want only interior faces
  cellBox.grow(a_idir, 1);
  cellBox &= domain;
  cellBox.grow(a_idir,-1);

  Box faceBox = surroundingNodes(cellBox, a_idir);

  //make a EBFaceFAB (including ghost cells) that will hold centered gradients
  BaseFab<Real>& regFlux = a_fluxCentroid.getSingleValuedFAB();
  const BaseFab<Real>& regPhi = a_phi.getSingleValuedFAB();
  const EBFaceFAB& bcoebff = (*m_bcoef)[a_datInd][a_idir];
  const FArrayBox& regBCo = (const FArrayBox&)(bcoebff.getSingleValuedFAB());

  FORT_DARCYGETFLUXEBCO(CHF_FRA1(regFlux,0),
                   CHF_CONST_FRA1(regBCo,0),
                   CHF_CONST_FRA1(regPhi, 0),
                   CHF_BOX(faceBox),
                   CHF_CONST_REAL(a_vectDx[a_idir]),
                   CHF_CONST_INT(a_idir));

  a_fluxCentroid *= m_beta;
}
//-----------------------------------------------------------------------
void
EBDarcyOp::
getFlux(FArrayBox&                    a_flux,
        const FArrayBox&              a_phi,
        const Box&                    a_faceBox,
        const int&                    a_idir,
        const RealVect&               a_vectDx,
        const DataIndex&              a_datInd)
{
  CH_TIME("ebco::getflux3");
  const EBFaceFAB& bcoebff = (*m_bcoef)[a_datInd][a_idir];
  const FArrayBox& regBCo = (const FArrayBox&)(bcoebff.getSingleValuedFAB());
  FORT_DARCYGETFLUXEBCO(CHF_FRA1(a_flux,0),
                   CHF_CONST_FRA1(regBCo,0),
                   CHF_CONST_FRA1(a_phi, 0),
                   CHF_BOX(a_faceBox),
                   CHF_CONST_REAL(a_vectDx[a_idir]),
                   CHF_CONST_INT(a_idir));

  a_flux *= m_beta;
}

void
EBDarcyOp::
getFluxIrregular(EBCellFAB&                   a_lhs,
                 const LevelData<EBCellFAB>&  a_phi,
                 const DataIndex&             a_datInd)
{
  CH_TIME("ebco::applyOpIrr");
  const EBISBox& ebisBox = m_eblg.getEBISL()[a_datInd];

  VoFIterator& vofit = m_vofIterIrreg[a_datInd];
  for (vofit.reset(); vofit.ok(); ++vofit)
  {
    Real flux;
    const VolIndex& vof = vofit();
    if (ebisBox.bndryArea(vof)>1.e-13)
    {
	const RealVect&   normal = ebisBox.normal(vof);
	LayoutData<IntVectSet> dummy_ivs;
	m_ebBC->getEBFlux(flux,vof,a_phi, dummy_ivs, a_datInd, RealVect::Zero, m_vectDx, false, m_time);
	
	Real norm_denom=0.0;
	for (int idir = 0; idir < SpaceDim; idir++)
	{
	    norm_denom += normal[idir]*normal[idir] 
		/(m_vectDx[idir]*m_vectDx[idir]);
	}
	for (int idir = 0; idir < SpaceDim; idir++)
	{
	    a_lhs(vof,idir) +=flux*normal[idir]/m_vectDx[idir]/sqrt(norm_denom);
	}
    }
  }
}
void
EBDarcyOp::
getFluxEBCF(EBFaceFAB&                    a_flux,
            const EBCellFAB&              a_phi,
            const Box&                    a_ghostedBox,
            Vector<FaceIndex>&            a_faceitEBCF,
            Vector<VoFStencil>&           a_stenEBCF)
{
  CH_TIME("EBDarcyOp::getFluxEBCF");

  //only do the evil stuff if you have a coarse-fine /  EB crossing situation

  if (m_hasEBCF)
    {
      CH_TIME("EBCF stuff");
      for (int iface = 0; iface < a_faceitEBCF.size(); iface++)
        {
          const FaceIndex& face =     a_faceitEBCF[iface];
          const VoFStencil& stencil   = a_stenEBCF[iface];
          Real fluxval = 0;
          for (int isten = 0; isten < stencil.size(); isten++)
            {
              fluxval += stencil.weight(isten)*(a_phi(stencil.vof(isten), 0));
            }
          //note the last minute beta
          a_flux(face, 0) = m_beta*fluxval;
        }
    }
}

#include "NamespaceFooter.H"
